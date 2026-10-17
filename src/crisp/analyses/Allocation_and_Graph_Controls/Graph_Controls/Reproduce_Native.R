arguments = commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4L)
work = normalizePath(arguments[[1]], mustWork = TRUE)
base = normalizePath(arguments[[2]], mustWork = TRUE)
prepared_path = normalizePath(arguments[[3]], mustWork = TRUE)
code = normalizePath(arguments[[4]], mustWork = TRUE)
stopifnot(!dir.exists(file.path(work, 'native')), !file.exists(file.path(work, 'NATIVE_EXIT.json')))
native = file.path(code, 'runtime/native')
source(file.path(native, 'functions.R'))
source(file.path(native, 'original_adapter.R'))
require_runtime()
stopifnot(file.exists(file.path(work, 'GRAPH_MANIFEST.json')))

certificate = function(transition, damping, q, values) {
    u = .Machine$double.eps / 2
    gamma = function(count) count * u / (1 - count * u)
    rows = Matrix::rowSums(transition != 0)
    columns = Matrix::colSums(transition != 0)
    propagated = as.numeric(transition %*% values)
    constant = (1 - damping) * q
    residual = damping * propagated + constant - values
    allowance = gamma(rows + 3) / (1 - gamma(rows + 3)) * abs(damping * propagated)
    allowance = allowance + gamma(2) / (1 - gamma(2)) * abs(constant)
    allowance = allowance + gamma(3) / (1 - gamma(3)) * (abs(damping * propagated) + abs(constant) + abs(values))
    roundoff = sum(allowance) / (1 - gamma(length(values)))
    norm_upper = max(Matrix::colSums(Matrix::Diagonal(x = damping) %*% transition)) / (1 - gamma(max(columns) + 3))
    stopifnot(norm_upper < 1, all(is.finite(values)), all(values >= 0), all(damping >= 0), all(damping < 1), all(q >= 0))
    positive = function(x) if (any(x > 0)) min(x[x > 0]) else 1
    minimum_product = min(positive(transition@x) * positive(values) * positive(damping), positive(1 - damping) * positive(q))
    stopifnot(minimum_product > .Machine$double.xmin)
    residual_upper = sum(abs(residual)) / (1 - gamma(length(values))) + roundoff
    list(residual_l1 = sum(abs(residual)), roundoff_l1 = roundoff, operator_l1_upper = norm_upper, error_bound = residual_upper / (1 - norm_upper), minimum_positive_product = minimum_product, maximum_column_nnz = max(columns))
}

fixture_adjacency = Matrix::Matrix(matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), 3, 3), sparse = TRUE)
fixture_transition = dawn_sparse_transition(fixture_adjacency)
fixture_damping = dawn_sparse_damping(fixture_adjacency, 3)
fixture_q = c(0.6, 0.3, 0.1)
fixture_result = dawn_instrumented(fixture_transition, fixture_q, rep(1, 3), fixture_damping, maxit = 10000L, epsilon = 1e-12)
fixture_scores = as.numeric(fixture_result$summaryOutput$Rank)
fixture_exact = solve(diag(3) - diag(fixture_damping) %*% as.matrix(fixture_transition), (1 - fixture_damping) * fixture_q)
stopifnot(max(abs(fixture_scores - fixture_exact)) <= certificate(fixture_transition, fixture_damping, fixture_q, fixture_scores)$error_bound, fixture_scores[3] == 0.1)
directed = Matrix::Matrix(matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 0), 3, 3), sparse = TRUE)
directed_transition = dawn_sparse_transition(directed)
stopifnot(max(abs(Matrix::colSums(directed_transition) - 1)) < 1e-15)
write_json(list(status = 'PASS', fixtures = c('isolate', 'sparse_dense_direct_solve', 'column_orientation')), file.path(work, 'NATIVE_FIXTURES.json'))

prepared = readRDS(prepared_path)
stopifnot(length(prepared$genes) == 7399L, length(prepared$models) == 85L, nrow(prepared$aliases) == 340L)
stopifnot(identical(as.integer(prepared$seeds), c(104729L, 130363L, 155921L, 196613L, 228017L)))
stopifnot(identical(prepared$genes, readLines(file.path(base, 'reference_allocations/numerical/genes.txt'))))
metadata = read.csv(file.path(work, 'graph_metadata.csv'), stringsAsFactors = FALSE)
stopifnot(all(metadata$new_fit %in% c('True', 'False')))
metadata$new_fit = metadata$new_fit == 'True'
stopifnot(nrow(metadata) == 12L, sum(metadata$new_fit) == 10L)
dir.create(file.path(work, 'native'), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(work, 'native/R_session.txt'))
for (h in 0:1) {
    evidence = prepared$evidence[[paste0('H', h)]]$values
    q = sweep(evidence, 2, colSums(evidence), '/')
    writeBin(as.double(q), file.path(work, paste0('native/H', h, '_q.bin')), size = 8L, endian = 'little')
}
write_json(list(status = 'PASS', prepared_path = prepared_path, genes = 7399L, models = 85L, seeds = as.integer(prepared$seeds), labels_in_fitting = FALSE), file.path(work, 'INPUT_INTAKE.json'))
write_json(list(status = 'RUNNING', pid = Sys.getpid(), utc = utc_now(), planned_new_states = 1700L), file.path(work, 'LAUNCH.json'))
started = utc_now()
diagnostics = list()
for (g in seq_len(nrow(metadata))) {
    name = metadata$graph[g]
    edges = read.delim(file.path(work, 'graphs', paste0(name, '.tsv')), stringsAsFactors = FALSE)
    left = match(edges$gene1, prepared$genes)
    right = match(edges$gene2, prepared$genes)
    stopifnot(!anyNA(left), !anyNA(right), all(left != right), nrow(edges) == metadata$edges[g])
    adjacency = Matrix::sparseMatrix(i = c(left, right), j = c(right, left), x = 1, dims = c(7399L, 7399L), dimnames = list(prepared$genes, prepared$genes))
    transition = dawn_sparse_transition(adjacency)
    damping = dawn_sparse_damping(adjacency, 3)
    if (!metadata$new_fit[g]) {
        label = metadata$release[g]
        stopifnot(identical(adjacency, prepared$graphs[[label]]$adjacency), identical(transition, prepared$graphs[[label]]$transition), identical(damping, prepared$graphs[[label]]$damping))
        next
    }
    saveRDS(list(transition = transition, damping = damping), file.path(work, 'native', paste0(name, '_operator.rds')), compress = FALSE)
    writeBin(as.integer(transition@p), file.path(work, 'native', paste0(name, '_p.bin')), size = 4L, endian = 'little')
    writeBin(as.integer(transition@i), file.path(work, 'native', paste0(name, '_i.bin')), size = 4L, endian = 'little')
    writeBin(as.double(transition@x), file.path(work, 'native', paste0(name, '_x.bin')), size = 8L, endian = 'little')
    writeBin(as.double(damping), file.path(work, 'native', paste0(name, '_d.bin')), size = 8L, endian = 'little')
    for (h in 0:1) {
        state = paste0(name, '_H', h)
        evidence = prepared$evidence[[paste0('H', h)]]$values
        vectors = matrix(NA_real_, nrow = 7399L, ncol = 85L, dimnames = list(prepared$genes, prepared$models))
        for (index in seq_along(prepared$models)) {
            sample = prepared$models[[index]]
            values = evidence[, sample]
            result = dawn_instrumented(transition, values, prepared$mutation[, sample], damping, maxit = 10000L, epsilon = 1e-12, patientTag = sample)
            stopifnot(is.finite(result$terminal_update_norm), result$terminal_update_norm < 1e-12)
            scores = as.numeric(result$summaryOutput$Rank)
            numerical = certificate(transition, damping, values / sum(values), scores)
            eligible = which(prepared$mutation[, sample] == 1)
            ordered = eligible[order(-scores[eligible], method = 'radix')]
            complete = length(ordered) <= 10L
            gap = if (complete) NA_real_ else scores[ordered[10L]] - scores[ordered[11L]]
            vectors[, index] = scores
            diagnostics[[length(diagnostics) + 1L]] = data.frame(state = state, sample_id = sample, iterations = result$convergenceIterations, terminal_update_l2 = result$terminal_update_norm, residual_l1 = numerical$residual_l1, roundoff_l1 = numerical$roundoff_l1, operator_l1_upper = numerical$operator_l1_upper, error_bound = numerical$error_bound, maximum_column_nnz = numerical$maximum_column_nnz, minimum_positive_product = numerical$minimum_positive_product, cutoff_gap = gap, complete_menu = complete, certified = complete || gap > 2 * numerical$error_bound)
        }
        writeBin(as.double(vectors), file.path(work, 'native', paste0(state, '_scores.bin')), size = 8L, endian = 'little')
        write.csv(do.call(rbind, diagnostics), file.path(work, 'cutoff_certificates.csv'), row.names = FALSE)
        write_json(list(status = 'RUNNING', completed_new_states = length(diagnostics), planned_new_states = 1700L, last_state = state, utc = utc_now()), file.path(work, 'PROGRESS.json'))
    }
}
write_json(list(status = 'NUMERICAL_COMPLETE_PENDING_INDEPENDENT_CHECK', started_utc = started, finished_utc = utc_now(), states = length(diagnostics), certified = sum(vapply(diagnostics, function(row) row$certified, logical(1))), baseline_graphs_exact = TRUE, baseline_fits_reused = 340L), file.path(work, 'NATIVE_EXIT.json'))
