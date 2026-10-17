# DawnRank-derived routines, adapted for CRISP. Original method: Jack Hou.
# SPDX-License-Identifier: GPL-3.0-only
# Source and licence details: licenses/DawnRank.md in the repository root.
require_runtime = function() {
    suppressPackageStartupMessages(library(Matrix))
    stopifnot(getRversion() == "4.5.2", as.character(packageVersion("Matrix")) == "1.7.4")
}

json_text = function(value) {
    escape = function(x) {
        x = gsub("\\", "\\\\", x, fixed = TRUE)
        x = gsub('"', '\\"', x, fixed = TRUE)
        x = gsub("\n", "\\n", x, fixed = TRUE)
        x = gsub("\r", "\\r", x, fixed = TRUE)
        x = gsub("\t", "\\t", x, fixed = TRUE)
        paste0('"', x, '"')
    }
    if (is.null(value)) return("null")
    if (is.list(value)) {
        contents = vapply(value, json_text, character(1))
        if (!is.null(names(value))) return(paste0("{", paste0(escape(names(value)), ":", contents, collapse = ","), "}"))
        return(paste0("[", paste(contents, collapse = ","), "]"))
    }
    if (length(value) != 1L) return(paste0("[", paste(vapply(as.list(unname(value)), json_text, character(1)), collapse = ","), "]"))
    if (is.na(value)) return("null")
    if (is.character(value)) return(escape(value))
    if (is.logical(value)) return(if (value) "true" else "false")
    if (is.numeric(value)) return(if (is.finite(value)) sprintf("%.17g", value) else "null")
    stop("Unsupported JSON value")
}

write_json = function(value, path) {
    writeLines(json_text(value), path, useBytes = TRUE)
}

utc_now = function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

project_inherited = function(original, native_genes, genes, models) {
    stopifnot(is.list(original), all(c("values", "evaluable", "sample_status") %in% names(original)))
    stopifnot(is.matrix(original$values), is.numeric(original$values), identical(rownames(original$values), native_genes), identical(colnames(original$values), models))
    stopifnot(anyDuplicated(native_genes) == 0L, anyDuplicated(genes) == 0L, all(genes %in% native_genes), identical(names(original$evaluable), models))
    stopifnot(is.logical(original$evaluable), !anyNA(original$evaluable), length(original$sample_status) == length(models))
    original_status = as.character(original$sample_status)
    if (!is.null(names(original$sample_status))) stopifnot(identical(names(original$sample_status), models))
    projected = original$values[genes, models, drop = FALSE]
    records = vector("list", length(models))
    consumed = vector("list", length(models))
    names(consumed) = models
    for (index in seq_along(models)) {
        sample = models[[index]]
        available = original$evaluable[[index]]
        full_finite = all(is.finite(original$values[, index]))
        stopifnot(!available || full_finite)
        values = projected[, index]
        finite = all(is.finite(values))
        nonnegative = finite && all(values >= 0)
        mass = if (finite) sum(values) else NA_real_
        ready = available && finite && nonnegative && is.finite(mass) && mass > 0
        reason = if (!available) "ORIGINAL_NORMALISATION_UNAVAILABLE" else if (!finite) "PROJECTED_NONFINITE" else if (!nonnegative) "PROJECTED_NEGATIVE" else if (!is.finite(mass)) "NONFINITE_L1_MASS" else if (mass <= 0) "NONPOSITIVE_L1_MASS" else "INPUT_READY"
        if (ready) {
            consumed[[sample]] = values / sum(values)
            stopifnot(all(is.finite(consumed[[sample]])))
        }
        records[[index]] = data.frame(sample_id = sample, original_evaluable = available, original_status = original_status[[index]], full_native_finite = full_finite, projected_finite = finite, projected_nonnegative = nonnegative, projected_l1_mass = mass, input_ready = ready, reason = reason, stringsAsFactors = FALSE)
    }
    return(list(values = projected, consumed = consumed, records = do.call(rbind, records)))
}

dawn_instrumented = function(transition, expression, mutation, damping, maxit = 100, epsilon = 0.0001, patientTag = 'defaultPatient') {
    ranking = matrix(expression / sum(expression), length(expression), 1)
    constant = (expression / sum(expression)) * (1 - damping)
    for (iteration in seq_len(maxit)) {
        updated = as.matrix((damping * (transition %*% ranking)) + constant)
        magnitude = sqrt(sum((updated - ranking) ^ 2))
        ranking = updated
        if (magnitude < epsilon) {
            break
        }
    }
    output = data.frame(Rank = ranking, PercentRank = (100 * rank(ranking) / (length(ranking) + 1)), isMutated = mutation)
    if (sum(output$isMutated) > 0) {
        mutated = output[(output$isMutated == 1), ]
        mutated = mutated[, -3]
        mutated = data.frame(Gene = rownames(mutated), Patient = patientTag, mutated)
        rownames(mutated) = NULL
    } else {
        mutated = logical(0)
    }
    return(list(summaryOutput = output, mutatedRanks = mutated, convergenceIterations = iteration, terminal_update_norm = magnitude, termination_reason = if (!is.finite(magnitude)) "NONFINITE_UPDATE" else if (magnitude < epsilon) "NATIVE_UPDATE_THRESHOLD" else "MAXIMUM_ITERATIONS_EXHAUSTED", fixed_point_error_certified = FALSE))
}
