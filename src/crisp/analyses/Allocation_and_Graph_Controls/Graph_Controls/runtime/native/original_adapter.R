# DawnRank-derived routines, adapted for CRISP. Original method: Jack Hou.
# SPDX-License-Identifier: GPL-3.0-only
# Source and licence details: licenses/DawnRank.md in the repository root.
dawn_sparse_transition = function(adjacency) {
    stopifnot(all(adjacency@x %in% c(0, 1)))
    degrees = Matrix::colSums(adjacency)
    transition = adjacency %*% Matrix::Diagonal(x = (1 / (degrees + 1e-16)))
    dimnames(transition) = dimnames(adjacency)
    return(transition)
}

dawn_sparse_damping = function(adjacency, mu) {
    degrees = Matrix::colSums(adjacency)
    return(degrees / (degrees + mu))
}

dawn_sparse = function(transition, expression, mutation, damping, maxit = 100, epsilon = 0.0001, patientTag = 'defaultPatient') {
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
    return(list(summaryOutput = output, mutatedRanks = mutated, convergenceIterations = iteration))
}
