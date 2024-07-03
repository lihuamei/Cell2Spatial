assign_cells_to_clusters <- function(prob_matrix, cell_types, target_proportions) {
    n_cells <- nrow(prob_matrix)
    n_clusters <- ncol(prob_matrix)
    n_cell_types <- nrow(target_proportions)

    if (length(cell_types) != n_cells) {
        stop("The length of cell_types must match the number of rows in prob_matrix")
    }

    if (ncol(target_proportions) != n_clusters) {
        stop("The number of columns in target_proportions must match the number of clusters")
    }

    assign_cluster <- function(probs, target_props, n_assign) {
        target_counts <- round(target_props * n_assign)
        target_counts[length(target_counts)] <- n_assign - sum(target_counts[-length(target_counts)])

        result <- numeric(nrow(probs))

        max_probs <- apply(probs, 1, max)
        best_clusters <- max.col(probs)

        order_idx <- order(max_probs, decreasing = TRUE)

        for (i in order_idx) {
            best_cluster <- best_clusters[i]
            if (target_counts[best_cluster] > 0) {
                result[i] <- best_cluster
                target_counts[best_cluster] <- target_counts[best_cluster] - 1
            } else {
                available_clusters <- which(target_counts > 0)
                if (length(available_clusters) > 0) {
                    next_best <- available_clusters[which.max(probs[i, available_clusters])]
                    result[i] <- next_best
                    target_counts[next_best] <- target_counts[next_best] - 1
                } else {
                    result[i] <- best_cluster
                }
            }
        }

        return(result)
    }

    final_clusters <- numeric(n_cells)
    for (cell_type in 1:n_cell_types) {
        cell_type_mask <- cell_types == cell_type
        n_cells_of_type <- sum(cell_type_mask)

        if (n_cells_of_type > 0) {
            final_clusters[cell_type_mask] <- assign_cluster(
                prob_matrix[cell_type_mask, ],
                target_proportions[cell_type, ],
                n_cells_of_type
            )
        }
    }

    actual_proportions <- matrix(0, nrow = n_cell_types, ncol = n_clusters)
    for (i in 1:n_cell_types) {
        for (j in 1:n_clusters) {
            actual_proportions[i, j] <- sum(cell_types == i & final_clusters == j) / sum(cell_types == i)
        }
    }

    mad <- mean(abs(actual_proportions - target_proportions))
    accuracy <- mean(final_clusters == max.col(prob_matrix))

    return(list(
        final_clusters = final_clusters,
        actual_proportions = actual_proportions,
        mad = mad,
        accuracy = accuracy
    ))
}
