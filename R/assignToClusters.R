#' @title assignCluster
#' 
#' @description Assigns clusters to items based on their probabilities and target proportions.
#' @param probs A matrix of probabilities, where each row represents an item and each column represents a cluster.
#' @param target.props A numeric vector of target proportions for each cluster, summing to 1.
#' @return A numeric vector indicating the assigned cluster for each item.

assignCluster <- function(probs, target.props) {
    target.counts <- round(target.props * nrow(probs)) 
    target.counts[which.max(target.counts)] <- target.counts[which.max(target.counts)] + nrow(probs) - sum(target.counts[-length(target.counts)])

    result <- numeric(nrow(probs))
    max.probs <- apply(probs, 1, max)
    best.clusters <- max.col(probs)
    order.idx <- order(max.probs, decreasing = TRUE)
    
    for (i in order.idx) {
        best.cluster <- best.clusters[i]
        if (target.counts[best.cluster] > 0) {
            result[i] <- best.cluster
            target.counts[best.cluster] <- target.counts[best.cluster] - 1
        } else {
            available.clusters <- which(target.counts > 0)
            if (length(available.clusters) > 0) {
                next.best <- available.clusters[which.max(probs[i, available.clusters])]
                result[i] <- next.best
                target.counts[next.best] <- target.counts[next.best] - 1
            } else {
                result[i] <- best.cluster
            }
        }
    }
    return(result)
}

#' @title assignCellsToClusters
#'
#' @description Assigns single cells to clusters based on a probability matrix and target proportions.
#' @param prob.matrix A matrix of probabilities, where each row represents a cell and each column represents a cluster.
#' @param target.proportions A matrix or data frame of target proportions for each cluster, where rows correspond to cell types and columns to clusters.
#' @param sc.obj A Seurat object of single-cell data.
#' @return A list where each element is a vector of cell names assigned to a particular cluster.

assignCellsToClusters <- function(prob.matrix, target.proportions, sc.obj) {
    final.clusters <- numeric(ncol(sc.obj)) %>% `names<-`(colnames(sc.obj))
    for (cell.type in levels(sc.obj)) {
        cell.type.mask <- Idents(sc.obj)[Idents(sc.obj) == cell.type] %>% names
        final.clusters[cell.type.mask] <- assignCluster(prob.matrix[cell.type.mask, ], target.proportions[cell.type, ]) %>% 
		colnames(prob.matrix)[.]
    }
    partioned.res <- split(names(final.clusters), final.clusters)[colnames(prob.matrix)]    
    return(partioned.res)
}
