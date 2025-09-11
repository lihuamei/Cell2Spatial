#' @title assignCellsGreedyClusters
#'
#' @description
#' Greedy assignment of single cells to clusters while respecting per-cluster and per-cell-type
#' target counts. Cells are processed in order of their maximum prediction score, and each cell
#' is placed into the most likely cluster that still has quota for its type.
#'
#' @param netx.pred A numeric matrix of predicted probabilities, with rows representing cells
#'   and columns representing clusters. Row names must be cell IDs, and column names must be cluster IDs.
#' @param cell.types A character or factor vector of cell type labels, one per cell.
#'   The length must equal the number of rows in netx.pred. Names should correspond to cell IDs.
#' @param target.cls.counts A numeric data frame or matrix of required cell counts,
#'   with rows representing clusters and columns representing cell types.
#'   Each entry gives how many cells of a specific type should be assigned to a specific cluster.
#'   Row names must include all cluster IDs from netx.pred.
#'
#' @return A list with the following elements:
#'   - assigned_cluster: a named list mapping each cluster ID to the vector of assigned cell IDs
#'   - assigned_types: a matrix of cluster by cell type showing how many cells of each type were assigned
#'   - target_quota: the original target count matrix aligned to clusters and cell types
#'   - leftover_cap: an integer vector giving remaining unfilled capacity per cluster
#'   - leftover_quota: a matrix of remaining unfilled per-type quotas (cluster by cell type)

assignCellsGreedyClusters <- function(netx.pred, cell.types, target.cls.counts) {
    P <- as.matrix(netx.pred)
    clusters <- colnames(P)

    stopifnot(!is.null(rownames(target.cls.counts)))
    stopifnot(all(clusters %in% rownames(target.cls.counts)))

    ctype <- factor(cell.types)
    if (is.null(names(ctype))) names(ctype) <- rownames(P)
    Tlev <- levels(ctype)

    tgt <- target.cls.counts[clusters, , drop = FALSE]
    K <- length(clusters)
    Tn <- length(Tlev)
    quota <- matrix(0L, nrow = K, ncol = Tn, dimnames = list(clusters, Tlev))
    common_types <- intersect(colnames(tgt), Tlev)
    if (length(common_types) > 0) {
        quota[, common_types] <- as.matrix(tgt[, common_types, drop = FALSE])
    }
    quota <- apply(quota, 2, as.integer)

    N <- nrow(P)
    stopifnot(sum(quota) == N)

    quota_left <- quota
    cap_left <- as.integer(rowSums(quota_left))

    best_score <- apply(P, 1, max)
    order_cells <- order(best_score, decreasing = TRUE)

    assign_idx <- integer(N)
    names(assign_idx) <- rownames(P)

    for (ii in order_cells) {
        t_idx <- as.integer(ctype[ii])

        cand1 <- which(cap_left > 0 & quota_left[, t_idx] > 0)
        if (length(cand1) > 0) {
            k_star <- cand1[which.max(P[ii, cand1])]
        } else {
            cand2 <- which(cap_left > 0)
            if (length(cand2) == 0) next
            k_star <- cand2[which.max(P[ii, cand2])]
        }

        assign_idx[ii] <- k_star
        cap_left[k_star] <- cap_left[k_star] - 1L
        quota_left[k_star, t_idx] <- max(0L, quota_left[k_star, t_idx] - 1L)
    }

    unassigned <- which(assign_idx == 0L)
    if (length(unassigned) > 0 && sum(cap_left) >= length(unassigned)) {
        for (ii in unassigned) {
            cand <- which(cap_left > 0)
            k_star <- cand[which.max(P[ii, cand])]
            assign_idx[ii] <- k_star
            cap_left[k_star] <- cap_left[k_star] - 1L
        }
    }

    cluster_names <- colnames(P)
    assigned_cluster <- lapply(seq_len(K), function(k) {
        cells <- names(assign_idx)[assign_idx == k]
        if (length(cells) == 0) character(0) else cells
    })
    names(assigned_cluster) <- cluster_names

    assigned_types <- do.call(rbind, lapply(seq_len(K), function(k) {
        cs <- assigned_cluster[[k]]
        tab <- if (length(cs)) table(ctype[cs]) else integer(0)
        out <- integer(Tn)
        names(out) <- Tlev
        if (length(tab)) out[names(tab)] <- as.integer(tab)
        out
    }))
    rownames(assigned_types) <- cluster_names

    return(list(
        assigned_cluster = assigned_cluster,
        assigned_types   = assigned_types,
        target_quota     = quota,
        leftover_cap     = cap_left,
        leftover_quota   = quota_left
    ))
}
