# param wrappers ---------------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Wrapper function for the node2vec parameters
#'
#' @param p Numeric. Return parameter for biased random walks. Defaults to
#' `1.0`.
#' @param q Numeric. In-out parameter for biased random walks. Defaults to
#' `1.0`.
#' @param walks_per_node Integer. Number of random walks per node. Defaults to
#' `40L`.
#' @param walk_length Integer. Length of each random walk. Defaults to `40L`.
#' @param num_workers Integer. Number of worker threads. Defaults to `4L`.
#' @param batch_size Integer. Batch size for training. Defaults to `256L`.
#' @param n_epochs Integer. Number of training epochs. Defaults to `20L`.
#' @param n_negatives Integer. Number of negative samples. Defaults to `5L`.
#' @param window_size Integer. Context window size. Defaults to `2L`.
#' @param lr Numeric. Learning rate. Defaults to `1e-2`.
#'
#' @returns A list with the node2vec parameters.
#'
#' @export
params_node2vec <- function(
  p = 1.0,
  q = 1.0,
  walks_per_node = 40L,
  walk_length = 40L,
  num_workers = 4L,
  batch_size = 256L,
  n_epochs = 20L,
  n_negatives = 5L,
  window_size = 2L,
  lr = 1e-2
) {
  # checkmate
  checkmate::qassert(p, "N1")
  checkmate::qassert(q, "N1")
  checkmate::qassert(walks_per_node, "I1")
  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(num_workers, "I1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(n_epochs, "I1")
  checkmate::qassert(n_negatives, "I1")
  checkmate::qassert(window_size, "I1")
  checkmate::qassert(lr, "N1")

  list(
    p = p,
    q = q,
    walks_per_node = walks_per_node,
    walk_length = walk_length,
    num_workers = num_workers,
    batch_size = batch_size,
    n_epochs = n_epochs,
    n_negatives = n_negatives,
    window_size = window_size,
    lr = lr
  )
}

## genewalk synthetic data -----------------------------------------------------

#' PPI parameters
#'
#' @param n_genes Integer. Number of genes. Defaults to `500L`.
#' @param ppi_m Integer. Number of edges to attach in preferential attachment.
#' Defaults to `3L`.
#' @param min_community_size Integer. Minimum size for gene communities.
#' Defaults to `20L`.
#'
#' @returns A list with PPI parameters.
#'
#' @export
params_ppi <- function(
  n_genes = 500L,
  ppi_m = 3L,
  min_community_size = 20L
) {
  checkmate::qassert(n_genes, "I1")
  checkmate::qassert(ppi_m, "I1")
  checkmate::qassert(min_community_size, "I1")

  list(
    n_genes = n_genes,
    ppi_m = ppi_m,
    min_community_size = min_community_size
  )
}

#' Pathway parameters
#'
#' @param n_pathways Integer. Maximum number of pathways. Defaults to `250L`.
#' You can undershoot the number here due to the building of the DAG.
#' @param pathway_depth Integer. Maximum depth of pathway hierarchy. Defaults
#' to `4L`.
#' @param pathway_branching Integer. Average branching factor. Defaults to `3L`.
#' @param n_focal_pathways Integer. Number of focal pathways per community.
#' Defaults to `15L`.
#' @param connections_per_gene Integer. Number of pathway connections per gene.
#' Defaults to `3L`.
#'
#' @returns A list with pathway parameters.
#'
#' @export
params_pathway <- function(
  n_pathways = 250L,
  pathway_depth = 4L,
  pathway_branching = 3L,
  n_focal_pathways = 15L,
  connections_per_gene = 3L
) {
  checkmate::qassert(n_pathways, "I1")
  checkmate::qassert(pathway_depth, "I1")
  checkmate::qassert(pathway_branching, "I1")
  checkmate::qassert(n_focal_pathways, "I1")
  checkmate::qassert(connections_per_gene, "I1")

  list(
    n_pathways = n_pathways,
    pathway_depth = pathway_depth,
    pathway_branching = pathway_branching,
    n_focal_pathways = n_focal_pathways,
    connections_per_gene = connections_per_gene
  )
}
