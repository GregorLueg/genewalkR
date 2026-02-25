# param wrappers ---------------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Wrapper function for the node2vec parameters
#'
#' @description
#' These are general node2vec parameters, designed for speed. They will default
#' to `available cores - 2 (min 1)`.
#'
#' @param p Numeric. Return parameter for biased random walks. Defaults to
#' `1.0`.
#' @param q Numeric. In-out parameter for biased random walks. Defaults to
#' `1.0`.
#' @param walks_per_node Integer. Number of random walks per node. Defaults to
#' `40L`.
#' @param walk_length Integer. Length of each random walk. Defaults to `40L`.
#' @param batch_size Integer. Batch size for training. Defaults to `256L`.
#' @param n_epochs Integer. Number of training epochs. Defaults to `20L`.
#' @param n_negatives Integer. Number of negative samples. Defaults to `5L`.
#' @param window_size Integer. Context window size. Defaults to `2L`.
#' @param lr Numeric. Learning rate. Defaults to `1e-2`.
#' @param num_workers Optional integer. If kept to `NULL`, it will default
#' to `available cores - 2 (min 1)`.
#'
#' @returns A list with the node2vec parameters.
#'
#' @export
params_node2vec <- function(
  p = 1.0,
  q = 1.0,
  walks_per_node = 40L,
  walk_length = 40L,
  batch_size = 256L,
  n_epochs = 20L,
  n_negatives = 5L,
  window_size = 2L,
  lr = 1e-2,
  num_workers = NULL
) {
  # checkmate
  checkmate::qassert(p, "N1")
  checkmate::qassert(q, "N1")
  checkmate::qassert(walks_per_node, "I1")
  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(n_epochs, "I1")
  checkmate::qassert(n_negatives, "I1")
  checkmate::qassert(window_size, "I1")
  checkmate::qassert(lr, "N1")
  checkmate::qassert(num_workers, c("I1", "0"))

  num_workers <- if (is.null(num_workers)) {
    max(1L, parallel::detectCores() - 2L)
  } else {
    num_workers
  }

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

## genewalk node2vec -----------------------------------------------------------

#' Wrapper function for the GeneWalk-specific node2vec parameters
#'
#' @param p Numeric. Return parameter for biased random walks. Defaults to
#' `1.0`.
#' @param q Numeric. In-out parameter for biased random walks. Defaults to
#' `1.0`.
#' @param walks_per_node Integer. Number of random walks per node. Defaults to
#' `100L`.
#' @param walk_length Integer. Length of each random walk. Defaults to `10L`.
#' @param batch_size Integer. Batch size for training. Defaults to `256L`.
#' @param n_epochs Integer. Number of training epochs. Defaults to `5L`.
#' @param n_negatives Integer. Number of negative samples. Defaults to `5L`.
#' @param window_size Integer. Context window size. Defaults to `1L`. For
#' GeneWalk you want embedding from direct neighbours.
#' @param lr Numeric. Learning rate. Defaults to `1e-2`.
#' @param num_workers Integer. Number of worker threads. Defaults to `1L` (for
#' determinism. If set to ≥1, the fitting is much faster at cost of
#' deterministic results.
#'
#' @returns A list with the (gene-walk) specific node2vec parameters.
#'
#' @export
params_genewalk <- function(
  p = 1.0,
  q = 1.0,
  walks_per_node = 100L,
  walk_length = 10L,
  batch_size = 256L,
  n_epochs = 5L,
  n_negatives = 5L,
  window_size = 1L,
  lr = 1e-2,
  num_workers = 1L
) {
  # checkmate
  checkmate::qassert(p, "N1")
  checkmate::qassert(q, "N1")
  checkmate::qassert(walks_per_node, "I1")
  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(n_epochs, "I1")
  checkmate::qassert(n_negatives, "I1")
  checkmate::qassert(window_size, "I1")
  checkmate::qassert(lr, "N1")
  checkmate::qassert(num_workers, "I1")

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

#' GeneWalk synthetic data params
#'
#' @param n_signal_genes Integer. Number of signal genes. Defaults to `300L`.
#' @param n_noise_genes Integer. Number of noise genes. Defaults to `300L`.
#' @param n_roots Integer. Number of ontology root nodes. Defaults to `10L`.
#' @param depth Integer. Depth of each ontology subtree. Defaults to `3L`.
#' @param branching Integer. Base branching factor per ontology node. Defaults
#' to `3L`.
#' @param p_lateral Numeric. Probability of a lateral edge between sibling
#' nodes. Defaults to `0.1`.
#' @param p_ppi Numeric. Probability of a PPI edge between any two genes in
#' the same group. Defaults to `0.05`.
#' @param min_annotations Integer. Minimum ontology term annotations per gene.
#' Defaults to `10L`.
#' @param max_annotations Integer. Maximum ontology term annotations per gene.
#' Defaults to `25L`.
#' @param min_noise_subtrees Integer. Minimum distinct ontology subtrees a noise
#' gene must span. Defaults to `10L`.
#'
#' @returns A list with GeneWalk parameters.
#'
#' @export
params_genewalk_data <- function(
  n_signal_genes = 150L,
  n_noise_genes = 150L,
  n_roots = 10L,
  depth = 3L,
  branching = 3L,
  p_lateral = 0.1,
  p_ppi = 0.05,
  min_annotations = 10L,
  max_annotations = 25L,
  min_noise_subtrees = 15L
) {
  # checks
  checkmate::qassert(n_signal_genes, "I1")
  checkmate::qassert(n_noise_genes, "I1")
  checkmate::qassert(n_roots, "I1")
  checkmate::qassert(depth, "I1")
  checkmate::qassert(branching, "I1")
  checkmate::qassert(p_lateral, "N1")
  checkmate::qassert(p_ppi, "N1")
  checkmate::qassert(min_annotations, "I1")
  checkmate::qassert(max_annotations, "I1")
  checkmate::qassert(min_noise_subtrees, "I1")

  # return
  list(
    n_signal_genes = n_signal_genes,
    n_noise_genes = n_noise_genes,
    n_roots = n_roots,
    depth = depth,
    branching = branching,
    p_lateral = p_lateral,
    p_ppi = p_ppi,
    min_annotations = min_annotations,
    max_annotations = max_annotations,
    min_noise_subtrees = min_noise_subtrees
  )
}
