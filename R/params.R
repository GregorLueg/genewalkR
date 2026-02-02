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
#' @param n_epochs Integer. Number of training epochs. Defaults to `15L`.
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
  n_epochs = 15L,
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
