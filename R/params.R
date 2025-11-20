# param wrappers ---------------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Wrapper function for the node2vec parameters
#'
#' @param p Numeric. Return parameter for biased random walks. Defaults to
#' `1.0`.
#' @param q Numeric. In-out parameter for biased random walks. Defaults to
#' `1.0`.
#' @param walks_per_node Integer. Number of random walks per node. Defaults to
#' `20L`.
#' @param walk_length Integer. Length of each random walk. Defaults to `20L`.
#' @param num_workers Integer. Number of worker threads. Defaults to `4L`.
#' @param batch_size Integer. Batch size for training. Defaults to `256L`.
#' @param num_epochs Integer. Number of training epochs. Defaults to `5L`.
#' @param num_negatives Integer. Number of negative samples. Defaults to `5L`.
#' @param window_size Integer. Context window size. Defaults to `2L`.
#' @param lr Numeric. Learning rate. Defaults to `1e-3`.
#'
#' @returns A list with the node2vec parameters.
#'
#' @export
params_node2vec <- function(
  p = 1.0,
  q = 1.0,
  walks_per_node = 20L,
  walk_length = 20L,
  num_workers = 4L,
  batch_size = 256L,
  num_epochs = 5L,
  num_negatives = 5L,
  window_size = 2L,
  lr = 1e-3
) {
  # checkmate
  checkmate::qassert(p, "N1")
  checkmate::qassert(q, "N1")
  checkmate::qassert(walks_per_node, "I1")
  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(num_workers, "I1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(num_epochs, "I1")
  checkmate::qassert(num_negatives, "I1")
  checkmate::qassert(window_size, "I1")
  checkmate::qassert(lr, "N1")

  list(
    p = p,
    q = q,
    walks_per_node = walks_per_node,
    walk_length = walk_length,
    num_workers = num_workers,
    batch_size = batch_size,
    num_epochs = num_epochs,
    num_negatives = num_negatives,
    window_size = window_size,
    lr = lr
  )
}
