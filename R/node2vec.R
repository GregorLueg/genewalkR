# node2vec ---------------------------------------------------------------------

#' Node2Vec
#'
#' @description
#' Runs the Node2Vec algorithm
#'
#' @param graph_dt data.table. The edge data.frame. Needs to have the columns
#' `"from"` and `"to"`. Can also optionally a `"weight"` column.
#' @param embd_dim Integer. Size of the embedding dimensions to create. Defaults
#' to `8L` in line with the authors recommendations.
#' @param node2vec_params Named list. Contains the node2vec parameters. The
#' list has the following elements:
#' \itemize{
#'  \item p - Numeric. Return parameter for biased random walks. Defaults to
#'  `1.0`.
#'  \item q - Numeric. In-out parameter for biased random walks. Defaults to
#'  `1.0.`.
#'  \item walks_per_node - Integer. Number of random walks per node. Defaults to
#'  ` 40L.`.
#'  \item walk_length - Integer. Length of each random walk. Defaults to `40L`.
#'  \item num_workers - Number of worker threads during batching. Defaults to
#'  `4L`. To note: during graph generation, Rust will use via the Rayon backend
#'  all available threads.
#'  \item batch_size - Integer. Batch size for training. Defaults to `256L`.
#'  \item n_epochs - Integer. Number of training epochs. Defaults to `15L`.
#'  \item n_negatives - Integer. Number of negative samples. Defaults to `5L`.
#'  \item window_size - Integer.  Context window size. Defaults to `2L`.
#'  \item lr - Numeric. Learning rate. Defaults to `1e-2`.
#' }
#' @param backend String. One of `c("ndarray", "tch-cpu")`. `"tch-cpu"` is
#' usually much faster and supported only on Unix-based systems.
#' @param directed Boolean. Indicates if this is a directed or undirected
#' network. Defaults to `FALSE`.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function
#'
#' @export
node2vec <- function(
  graph_dt,
  embd_dim = 8L,
  node2vec_params = params_node2vec(),
  backend = c("ndarray", "tch-cpu"),
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  backend <- match.arg(backend)

  # checks
  checkmate::assertDataTable(graph_dt)
  checkmate::assertNames(names(graph_dt), must.include = c("from", "to"))
  checkmate::qassert(embd_dim, "I1")
  assertNode2VecParam(node2vec_params)
  checkmate::qassert(directed, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  checkmate::assertChoice(backend, c("ndarray", "tch-cpu"))

  # funtion body
  weights <- if ("weight" %in% names(graph_dt)) {
    graph_dt[["weights"]]
  } else {
    NULL
  }

  nodes <- unique(c(graph_dt$from, graph_dt$to))
  from_idx <- match(graph_dt$from, nodes)
  to_idx <- match(graph_dt$to, nodes)

  embd <- rs_gene_walk(
    from = from_idx,
    to = to_idx,
    weights = weights,
    gene_walk_params = node2vec_params,
    backend = backend,
    embd_dim = embd_dim,
    directed = directed,
    seed = seed,
    verbose = .verbose
  )

  rownames(embd) <- nodes
  colnames(embd) <- sprintf("emb_dim_%i", 1:embd_dim)

  return(embd)
}
