# embed drift main functions ---------------------------------------------------

## generate embeddings ---------------------------------------------------------

#' Generate initial node embeddings for both graphs
#'
#' @description
#' Runs node2vec independently on both graphs using identical hyperparameters,
#' storing the resulting embedding matrices in the object. Graph 1 serves as
#' the reference for subsequent Procrustes alignment. For details on the
#' node2vec algorithm, see Grover & Leskovec. For the approach generally
#' speaking, see Amir, et al.
#'
#' @param object An `EmbedDrift` object, see [EmbedDrift()].
#' @param embd_dim Integer. Embedding dimension size. Defaults to `64L`.
#' @param node2vec_params Named list. Node2vec training parameters, see
#' [params_node2vec()].
#' @param directed Boolean. Whether to treat both graphs as directed. Defaults
#' to `FALSE`.
#' @param seed Integer. Random seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity. Defaults to `TRUE`.
#'
#' @return The `EmbedDrift` object with `embd_1`, `embd_2`, and `params`
#' populated.
#'
#' @export
#'
#' @references
#' Grover & Leskovec, arXiv, 2016 \cr
#' Amir, et al., Cancer Cell, 2025
generate_initial_embeddings <- S7::new_generic(
  name = "generate_initial_embeddings",
  dispatch_args = "object",
  fun = function(
    object,
    embd_dim = 64L,
    node2vec_params = params_node2vec(),
    directed = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_initial_embeddings EmbedDrift
#'
#' @export
S7::method(generate_initial_embeddings, EmbedDrift) <- function(
  object,
  embd_dim = 64L,
  node2vec_params = params_node2vec(),
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, EmbedDrift))
  assertNode2VecParam(node2vec_params)
  checkmate::qassert(embd_dim, "I1[1,)")
  checkmate::qassert(directed, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # internal helper
  .embed_graph <- function(graph_dt, seed) {
    nodes <- sort(unique(c(graph_dt$from, graph_dt$to)))
    from_idx <- match(graph_dt$from, nodes)
    to_idx <- match(graph_dt$to, nodes)
    weights <- if ("weight" %in% names(graph_dt)) graph_dt$weight else NULL

    m <- rs_node2vec(
      from = from_idx,
      to = to_idx,
      weights = weights,
      node2vec_params = node2vec_params,
      embd_dim = embd_dim,
      directed = directed,
      seed = seed,
      verbose = .verbose
    )
    rownames(m) <- nodes
    colnames(m) <- sprintf("dim_%i", seq_len(embd_dim))
    m
  }

  g1 <- S7::prop(object, "graph_dt_1")
  g2 <- S7::prop(object, "graph_dt_2")

  if (.verbose) {
    message("Embedding graph 1...")
  }
  embd_1 <- .embed_graph(g1, seed)

  if (.verbose) {
    message("Embedding graph 2...")
  }
  embd_2 <- .embed_graph(g2, seed + 1L)

  full_params <- node2vec_params
  full_params[["embd_dim"]] <- embd_dim
  full_params[["directed"]] <- directed

  S7::prop(object, "embd_1") <- embd_1
  S7::prop(object, "embd_2") <- embd_2
  S7::prop(object, "params")[["node2vec"]] <- full_params

  return(object)
}

## calculate drift -------------------------------------------------------------

#' Calculate embedding drift between two graphs
#'
#' @description
#' Aligns the embedding of graph 2 onto the reference embedding of graph 1
#' via orthogonal Procrustes, then computes per-node cosine similarities as a
#' measure of context drift. Nodes shared between both graphs receive a cosine
#' similarity in `[-1, 1]`. Nodes exclusive to graph 1 are assigned `-1.1` and
#' nodes exclusive to graph 2 are assigned `+1.1` as sentinel values indicating
#' absence rather than a computed similarity.
#'
#' @param object An `EmbedDrift` object, see [EmbedDrift()].
#'
#' @return The `EmbedDrift` object with `aligned_embd` and `stats` populated.
#' `stats` is a `data.table` with columns `node`, `cosine_similarity`, and
#' `node_status` (`"shared"`, `"g1_only"`, `"g2_only"`).
#'
#' @export
calculate_drift <- S7::new_generic(
  name = "calculate_drift",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method calculate_drift EmbedDrift
#'
#' @export
S7::method(calculate_drift, EmbedDrift) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, EmbedDrift))

  # get the objects
  embd_1 <- S7::prop(object, "embd_1")
  embd_2 <- S7::prop(object, "embd_2")

  # early return
  if (length(embd_1) == 0 || length(embd_2) == 0) {
    warning(paste(
      "Embeddings are missing. Returning object as is.",
      "Please run generate_initial_embeddings() first."
    ))
    return(object)
  }

  # get all of the nodes
  nodes_1 <- rownames(embd_1)
  nodes_2 <- rownames(embd_2)
  shared <- intersect(nodes_1, nodes_2)
  g1_only <- setdiff(nodes_1, nodes_2)
  g2_only <- setdiff(nodes_2, nodes_1)

  # calculate the procrustes
  zeallot::`%<-%`(
    c(aligned, sims),
    rs_procrustes_align(
      embd1 = embd_1[shared, , drop = FALSE],
      embd2 = embd_2[shared, , drop = FALSE]
    )
  )

  stats <- data.table::rbindlist(list(
    data.table::data.table(
      node = shared,
      cosine_similarity = sims,
      node_status = "shared"
    ),
    data.table::data.table(
      node = g1_only,
      cosine_similarity = -1.1,
      node_status = "g1_only"
    ),
    data.table::data.table(
      node = g2_only,
      cosine_similarity = 1.1,
      node_status = "g2_only"
    )
  ))

  S7::prop(object, "aligned_embd") <- aligned
  S7::prop(object, "stats") <- stats

  return(object)
}
