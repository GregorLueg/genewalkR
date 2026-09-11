# metapath2vec -----------------------------------------------------------------

#' Metapath2Vec
#'
#' @description
#' Runs metapath2vec over a heterogeneous graph. Random walks are constrained
#' to follow a metapath schema over node types, e.g.
#' `c("gene", "pathway", "gene")`, and a skip-gram model is trained on them.
#' With `metapath_plus = TRUE`, negative samples are drawn from the context
#' node's own type (metapath2vec++).
#'
#' Walks start only on nodes of the metapath's first type. Nodes no walk ever
#' reaches keep their random initialisation, so by default they are removed
#' and listed in the `unvisited_nodes` attribute. Check the `walk_stats`
#' attribute as well: a high share of truncated or dropped walks means the
#' schema does not fit the graph.
#'
#' @param graph_dt data.table. The edge table. Needs to have the columns
#' `"from"` and `"to"`, and can optionally have a `"weight"` column.
#' @param node_dt data.table. The node table with the columns `"id"` and
#' `"type"`. Every edge endpoint needs to be in `id`.
#' @param metapath Character vector. The metapath over node types, closing on
#' its starting type, e.g. `c("gene", "pathway", "gene")`.
#' @param embd_dim Integer. Size of the embedding dimensions to create.
#' Defaults to `8L`.
#' @param metapath2vec_params Named list. The training parameters, see
#' [params_metapath2vec()].
#' @param metapath_plus Boolean. Use per-type negative sampling
#' (metapath2vec++). Defaults to `FALSE`.
#' @param filter_unvisited Boolean. Remove the rows of nodes that no walk
#' visited. Defaults to `TRUE`.
#' @param directed Boolean. Indicates if this is a directed or undirected
#' network. Defaults to `FALSE`.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A numeric matrix of `n_nodes x embd_dim` with the node ids as
#' rownames. Rows are ordered by node type. Carries three attributes:
#' \itemize{
#'   \item node_type - Named factor with the type of each row.
#'   \item unvisited_nodes - Character vector with the ids of the nodes no
#'   walk visited. Removed from the matrix if `filter_unvisited = TRUE`.
#'   \item walk_stats - List with `start_nodes`, `attempted`, `truncated`,
#'   `dropped`, `mean_length` and `walk_length`.
#' }
#'
#' @references Dong, Chawla and Swami, metapath2vec: Scalable Representation
#' Learning for Heterogeneous Networks, KDD 2017.
#'
#' @export
metapath2vec <- function(
  graph_dt,
  node_dt,
  metapath,
  embd_dim = 8L,
  metapath2vec_params = params_metapath2vec(),
  metapath_plus = FALSE,
  filter_unvisited = TRUE,
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertDataTable(graph_dt)
  checkmate::assertNames(names(graph_dt), must.include = c("from", "to"))
  checkmate::assertDataTable(node_dt, min.rows = 1L)
  checkmate::assertNames(names(node_dt), must.include = c("id", "type"))
  checkmate::assertCharacter(
    as.character(node_dt$id),
    any.missing = FALSE,
    unique = TRUE,
    .var.name = "node_dt$id"
  )
  checkmate::assertCharacter(metapath, min.len = 3L, any.missing = FALSE)
  checkmate::qassert(embd_dim, "I1[1,)")
  assertMetapath2VecParams(metapath2vec_params)
  checkmate::qassert(metapath_plus, "B1")
  checkmate::qassert(filter_unvisited, "B1")
  checkmate::qassert(directed, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  node_ids <- as.character(node_dt$id)
  node_types <- as.character(node_dt$type)
  if (any(grepl("-", node_types, fixed = TRUE))) {
    stop("Node types must not contain '-', it separates the metapath.")
  }
  checkmate::assertSubset(metapath, unique(node_types))
  if (metapath[1] != metapath[length(metapath)]) {
    stop(sprintf(
      "The metapath must close on its starting type: starts on '%s', ends on '%s'.",
      metapath[1],
      metapath[length(metapath)]
    ))
  }

  from_idx <- match(as.character(graph_dt$from), node_ids)
  to_idx <- match(as.character(graph_dt$to), node_ids)
  if (anyNA(from_idx) || anyNA(to_idx)) {
    stop("All edge endpoints in `graph_dt` need to be present in `node_dt$id`.")
  }

  weights <- if ("weight" %in% names(graph_dt)) {
    checkmate::assertNumeric(graph_dt$weight, lower = 0, any.missing = FALSE)
    as.numeric(graph_dt$weight)
  } else {
    NULL
  }

  # function body
  res <- rs_metapath2vec(
    node_ids = node_ids,
    node_types = node_types,
    from = from_idx,
    to = to_idx,
    weights = weights,
    metapath = paste(metapath, collapse = "-"),
    metapath_plus = metapath_plus,
    metapath2vec_params = metapath2vec_params,
    embd_dim = embd_dim,
    directed = directed,
    seed = seed,
    verbose = .verbose
  )

  keep <- if (filter_unvisited) res$visited else rep(TRUE, length(res$visited))
  unvisited <- res$node_names[!res$visited]
  if (.verbose && length(unvisited) > 0L) {
    message(sprintf(
      "%i nodes were never visited by a walk%s.",
      length(unvisited),
      if (filter_unvisited) " and were removed" else ""
    ))
  }

  embd <- res$embedding[keep, , drop = FALSE]
  rownames(embd) <- res$node_names[keep]
  colnames(embd) <- sprintf("emb_dim_%i", seq_len(embd_dim))
  node_type <- factor(res$node_types)[keep]
  names(node_type) <- res$node_names[keep]
  attr(embd, "node_type") <- node_type
  attr(embd, "unvisited_nodes") <- unvisited
  attr(embd, "walk_stats") <- res$walk_stats

  return(embd)
}
