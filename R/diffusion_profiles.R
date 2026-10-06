# diffusion profiles -----------------------------------------------------------

## edge helpers ----------------------------------------------------------------

#' Deduplicate edges
#'
#' @description
#' Collapses duplicate edges into one. The constrained PageRank keeps parallel
#' edges, so a pair listed twice (e.g. the same gene pair from two interaction
#' sources, or `A-B` plus `B-A` in an undirected graph) gets twice the
#' transition weight. Run this before [DiffusionProfiles()] unless that is
#' what you want.
#'
#' @param graph_dt data.table. The edge table with `"from"`, `"to"` and
#' optionally `"weight"` (non-negative).
#' @param node_dt Optional data.table with the columns `"id"` and `"type"`. If
#' supplied, the removed duplicates are reported per node-type pair.
#' @param directed Boolean. If `FALSE`, `A-B` and `B-A` are the same edge.
#' Defaults to `FALSE`.
#' @param weight_agg String. How to combine the weights of duplicates. One of
#' `c("max", "sum", "mean")`. Ignored without a `"weight"` column. Defaults to
#' `"max"`, matching the deduplication in [node2vec()].
#' @param .verbose Boolean. Controls verbosity. Defaults to `TRUE`.
#'
#' @return A data.table with the columns `"from"`, `"to"` and, if present in
#' the input, `"weight"`. All other columns are dropped. For undirected graphs
#' each edge is stored once with `from <= to` (lexicographically). Endpoints
#' are returned as character.
#'
#' @export
dedup_edges <- function(
  graph_dt,
  node_dt = NULL,
  directed = FALSE,
  weight_agg = c("max", "sum", "mean"),
  .verbose = TRUE
) {
  weight_agg <- match.arg(weight_agg)

  checkmate::assertDataTable(graph_dt)
  checkmate::assertNames(names(graph_dt), must.include = c("from", "to"))
  checkmate::assertDataTable(node_dt, null.ok = TRUE)
  if (!is.null(node_dt)) {
    checkmate::assertNames(names(node_dt), must.include = c("id", "type"))
  }
  checkmate::qassert(directed, "B1")
  checkmate::assertChoice(weight_agg, c("max", "sum", "mean"))
  checkmate::qassert(.verbose, "B1")
  has_weight <- "weight" %in% names(graph_dt)
  if (has_weight) {
    checkmate::assertNumeric(
      graph_dt$weight,
      lower = 0,
      finite = TRUE,
      any.missing = FALSE,
      .var.name = "graph_dt$weight"
    )
  }

  from <- as.character(graph_dt$from)
  to <- as.character(graph_dt$to)
  edges <- if (directed) {
    data.table::data.table(from = from, to = to)
  } else {
    data.table::data.table(from = pmin(from, to), to = pmax(from, to))
  }

  if (.verbose) {
    dups <- edges[duplicated(edges, by = c("from", "to"))]
    if (nrow(dups) == 0L) {
      message("No duplicate edges found.")
    } else if (is.null(node_dt)) {
      message(sprintf("Removed %i duplicate edge(s).", nrow(dups)))
    } else {
      node_ids <- as.character(node_dt$id)
      node_types <- as.character(node_dt$type)
      type_from <- node_types[match(dups$from, node_ids)]
      type_to <- node_types[match(dups$to, node_ids)]
      type_pair <- if (directed) {
        paste(type_from, type_to, sep = "-")
      } else {
        paste(pmin(type_from, type_to), pmax(type_from, type_to), sep = "-")
      }
      counts <- sort(table(type_pair), decreasing = TRUE)
      message(sprintf(
        "Removed %i duplicate edge(s):\n%s",
        nrow(dups),
        paste(sprintf("  %s: %i", names(counts), counts), collapse = "\n")
      ))
    }
  }

  if (!has_weight) {
    return(unique(edges, by = c("from", "to")))
  }

  edges[, weight := as.numeric(graph_dt$weight)]
  agg_fun <- switch(weight_agg, max = max, sum = sum, mean = mean)
  res <- edges[, .(weight = agg_fun(weight)), by = .(from, to)]

  return(res)
}

## profile generation ----------------------------------------------------------

#' Generate diffusion profiles
#'
#' @description
#' Runs one constrained personalised PageRank per seed (Ruiz et al.). The graph
#' is built once in Rust and the seeds run in parallel. Each call replaces the
#' stored profiles, so pass drugs and diseases together if you want to compare
#' them afterwards.
#'
#' @param object A `DiffusionProfiles` object, see [DiffusionProfiles()].
#' @param seeds Character vector. Node ids to start the walks from; one
#' profile per seed.
#' @param diffusion_profile_params Named list. The PageRank parameters, see
#' [params_diffusion_profiles()].
#' @param .verbose Boolean. Controls verbosity. Defaults to `TRUE`.
#'
#' @return The `DiffusionProfiles` object with `profiles` populated as a
#' numeric matrix of nodes x seeds. Each column sums to 1.
#'
#' @export
generate_profiles <- S7::new_generic(
  name = "generate_profiles",
  dispatch_args = "object",
  fun = function(
    object,
    seeds,
    diffusion_profile_params = params_diffusion_profiles(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_profiles DiffusionProfiles
#'
#' @export
S7::method(generate_profiles, DiffusionProfiles) <- function(
  object,
  seeds,
  diffusion_profile_params = params_diffusion_profiles(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, DiffusionProfiles))
  node_dt <- S7::prop(object, "node_dt")
  node_ids <- as.character(node_dt$id)
  checkmate::assertCharacter(
    seeds,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE
  )
  checkmate::assertSubset(seeds, node_ids)
  assertDiffusionProfilesParams(diffusion_profile_params)
  checkmate::qassert(.verbose, "B1")

  graph_dt <- S7::prop(object, "graph_dt")
  type_weights <- S7::prop(object, "type_weights")

  profiles <- rs_diffusion_profiles(
    node_types = as.character(node_dt$type),
    from = match(as.character(graph_dt$from), node_ids),
    to = match(as.character(graph_dt$to), node_ids),
    weights = if ("weight" %in% names(graph_dt)) {
      as.numeric(graph_dt$weight)
    } else {
      NULL
    },
    type_weight_names = names(type_weights),
    type_weight_values = if (is.null(type_weights)) {
      NULL
    } else {
      as.numeric(type_weights)
    },
    sink_types = S7::prop(object, "sink_types"),
    seeds = as.list(match(seeds, node_ids)),
    directed = S7::prop(object, "params")[["directed"]],
    diffusion_profile_params = diffusion_profile_params
  )
  dimnames(profiles) <- list(node_ids, seeds)

  if (.verbose) {
    message(sprintf(
      "Generated %i diffusion profile(s) over %i nodes.",
      length(seeds),
      length(node_ids)
    ))
  }

  S7::prop(object, "profiles") <- profiles
  S7::prop(object, "params")[["profiles"]] <- diffusion_profile_params

  return(object)
}

## profile distances -----------------------------------------------------------

#' Calculate distances between diffusion profiles
#'
#' @description
#' Distances between two sets of stored profiles, e.g. drugs against diseases.
#' Lower means more similar. Ruiz et al. rank drug-disease pairs by the
#' correlation distance on the multiscale interactome and by Canberra when
#' comparing against gene expression signatures.
#'
#' @details
#' Every profile sums to 1 and the seed keeps roughly `1 - alpha` of the mass
#' on itself. What each metric weights:
#' \itemize{
#'   \item `"correlation"`: Pearson on the full vector. Driven by the large
#'   entries, i.e. the immediate neighbourhood of the seed and hubs.
#'   \item `"cosine"`: Without centring. Since all profiles sum to 1 and are
#'   mostly near zero, it ranks almost identically to `"correlation"`.
#'   \item `"l1"`: Twice the total variation distance between the two
#'   distributions. Less dominated by the top entries than the above.
#'   \item `"l2"`: Dominated by the largest entries; favours big,
#'   well-connected nodes.
#'   \item `"canberra"`: Relative differences per node, so the many
#'   near-zero entries count as much as the large ones. Useful when the
#'   profile is compared against a signature on a different scale.
#' }
#'
#' @param object A `DiffusionProfiles` object with generated profiles.
#' @param from Character vector. Seeds for the rows of the result.
#' @param to Character vector. Seeds for the columns of the result.
#' @param metric Character. One of `"correlation"`, `"canberra"`, `"l1"`,
#' `"l2"` or `"cosine"`. Defaults to `"correlation"`.
#'
#' @return A numeric matrix of `length(from) x length(to)` distances.
#'
#' @references Ruiz, Zitnik and Leskovec, Identification of disease treatment
#' mechanisms through the multiscale interactome, Nat Commun 2021.
#'
#' @export
calculate_profile_distances <- S7::new_generic(
  name = "calculate_profile_distances",
  dispatch_args = "object",
  fun = function(
    object,
    from,
    to,
    metric = c("correlation", "canberra", "l1", "l2", "cosine")
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_profile_distances DiffusionProfiles
#'
#' @export
S7::method(calculate_profile_distances, DiffusionProfiles) <- function(
  object,
  from,
  to,
  metric = c("correlation", "canberra", "l1", "l2", "cosine")
) {
  metric <- match.arg(metric)

  checkmate::assertTRUE(S7::S7_inherits(object, DiffusionProfiles))
  checkmate::assertChoice(
    metric,
    choices = c("correlation", "canberra", "l1", "l2", "cosine")
  )

  profiles <- S7::prop(object, "profiles")
  if (is.null(profiles)) {
    stop("No profiles found. Run generate_profiles() first.")
  }
  checkmate::assertCharacter(from, min.len = 1L, any.missing = FALSE)
  checkmate::assertSubset(from, colnames(profiles))
  checkmate::assertCharacter(to, min.len = 1L, any.missing = FALSE)
  checkmate::assertSubset(to, colnames(profiles))

  res <- rs_profile_distances(
    mat_a = profiles[, from, drop = FALSE],
    mat_b = profiles[, to, drop = FALSE],
    metric = metric
  )
  dimnames(res) <- list(from, to)

  return(res)
}
