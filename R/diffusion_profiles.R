# diffusion profiles -----------------------------------------------------------

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
    seeds = match(seeds, node_ids),
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
