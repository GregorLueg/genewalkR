# network diffusion methods ----------------------------------------------------

## kernel building -------------------------------------------------------------

#' Build the diffusion kernel
#'
#' @description
#' Constructs the graph Laplacian and runs the spectral decomposition,
#' then applies the chosen kernel transform. The result is stored as an
#' external pointer for subsequent calls to `diffuse()`.
#'
#' @param object A `DiffusionScores` object, see [DiffusionScores()].
#' @param kernel Character. One of `"regularised_laplacian"`,
#' `"commute_time"`, `"diffusion"`, `"inverse_cosine"`, `"pstep"`.
#' Defaults to `"regularised_laplacian"`.
#' @param kernel_params Named list. Kernel-specific parameters, see
#' [params_kernel()].
#' @param normalised Logical. Use the normalised Laplacian. Defaults to
#' `TRUE`.
#' @param strategy Character. `"full"` or `"truncated"`. Defaults to
#' `"truncated"`.
#' @param k Integer. Number of eigenvalues for the truncated strategy.
#' Ignored if `strategy = "full"`. Defaults to `100L`.
#' @param .verbose Logical. Controls verbosity. Defaults to `TRUE`.
#'
#' @return The `DiffusionScores` object with `kernel_pointer` and `params`
#' populated.
#'
#' @export
build_kernel <- S7::new_generic(
  name = "build_kernel",
  dispatch_args = "object",
  fun = function(
    object,
    kernel = c(
      "regularised_laplacian",
      "commute_time",
      "diffusion",
      "inverse_cosine",
      "pstep"
    ),
    kernel_params = params_kernel(),
    normalised = TRUE,
    strategy = "truncated",
    k = 100L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method build_kernel DiffusionScores
#'
#' @export
S7::method(build_kernel, DiffusionScores) <- function(
  object,
  kernel = c(
    "regularised_laplacian",
    "commute_time",
    "diffusion",
    "inverse_cosine",
    "pstep"
  ),
  kernel_params = params_kernel(),
  normalised = TRUE,
  strategy = "truncated",
  k = 100L,
  .verbose = TRUE
) {
  kernel <- match.arg(kernel)

  checkmate::assertTRUE(S7::S7_inherits(object, DiffusionScores))
  checkmate::assertChoice(
    kernel,
    choices = c(
      "regularised_laplacian",
      "commute_time",
      "diffusion",
      "inverse_cosine",
      "pstep"
    )
  )
  assertKernelParams(kernel_params)
  checkmate::qassert(normalised, "B1")
  checkmate::assertChoice(strategy, choices = c("full", "truncated"))
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(.verbose, "B1")

  g <- S7::prop(object, "graph")
  adj <- igraph::as_adjacency_matrix(g, sparse = TRUE)

  ptr <- rs_diffusion_kernel(
    i = adj@i,
    p = adj@p,
    x = adj@x,
    n = nrow(adj),
    kernel = kernel,
    kernel_params = kernel_params,
    normalised = normalised,
    strategy = strategy,
    k = k,
    node_names = rownames(adj),
    verbose = .verbose
  )

  S7::prop(object, "kernel_pointer") <- ptr
  S7::prop(object, "params")[["kernel"]] <- list(
    kernel = kernel,
    kernel_params = kernel_params,
    normalised = normalised,
    strategy = strategy,
    k = k
  )

  return(object)
}

## diffuse scores --------------------------------------------------------------

#' Run diffusion scoring
#'
#' @description
#' Diffuses input scores over the precomputed kernel. Requires
#' `build_kernel()` to have been called first.
#'
#' @param object A `DiffusionScores` object with a built kernel.
#' @param input Named numeric vector or named matrix. Names must correspond
#' to node names in the graph. For a vector, treated as a single input
#' column. For a matrix, rows are background nodes and columns are
#' separate input sets.
#' @param background Character vector. Node names forming the background
#' set. Defaults to all nodes in `input`.
#' @param method Character. One of `"raw"`, `"z"`, or `"mc"`. Defaults to
#' `"raw"`.
#' @param n_perm Integer. Number of permutations for `method = "mc"`.
#' Defaults to `1000L`.
#' @param seed Integer. RNG seed for `method = "mc"`. Defaults to `42L`.
#'
#' @return The `DiffusionScores` object with `scores` populated as a
#' `data.table` with columns `node` and one score column per input.
#'
#' @export
diffuse_scores <- S7::new_generic(
  name = "diffuse_scores",
  dispatch_args = "object",
  fun = function(
    object,
    input,
    background = NULL,
    method = "raw",
    n_perm = 1000L,
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#' @method diffuse_scores DiffusionScores
#'
#' @export
S7::method(diffuse_scores, DiffusionScores) <- function(
  object,
  input,
  background = NULL,
  method = "raw",
  n_perm = 1000L,
  seed = 42L
) {
  checkmate::assertTRUE(S7::S7_inherits(object, DiffusionScores))
  checkmate::assertChoice(method, choices = c("raw", "z", "mc"))
  checkmate::qassert(n_perm, "I1[1,)")
  checkmate::qassert(seed, "I1")

  ptr <- S7::prop(object, "kernel_pointer")
  if (is.null(ptr)) {
    stop("Kernel not built. Call build_kernel() first.")
  }

  # coerce vector to single-column matrix
  if (is.numeric(input) && is.null(dim(input))) {
    nms <- names(input)
    input <- matrix(input, ncol = 1L)
    rownames(input) <- nms
  }
  checkmate::assertMatrix(input, mode = "numeric", min.rows = 1)
  checkmate::assertCharacter(rownames(input), min.len = 1L)

  node_names <- ptr$node_names # from the Rust struct
  if (is.null(background)) {
    background <- rownames(input)
  }
  checkmate::assertCharacter(background, min.len = 2L)
  checkmate::assertSubset(background, node_names)
  checkmate::assertSubset(rownames(input), background)

  # map background node names to 0-based indices
  bkgd_idx <- match(background, node_names) - 1L

  # build the score matrix: rows = background, cols = inputs
  # ensure input rows are ordered to match background
  score_mat <- matrix(0, nrow = length(background), ncol = ncol(input))
  input_pos <- match(rownames(input), background)
  score_mat[input_pos, ] <- input

  n_bkgd <- length(background)
  n_inputs <- ncol(input)

  raw_out <- rs_diffuse(
    kernel = ptr,
    scores = as.vector(score_mat), # column-major
    n_bkgd = n_bkgd,
    n_inputs = n_inputs,
    bkgd_indices = bkgd_idx,
    method = method,
    n_perm = n_perm,
    seed = seed
  )

  # reshape: Rust returns column-major (n x n_inputs)
  n <- length(node_names)
  res_mat <- matrix(raw_out, nrow = n, ncol = n_inputs)

  col_names <- if (!is.null(colnames(input))) {
    colnames(input)
  } else {
    sprintf("input_%i", seq_len(n_inputs))
  }

  dt <- data.table::data.table(node = node_names)
  for (j in seq_len(n_inputs)) {
    data.table::set(dt, j = col_names[j], value = res_mat[, j])
  }

  S7::prop(object, "scores") <- dt
  S7::prop(object, "params")[["diffusion"]] <- list(
    method = method,
    n_perm = n_perm,
    seed = seed,
    background = background
  )

  return(object)
}
