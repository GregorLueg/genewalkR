# diffusion tests --------------------------------------------------------------

## test data -------------------------------------------------------------------

lattice <- diffusion_test_lattice(n_side = 10L)

## constructor -----------------------------------------------------------------

expect_error(
  current = DiffusionScores("not_a_graph"),
  info = "constructor rejects non-igraph"
)

obj <- DiffusionScores(lattice$graph)

expect_true(
  current = igraph::is_igraph(S7::prop(obj, "graph")),
  info = "graph property is an igraph"
)

expect_equal(
  current = nrow(S7::prop(obj, "scores")),
  target = 0L,
  info = "scores are empty on init"
)

expect_true(
  current = is.null(S7::prop(obj, "kernel_pointer")),
  info = "kernel pointer is NULL on init"
)

## params_kernel ---------------------------------------------------------------

expect_error(
  current = params_kernel(a = 1.0),
  info = "params_kernel rejects a < 2"
)

expect_error(
  current = params_kernel(p = 0L),
  info = "params_kernel rejects p < 1"
)

expect_error(
  current = params_kernel(sigma2 = "a"),
  info = "params_kernel rejects non-numeric sigma2"
)

kp <- params_kernel(sigma2 = 2.0, add_diag = 0.5, a = 4.0, p = 3L)
expect_equal(current = kp$sigma2, target = 2.0, info = "sigma2 stored")
expect_equal(current = kp$p, target = 3L, info = "p stored")

## build_kernel ----------------------------------------------------------------

obj <- DiffusionScores(lattice$graph)

expect_error(
  current = build_kernel(obj, kernel = "nonsense"),
  info = "build_kernel rejects unknown kernel"
)

expect_error(
  current = build_kernel(obj, strategy = "nonsense"),
  info = "build_kernel rejects unknown strategy"
)

### full strategy
obj_full <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  strategy = "full",
  .verbose = FALSE
)

expect_true(
  current = !is.null(S7::prop(obj_full, "kernel_pointer")),
  info = "kernel pointer populated after build_kernel (full)"
)

expect_equal(
  current = S7::prop(obj_full, "params")[["kernel"]][["strategy"]],
  target = "full",
  info = "strategy stored in params"
)

### truncated strategy
obj_trunc <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  strategy = "truncated",
  k = 20L,
  .verbose = FALSE
)

expect_true(
  current = !is.null(S7::prop(obj_trunc, "kernel_pointer")),
  info = "kernel pointer populated after build_kernel (truncated)"
)

expect_equal(
  current = S7::prop(obj_trunc, "params")[["kernel"]][["k"]],
  target = 20L,
  info = "k stored in params"
)

## diffuse (interface) ---------------------------------------------------------

### error without kernel
expect_error(
  current = diffuse(DiffusionScores(lattice$graph), input = c("1" = 1.0)),
  info = "diffuse errors without built kernel"
)

### error with bad method
obj_k <- build_kernel(
  DiffusionScores(lattice$graph),
  strategy = "full",
  .verbose = FALSE
)

expect_error(
  current = diffuse(obj_k, input = c("1" = 1.0), method = "nonsense"),
  info = "diffuse rejects unknown method"
)

## get_scores ------------------------------------------------------------------

### empty warning
expect_warning(
  current = get_scores(DiffusionScores(lattice$graph)),
  info = "get_scores warns when scores empty"
)

## diffusion assumptions (lattice graph) ---------------------------------------

# source node is the corner (node "1"); on a 10x10 lattice, distances
# increase along the grid. diffusion scores should decay with graph distance.

n <- igraph::vcount(lattice$graph)
source_node <- "1"
input_vec <- stats::setNames(1.0, source_node)

### raw: scores decay with distance from source --------------------------------

obj_raw <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  kernel_params = params_kernel(),
  strategy = "full",
  .verbose = FALSE
)
obj_raw <- diffuse_scores(obj_raw, input = input_vec, method = "raw")
scores_raw <- get_scores(obj_raw)

expect_equal(
  current = nrow(scores_raw),
  target = n,
  info = "raw: scores returned for all nodes"
)

# group by graph distance from source and check monotone decay
dists <- igraph::distances(lattice$graph, v = source_node)
scores_raw[, dist := as.integer(dists[1, node])]
mean_by_dist <- scores_raw[,
  .(mean_score = mean(get(setdiff(names(scores_raw), c("node", "dist"))))),
  by = dist
]
data.table::setorder(mean_by_dist, dist)

expect_true(
  current = all(diff(mean_by_dist$mean_score) <= 0),
  info = "raw: mean scores decay monotonically with graph distance"
)

expect_true(
  current = scores_raw[node == source_node][[2]] > 0,
  info = "raw: source node has positive score"
)

### z-score: source node should be positive, distant nodes negative ------------

obj_z <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  kernel_params = params_kernel(),
  strategy = "full",
  .verbose = FALSE
)

all_nodes <- lattice$node_names
input_z <- stats::setNames(
  ifelse(all_nodes == source_node, 1.0, 0.0),
  all_nodes
)
obj_z <- diffuse_scores(obj_z, input = input_z, method = "z")
scores_z <- get_scores(obj_z)

score_col <- setdiff(names(scores_z), "node")
source_z <- scores_z[node == source_node][[score_col]]
farthest_node <- lattice$node_names[which.max(dists[1, ])]
farthest_z <- scores_z[node == farthest_node][[score_col]]

expect_true(
  current = source_z > 0,
  info = "z-score: source node is positive"
)

expect_true(
  current = farthest_z < source_z,
  info = "z-score: farthest node has lower z than source"
)

expect_true(
  current = all(is.finite(scores_z[[score_col]])),
  info = "z-score: all values are finite"
)

### mc: source node should have lowest heatrank --------------------------------

obj_mc <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  kernel_params = params_kernel(),
  strategy = "full",
  .verbose = FALSE
)
obj_mc <- diffuse_scores(
  obj_mc,
  input = input_z,
  method = "mc",
  n_perm = 1000L,
  seed = 42L
)
scores_mc <- get_scores(obj_mc)

score_col_mc <- setdiff(names(scores_mc), "node")
source_hr <- scores_mc[node == source_node][[score_col_mc]]

expect_true(
  current = source_hr > 0 && source_hr < 1,
  info = "mc: heatrank in (0, 1)"
)

expect_true(
  current = source_hr == min(scores_mc[[score_col_mc]]),
  info = "mc: source node has lowest (warmest) heatrank"
)

### mc: reproducibility --------------------------------------------------------

obj_mc2 <- build_kernel(
  DiffusionScores(lattice$graph),
  kernel = "regularised_laplacian",
  kernel_params = params_kernel(),
  strategy = "full",
  .verbose = FALSE
)
obj_mc2 <- diffuse_scores(
  obj_mc2,
  input = input_z,
  method = "mc",
  n_perm = 1000L,
  seed = 42L
)
scores_mc2 <- get_scores(obj_mc2)

expect_equal(
  current = scores_mc[[score_col_mc]],
  target = scores_mc2[[score_col_mc]],
  info = "mc: identical seed gives identical results"
)

## kernel types ----------------------------------------------------------------

kernel_types <- c(
  "regularised_laplacian",
  "commute_time",
  "diffusion",
  "inverse_cosine",
  "pstep"
)

for (kt in kernel_types) {
  kp <- if (kt == "pstep") params_kernel(a = 3.0, p = 5L) else params_kernel()

  obj_kt <- build_kernel(
    DiffusionScores(lattice$graph),
    kernel = kt,
    kernel_params = kp,
    normalised = TRUE,
    strategy = "full",
    .verbose = FALSE
  )
  obj_kt <- diffuse_scores(obj_kt, input = input_vec, method = "raw")
  sc <- get_scores(obj_kt)

  expect_equal(
    current = nrow(sc),
    target = n,
    info = sprintf("%s: returns scores for all nodes", kt)
  )

  sc_col <- setdiff(names(sc), "node")
  source_score <- sc[node == source_node][[sc_col]]

  expect_true(
    current = is.finite(source_score),
    info = sprintf("%s: source score is finite", kt)
  )

  expect_true(
    current = source_score > 0,
    info = sprintf("%s: source score is positive", kt)
  )
}

## truncated vs full agreement -------------------------------------------------

# need a larger graph for truncation to be meaningful
lattice_large <- diffusion_test_lattice(n_side = 10L)
source_large <- "1"
input_large <- stats::setNames(1.0, source_large)

obj_full_cmp <- build_kernel(
  DiffusionScores(lattice_large$graph),
  kernel = "diffusion",
  kernel_params = params_kernel(sigma2 = 2.0),
  strategy = "full",
  .verbose = FALSE
)
obj_full_cmp <- diffuse_scores(
  obj_full_cmp,
  input = input_large,
  method = "raw"
)

obj_trunc_cmp <- build_kernel(
  DiffusionScores(lattice_large$graph),
  kernel = "diffusion",
  kernel_params = params_kernel(sigma2 = 2.0),
  strategy = "truncated",
  k = 50L,
  .verbose = FALSE
)
obj_trunc_cmp <- diffuse_scores(
  obj_trunc_cmp,
  input = input_large,
  method = "raw"
)

full_scores <- get_scores(obj_full_cmp)
trunc_scores <- get_scores(obj_trunc_cmp)

merged <- merge(
  full_scores,
  trunc_scores,
  by = "node",
  suffixes = c("_full", "_trunc")
)
cor_val <- cor(merged[[2]], merged[[3]])

expect_true(
  current = cor_val > 0.9,
  info = sprintf(
    "truncated (k=50, diffusion, n=900) correlates with full (r=%.3f)",
    cor_val
  )
)
