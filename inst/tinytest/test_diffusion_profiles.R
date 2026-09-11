# diffusion profiles -----------------------------------------------------------

## test data -------------------------------------------------------------------

node_dt <- data.table::data.table(
  id = c("d0", "d1", "p0", "p1", "p2", "g0", "s0"),
  type = c("drug", "drug", "protein", "protein", "protein", "go", "disease")
)
graph_dt <- data.table::data.table(
  from = c("d0", "d0", "d1", "p0", "p1", "p0", "p2", "s0"),
  to = c("p0", "p1", "p1", "p1", "p2", "g0", "g0", "p2")
)
type_weights <- c(drug = 3, protein = 2, go = 1.5, disease = 4)
sink_types <- c("drug", "disease")
alpha <- 0.85

# dense power iteration of eqs 1-3 in Ruiz et al., incl. the sink teleport
dense_profile <- function(seed) {
  ids <- node_dt$id
  types <- node_dt$type
  n <- length(ids)
  adj <- matrix(FALSE, n, n, dimnames = list(ids, ids))
  adj[cbind(graph_dt$from, graph_dt$to)] <- TRUE
  adj[cbind(graph_dt$to, graph_dt$from)] <- TRUE
  is_sink <- types %in% sink_types & ids != seed
  adj[, seed] <- FALSE
  adj[is_sink, ] <- FALSE

  m <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nb <- which(adj[i, ])
    if (!length(nb)) {
      next
    }
    counts <- table(types[nb])
    denom <- sum(type_weights[names(counts)])
    for (j in nb) {
      m[i, j] <- type_weights[[types[j]]] / denom / counts[[types[j]]]
    }
  }

  s <- as.numeric(ids == seed)
  r <- s
  for (k in seq_len(5000L)) {
    r <- (1 - alpha) *
      s +
      alpha * s * sum(r[is_sink]) +
      alpha * as.vector(r %*% m)
  }
  stats::setNames(r / sum(r), ids)
}

ref_dist <- list(
  l2 = \(a, b) sqrt(sum((a - b)^2)),
  l1 = \(a, b) sum(abs(a - b)),
  canberra = \(a, b) {
    den <- abs(a) + abs(b)
    sum(ifelse(den > 0, abs(a - b) / den, 0))
  },
  cosine = \(a, b) 1 - abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2)),
  correlation = \(a, b) 1 - stats::cor(a, b)
)

## tests -----------------------------------------------------------------------

### profiles -------------------------------------------------------------------

obj <- DiffusionProfiles(
  graph_dt = graph_dt,
  node_dt = node_dt,
  sink_types = sink_types,
  type_weights = type_weights
)

expect_warning(get_profiles(obj), info = "no profiles yet")
expect_error(
  calculate_profile_distances(obj, from = "d0", to = "s0"),
  info = "distances before profiles"
)

obj <- generate_profiles(
  obj,
  seeds = node_dt$id,
  diffusion_profile_params = params_diffusion_profiles(
    alpha = alpha,
    max_iter = 10000L,
    tol = 1e-14
  ),
  .verbose = FALSE
)
profiles <- get_profiles(obj)

expect_equal(dim(profiles), c(7L, 7L), info = "nodes x seeds")
expect_equal(colnames(profiles), node_dt$id, info = "seed order")
expect_equal(
  unname(colSums(profiles)),
  rep(1, 7L),
  tolerance = 1e-10,
  info = "profiles sum to 1"
)
expect_true(profiles["d0", "d0"] < 0.5, info = "seed is a source, not a sink")
expect_true(profiles["d1", "d0"] > 0, info = "other drugs receive mass")

for (seed in c("d0", "d1", "p1", "s0")) {
  expect_equal(
    profiles[, seed],
    dense_profile(seed),
    tolerance = 1e-8,
    info = sprintf("dense reference, seed %s", seed)
  )
}

### distances ------------------------------------------------------------------

from <- c("d0", "d1")
to <- c("p0", "p1", "p2", "g0", "s0")

for (metric in names(ref_dist)) {
  res <- calculate_profile_distances(obj, from = from, to = to, metric = metric)
  expected <- outer(
    from,
    to,
    Vectorize(\(f, t) ref_dist[[metric]](profiles[, f], profiles[, t]))
  )
  dimnames(expected) <- list(from, to)
  expect_equal(res, expected, tolerance = 1e-10, info = metric)
}

### plain random walk ----------------------------------------------------------

plain <- generate_profiles(
  DiffusionProfiles(graph_dt, node_dt, sink_types = character(0)),
  seeds = "d0",
  .verbose = FALSE
)
expect_equal(
  sum(get_profiles(plain)),
  1,
  tolerance = 1e-10,
  info = "no type weights, no sinks"
)

### input errors ---------------------------------------------------------------

expect_error(
  DiffusionProfiles(graph_dt, node_dt, sink_types = "cell"),
  info = "unknown sink type"
)
expect_error(
  DiffusionProfiles(graph_dt, node_dt, sink_types, type_weights = c(drug = 1)),
  info = "type weights missing types"
)
expect_error(
  DiffusionProfiles(
    data.table::data.table(from = "d0", to = "zz"),
    node_dt,
    sink_types
  ),
  info = "edge endpoint not in node_dt"
)
expect_error(
  generate_profiles(obj, seeds = "nope", .verbose = FALSE),
  info = "unknown seed"
)
expect_error(
  params_diffusion_profiles(alpha = 1),
  info = "alpha must be below 1"
)
