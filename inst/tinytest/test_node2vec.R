# node2vec tests ---------------------------------------------------------------

## test data -------------------------------------------------------------------

barbell_data <- node2vec_test_data(
  test_type = "barbell",
  n_nodes_per_cluster = 20L
)

cavemen_data <- node2vec_test_data(
  test_type = "cavemen",
  n_clusters = 5L,
  p_between = 0.5
)

stochastic_data <- node2vec_test_data(
  test_type = "stochastic_block"
)

## tests -----------------------------------------------------------------------

### check the underlying data --------------------------------------------------

expect_equal(
  current = unique(barbell_data$node_labels$cluster),
  target = c(1, 2),
  info = "barbell returns 2 clusters"
)

expect_equal(
  current = unique(cavemen_data$node_labels$cluster),
  target = 1:5,
  info = "barbell returns 5 clusters"
)

expect_equal(
  current = unique(stochastic_data$node_labels$cluster),
  target = 1:3,
  info = "barbell returns 3 clusters"
)

### barbell graph --------------------------------------------------------------

barbell_res <- node2vec(
  graph_dt = barbell_data$edges,
  .verbose = FALSE
)

barbell_metrics <- evaluate_node2vec_test(
  embeddings = barbell_res,
  node_labels = barbell_data$node_labels
)

expect_true(
  current = barbell_metrics$within_cluster_sim >= 0.3,
  info = "barbell graph - within cluster similarity"
)

expect_true(
  current = barbell_metrics$between_cluster_sim <
    barbell_metrics$within_cluster_sim,
  info = "barbell graph - between cluster similarity"
)

### cavemn graph --------------------------------------------------------------

cavemen_res <- node2vec(
  graph_dt = cavemen_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L),
  .verbose = FALSE
)

cavemen_metrics <- evaluate_node2vec_test(
  embeddings = cavemen_res,
  node_labels = cavemen_data$node_labels
)

expect_true(
  current = cavemen_metrics$within_cluster_sim >= 0.3,
  info = "cavemen graph - within cluster similarity"
)

expect_true(
  current = cavemen_metrics$between_cluster_sim <
    cavemen_metrics$within_cluster_sim,
  info = "cavemen graph - between cluster similarity"
)

### stochastic data ------------------------------------------------------------

stochastic_res <- node2vec(
  graph_dt = stochastic_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L),
  .verbose = FALSE
)

stochastic_metrics <- evaluate_node2vec_test(
  embeddings = stochastic_res,
  node_labels = stochastic_data$node_labels
)

expect_true(
  current = stochastic_metrics$within_cluster_sim >= 0.3,
  info = "stochastic graph - within cluster similarity"
)

expect_true(
  current = stochastic_metrics$between_cluster_sim <=
    stochastic_metrics$within_cluster_sim,
  info = "stochastic graph - between cluster similarity"
)

### seed reproducibility -------------------------------------------------------

# due to the race condition implemented in the gradient, I need to reduce
# to one thread

cavemen_res_1 <- node2vec(
  graph_dt = cavemen_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L, num_workers = 1L),
  .verbose = FALSE
)

cavemen_res_2 <- node2vec(
  graph_dt = cavemen_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L, num_workers = 1L),
  .verbose = FALSE
)

expect_equal(
  current = cavemen_res_1,
  target = cavemen_res_2,
  info = "reproducibility with n_workers is set to 1"
)
