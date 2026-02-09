# node2vec tests ---------------------------------------------------------------

## test data -------------------------------------------------------------------

barbell_data <- node2vec_test_data(
  test_type = "barbell",
  n_nodes_per_cluster = 20L
)

caveman_data <- node2vec_test_data(
  test_type = "caveman"
)

stochastic_data <- node2vec_test_data(
  test_type = "stochastic_block"
)

## tests -----------------------------------------------------------------------

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

### caveman graph --------------------------------------------------------------

caveman_res <- node2vec(
  graph_dt = caveman_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L),
  .verbose = FALSE
)

caveman_metrics <- evaluate_node2vec_test(
  embeddings = caveman_res,
  node_labels = caveman_data$node_labels
)

expect_true(
  current = caveman_metrics$within_cluster_sim >= 0.3,
  info = "caveman graph - within cluster similarity"
)

expect_true(
  current = caveman_metrics$between_cluster_sim <
    caveman_metrics$within_cluster_sim,
  info = "caveman graph - between cluster similarity"
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

caveman_res_2 <- node2vec(
  graph_dt = caveman_data$edges,
  node2vec_params = params_node2vec(n_epochs = 25L),
  .verbose = FALSE
)

# there should be a decent'ish mean absolute correlation despite the race
# conditions

expect_true(
  current = mean(abs(diag(cor(caveman_res, caveman_res_2)))) >= 0.7,
  info = "reproducibility of the seeds"
)
