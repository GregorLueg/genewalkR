# embd drift -------------------------------------------------------------------

## synthetic data --------------------------------------------------------------

embd_drift_data <- differential_graph_test_data()

## class generation ------------------------------------------------------------

embd_drift_obj <- EmbedDrift(
  graph_dt_1 = embd_drift_data$g1_edges,
  graph_dt_2 = embd_drift_data$g2_edges
)

expect_warning(
  current = get_embeddings(embd_drift_obj),
  info = "embd drift: warning with get embeddings"
)

expect_warning(
  current = get_aligned_embd(embd_drift_obj),
  info = "embd drift: warning with get aligned embedding"
)

expect_warning(
  current = get_stats(embd_drift_obj),
  info = "embd drift: warning with get statistics"
)

## test embeddings -------------------------------------------------------------

embd_drift_obj <- generate_initial_embeddings(
  object = embd_drift_obj,
  node2vec_params = params_node2vec(window_size = 2L),
  embd_dim = 8L,
  .verbose = FALSE
)

embd_ls <- get_embeddings(embd_drift_obj)

expect_true(
  current = checkmate::testList(embd_ls, types = "matrix", len = 2),
  info = "embd drift: list of matrices returned"
)

embd_1 <- get_embeddings(embd_drift_obj, which = 1L)
embd_2 <- get_embeddings(embd_drift_obj, which = 2L)

expect_true(
  current = checkmate::testMatrix(
    embd_1,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "embd drift: first matrix returned"
)

expect_true(
  current = checkmate::testMatrix(
    embd_2,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "embd drift: second matrix returned"
)

## calculate statistics --------------------------------------------------------

embd_drift_obj <- calculate_drift(
  object = embd_drift_obj
)

stat_results <- get_stats(embd_drift_obj)

expect_true(
  current = checkmate::testDataTable(stat_results),
  info = "embd drift: stats data.table returned"
)

expect_equal(
  current = colnames(stat_results),
  target = c("node", "cosine_similarity", "node_status"),
  info = "correct: stats data.table has expected column names"
)

aligned_embd <- get_aligned_embd(embd_drift_obj)

expect_true(
  current = checkmate::testMatrix(
    aligned_embd,
    mode = "numeric"
  ),
  info = "embd drift: aligned_matrix returned"
)

expected_results <- merge(
  stat_results,
  embd_drift_data$data_info[, c("node", "is_differential")],
  by = "node"
)

similarity_neg <- expected_results[
  !(is_differential) & node_status == "shared",
  mean(cosine_similarity)
]

similarity_pos <- expected_results[
  (is_differential) & node_status == "shared",
  mean(cosine_similarity)
]

expect_true(
  current = similarity_pos < similarity_neg,
  info = "embd drift: genes with different context have lower cosine similarity"
)

expect_true(
  current = all(
    expected_results[node_status == "g1_only", cosine_similarity] == -1.1
  ),
  info = "embd drift: g1 genes have expected values"
)

expect_true(
  current = all(
    expected_results[node_status == "g2_only", cosine_similarity] == 1.1
  ),
  info = "embd drift: g2 genes have expected values"
)
