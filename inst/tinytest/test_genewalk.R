# integration tests for gene walk ----------------------------------------------

## synthetic data --------------------------------------------------------------

# generate synthetic data
gene_walk_test_data <- synthetic_genewalk_data(
  ppi_params = params_ppi(),
  pathway_params = params_pathway(),
  n_communities = 3L
)

## tests -----------------------------------------------------------------------

### object generation and getters ----------------------------------------------

genewalk_obj <- GeneWalk(
  graph_dt = gene_walk_test_data$edges,
  graph_gene_params = list()
)

expect_warning(
  current = get_embedding(genewalk_obj),
  info = "get_embedding() warning"
)

expect_warning(
  current = get_stats(genewalk_obj),
  info = "get_stats() warning"
)

### embedding generation -------------------------------------------------------

genewalk_obj <- generate_initial_emb(
  genewalk_obj,
  node2vec_params = params_node2vec(
    window_size = 2L,
    n_epochs = 20L,
    lr = 1e-2,
    walks_per_node = 25L,
    walk_length = 25L
  ),
  .verbose = FALSE
)

embeddings <- get_embedding(genewalk_obj)

expect_true(
  current = checkmate::testMatrix(
    embeddings,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "initial embedding generation working"
)

### permutation and statistical testing ----------------------------------------

# generate random permutations
genewalk_obj <- generate_permuted_emb(genewalk_obj, .verbose = FALSE)

expect_true(
  current = checkmate::testList(
    genewalk_obj@permuted_embd,
    types = "numeric",
    len = 3
  ),
  info = "permutation generation working"
)

# calculate the test statistics for gene <> pathway pairs
genewalk_obj <- calculate_genewalk_stats(
  genewalk_obj,
  gene_nodes = gene_walk_test_data$ground_truth$gene,
  pathway_nodes = gene_walk_test_data$pathway_metadata$pathway,
  .verbose = FALSE
)

gw_stats <- get_stats(genewalk_obj)

expect_true(
  current = checkmate::testDataTable(gw_stats),
  info = "gw stats is a data.table"
)

expect_equal(
  current = colnames(gw_stats),
  target = c(
    "gene",
    "pathway",
    "similarity",
    "avg_pval",
    "pval_ci_lower",
    "pval_ci_upper",
    "avg_global_fdr",
    "global_fdr_ci_lower",
    "global_fdr_ci_upper",
    "avg_gene_fdr",
    "gene_fdr_ci_lower",
    "gene_fdr_ci_upper"
  ),
  info = "gw_stats expected columns"
)

### assumption testing ---------------------------------------------------------

ground_truth_data <- get_expected_associations(
  synthetic_data = gene_walk_test_data,
  return_negatives = TRUE
)

positive_examples <- merge(
  ground_truth_data[(expected_signal)],
  gw_stats,
  by = c("gene", "pathway")
)

negative_examples <- merge(
  ground_truth_data[!(expected_signal)],
  gw_stats,
  by = c("gene", "pathway")
)

expect_true(
  current = mean(positive_examples$similarity) >
    mean(negative_examples$similarity),
  info = "gw stats - positive examples have higher similiarities"
)

expect_true(
  current = mean(-log10(positive_examples$avg_pval)) >
    mean(-log10(negative_examples$avg_pval)),
  info = "gw stats - positive examples have lower significance"
)
