# integration tests for gene walk ----------------------------------------------

## synthetic data --------------------------------------------------------------

# generate synthetic data
gene_walk_test_data <- synthetic_genewalk_data(
  data_params = params_genewalk_data(
    n_signal_genes = 150L,
    n_noise_genes = 150L
  )
)

## tests -----------------------------------------------------------------------

### object generation and getters ----------------------------------------------

genewalk_obj <- GeneWalk(
  graph_dt = gene_walk_test_data$full_data,
  gene_to_pathway_dt = gene_walk_test_data$gene_to_pathways,
  gene_ids = gene_walk_test_data$gene_ids,
  pathway_ids = gene_walk_test_data$pathway_ids
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
  genewalk_params = params_genewalk(),
  .verbose = FALSE
)

embeddings <- get_embedding(genewalk_obj)

expect_true(
  current = checkmate::testList(
    embeddings,
    types = "matrix"
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
    "sem_sim",
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

signal_genes_res <- gw_stats[grepl("signal", gene)]

noise_genes_res <- gw_stats[grepl("noise", gene)]

expect_true(
  current = mean(signal_genes_res$similarity) >
    mean(noise_genes_res$similarity),
  info = "similarities of the signal genes is better"
)

expect_true(
  current = mean(signal_genes_res$avg_pval) < mean(noise_genes_res$avg_pval),
  info = "p-values of the signal genes is better"
)
