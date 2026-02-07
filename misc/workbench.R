devtools::load_all()

library(data.table)
library(magrittr)
library(ggplot2)

devtools::load_all()

test_data <- synthetic_genewalk_data(
  n_genes = 500L,
  n_pathways = 250L,
  n_communities = 3L,
  ppi_m = 4L
)

test <- genewalkR_class(test_data$edges, list())

# run node2vec to generate the initial embeddings
test <- generate_initial_emb(
  test,
  node2vec_params = params_node2vec(
    window_size = 2L,
    n_epochs = 20L,
    lr = 1e-2,
    walks_per_node = 25L,
    walk_length = 25L
  ),
  backend = "tch-cpu"
)

# generate random permutations
test <- generate_permuted_emb(test)

# calculate the test statistics for gene <> pathway pairs
test <- calculate_genewalk_stats(
  test,
  gene_nodes = intersect(test_data$ground_truth$gene, rownames(test@embd)),
  pathway_nodes = intersect(
    test_data$pathway_metadata$pathway,
    rownames(test@embd)
  )
)

stats <- get_stats(test)

ground_truth_data <- get_expected_associations(
  synthetic_data = test_data,
  return_negatives = TRUE
)

positive_examples <- merge(
  ground_truth_data[(expected_signal)],
  stats,
  by = c("gene", "pathway")
)

hist(
  positive_examples$similarity,
  xlab = "Cosine Similarity",
  main = "Positive pathways"
)
hist(
  positive_examples$avg_pval,
  xlab = "p-value",
  main = "Positive pathways"
)

negative_examples <- merge(
  ground_truth_data[!(expected_signal)],
  stats,
  by = c("gene", "pathway")
)

hist(
  negative_examples$similarity,
  xlab = "Cosine Similarity",
  main = "Negative pathways"
)
hist(negative_examples$avg_pval, xlab = "p-value", main = "Negative pathways")

# generate the internal db -----------------------------------------------------

devtools::load_all()

gw_factory <- GeneWalkGenerator$new()

gw_factory$add_ppi(source = "combined")
gw_factory$add_ppi(source = "intact")
gw_factory$add_pathways()

gw_factory$ppi_sources

gw_factory$reset_choices()

gw_factory$build()

gois <- as.character(gw_factory$return_full_network_dt()[, from][1:100])

gw_obj <- gw_factory$create_for_genes(genes = gois)

gw_obj@graph_dt

devtools::install(".")

library(genewalkR)

system("otool -L $(find . -name 'genewalkR.so' | head -1)")
