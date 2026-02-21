# new synthetic data -----------------------------------------------------------

gene_walk_syn_data <- synthetic_genewalk_data()

genewalk_obj <- GeneWalk(
  graph_dt = gene_walk_syn_data$full_data,
  gene_to_pathway_dt = gene_walk_syn_data$gene_to_pathways,
  gene_ids = gene_walk_syn_data$gene_ids,
  pathway_ids = gene_walk_syn_data$pathway_ids
)

genewalk_obj <- generate_initial_emb(
  genewalk_obj,
  node2vec_params = params_node2vec(),
  .verbose = TRUE
)

genewalk_obj <- generate_permuted_emb(genewalk_obj, .verbose = TRUE)

genewalk_obj <- calculate_genewalk_stats(
  genewalk_obj,
  .verbose = FALSE
)

gw_stats <- get_stats(genewalk_obj)

hist(gw_stats[grepl("signal", gene), similarity], xlim = c(-1, 1))

hist(gw_stats[grepl("noise", gene), similarity])

summary(gw_stats[grepl("signal", gene), similarity])

summary(gw_stats[grepl("noise", gene), similarity])

summary(gw_stats[grepl("signal", gene), avg_pval])
summary(gw_stats[grepl("noise", gene), avg_pval])

hist(gw_stats[grepl("signal", gene), avg_pval])

hist(gw_stats[grepl("noise", gene), avg_pval])

summary(gw_stats[grepl("signal", gene), avg_global_fdr])
summary(gw_stats[grepl("noise", gene), avg_global_fdr])

t.test(
  gw_stats[grepl("signal", gene), similarity],
  gw_stats[grepl("noise", gene), similarity]
)

table(
  grepl("signal", gw_stats$gene),
  gw_stats$avg_global_fdr < 0.05
)

# test on real data ------------------------------------------------------------

go_info <- get_gene_ontology_info()

go_translator <- setNames(go_info$go_name, go_info$go_id)

go_ontology <- get_gene_ontology_hierarchy()

go_genes <- get_gene_to_go()

combined_network <- get_interactions_combined()

library(msigdbr)

h_gene_sets <- msigdbr(species = "human", collection = "H")

gene_translator <- setNames(h_gene_sets$gene_symbol, h_gene_sets$ensembl_gene)

setDT(h_gene_sets)

myc_genes <- h_gene_sets[gs_name == "HALLMARK_MYC_TARGETS_V1"][]

go_genes[from %in% myc_genes]

final_network <- rbindlist(
  list(
    go_genes[from %in% myc_genes],
    go_ontology[, c("from", "to")],
    combined_network[, c("from", "to")][from %in% myc_genes & to %in% myc_genes]
  )
) %>%
  unique() %>%
  .[, `:=`(from = as.character(from), to = as.character(to))]

genes <- unique(c(
  combined_network[, c("from", "to")][
    from %in% myc_genes & to %in% myc_genes
  ] %$%
    c(from, to),
  go_genes[from %in% myc_genes, from]
))

pathways <- unique(c(
  go_ontology[, c("from", "to")] %$% c(from, to),
  go_genes[from %in% myc_genes, to]
))

go_genes[from %in% myc_genes] %$% unique(c(from, to))

genewalk_obj <- GeneWalk(
  graph_dt = final_network,
  gene_to_pathway_dt = go_genes[from %in% myc_genes],
  gene_ids = myc_genes,
  pathway_ids = pathways
)

genewalk_obj <- generate_initial_emb(
  genewalk_obj,
  node2vec_params = params_node2vec(n_epochs = 5L, num_workers = 4L),
  .verbose = TRUE
)

genewalk_obj <- generate_permuted_emb(genewalk_obj, .verbose = TRUE)

genewalk_obj <- calculate_genewalk_stats(
  genewalk_obj,
  .verbose = TRUE
)

res <- get_stats(genewalk_obj)[, `:=`(
  gene_symbol = gene_translator[gene],
  go_name = go_translator[pathway]
)]


devtools::load_all()

builder <- GeneWalkGenerator$new()

builder$add_ppi(source = "combined")

builder$add_pathways()

builder$build()

builder$create_for_genes(genes = myc_genes)

builder[['.__enclos_env__']]$private$gene_walk_data

nrow(builder[['.__enclos_env__']]$private$gene_walk_data)
n_nodes(builder[['.__enclos_env__']]$private$gene_walk_data)

class(builder[['.__enclos_env__']]$private$gene_walk_data)

get_gw_data(builder[['.__enclos_env__']]$private$gene_walk_data)

x <- get_gw_data_filtered(
  builder[['.__enclos_env__']]$private$gene_walk_data,
  gene_ids = myc_genes
)

str(x)
