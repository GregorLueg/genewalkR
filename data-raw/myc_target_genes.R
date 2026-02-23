# get the MYC target genes -----------------------------------------------------

library(msigdbr)
library(data.table)

# get the hallmark sets --------------------------------------------------------

h_gene_sets <- msigdbr(species = "human", collection = "H")
setDT(h_gene_sets)

myc_genes <- h_gene_sets[gs_name == "HALLMARK_MYC_TARGETS_V1"][, c(
  "ensembl_gene",
  "gene_symbol"
)][, gene_set := "MYC target genes (human)"]

usethis::use_data(myc_genes, overwrite = TRUE)
