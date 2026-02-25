# test internal data -----------------------------------------------------------

## main tables -----------------------------------------------------------------

### genes ----------------------------------------------------------------------

gene_info <- get_gene_info(head_only = TRUE)

expect_true(
  checkmate::testDataTable(gene_info),
  info = "gene info correctly returned"
)

expect_true(
  checkmate::testNames(
    names(gene_info),
    must.include = c("ensembl_id", "symbol", "biotype")
  ),
  info = "correct column for gene info"
)

### reactome -------------------------------------------------------------------

reactome_info <- get_reactome_info(head_only = TRUE)

expect_true(
  checkmate::testDataTable(reactome_info),
  info = "reactome info correctly returned"
)

expect_true(
  checkmate::testNames(
    names(reactome_info),
    must.include = c("reactome_id", "reactome_name")
  ),
  info = "correct colums for reactome info"
)

### gene ontology --------------------------------------------------------------

go_info <- get_gene_ontology_info(head_only = TRUE)

expect_true(
  checkmate::testDataTable(go_info),
  info = "gene ontology info correctly returned"
)

expect_true(
  checkmate::testNames(
    names(go_info),
    must.include = c("go_id", "go_name", "name_space")
  ),
  info = "correct colums for gene ontology info"
)

## pathways --------------------------------------------------------------------

### gene ontology --------------------------------------------------------------

go_to_gene <- get_gene_to_go(head_only = TRUE)

expect_true(
  checkmate::testDataTable(go_to_gene),
  info = "gene ontology to genes correctly returned"
)

expect_true(
  checkmate::testNames(
    names(go_to_gene),
    must.include = c("from", "to")
  ),
  info = "correct colums for gene ontology to genes"
)

go_to_go <- get_gene_ontology_hierarchy(head_only = TRUE)

expect_true(
  checkmate::testDataTable(go_to_go),
  info = "gene ontology to genes correctly returned"
)

expect_true(
  checkmate::testNames(
    names(go_to_go),
    must.include = c("from", "to", "relationship")
  ),
  info = "correct colums for gene ontology hierarchy"
)

### pathways -------------------------------------------------------------------

reactome_to_gene <- get_gene_to_reactome(head_only = TRUE)

expect_true(
  checkmate::testDataTable(reactome_to_gene),
  info = "reactome to genes correctly returned"
)

expect_true(
  checkmate::testNames(
    names(reactome_to_gene),
    must.include = c("from", "to")
  ),
  info = "correct colums for reactome to genes"
)

reactome_to_reactome <- get_reactome_hierarchy(head_only = TRUE)

reactome_to_reactome_filtered <- get_reactome_hierarchy(
  head_only = TRUE,
  relationship = "child_of"
)

expect_true(
  checkmate::testDataTable(reactome_to_reactome),
  info = "reactome to reactome correctly returned"
)

expect_true(
  checkmate::testNames(
    names(reactome_to_reactome),
    must.include = c("from", "to", "relationship")
  ),
  info = "correct colums for reactome hierarchy"
)

expect_true(
  current = length(unique(reactome_to_reactome_filtered$relationship)) == 1,
  info = "only child_of returned for reactome"
)

## interactions ----------------------------------------------------------------

### string ---------------------------------------------------------------------

string <- get_interactions_string(head_only = TRUE)

string_hq <- get_interactions_string(head_only = TRUE, threshold = 0.85)

expect_true(
  checkmate::testDataTable(string),
  info = "string interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(string),
    must.include = c("from", "to", "string_score", "string_high_quality")
  ),
  info = "correct colums for string"
)

expect_true(
  all(string_hq$string_score >= 0.85),
  info = "string filtering behaving"
)

### intact ---------------------------------------------------------------------

intact <- get_interactions_intact(head_only = TRUE)

intact_filtered <- get_interactions_intact(head_only = TRUE, threshold = 0.85)

intact_physical <- get_interactions_intact(
  head_only = TRUE,
  physical_interactions = TRUE
)

expect_true(
  checkmate::testDataTable(intact_filtered),
  info = "intact interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(intact),
    must.include = c("from", "to", "intact_score", "physical_association")
  ),
  info = "correct colums for string"
)

expect_true(
  all(intact_filtered$intact_score >= 0.85),
  info = "intact filtering behaving"
)

expect_true(
  all(intact_physical$physical_association),
  info = "intact filtering for physical behaving"
)

### other networks -------------------------------------------------------------

reactome_network <- get_interactions_reactome(head_only = TRUE)

expect_true(
  checkmate::testDataTable(reactome_network),
  info = "reactome interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(reactome_network),
    must.include = c("from", "to")
  ),
  info = "correct colums for reactome network"
)


signor_network <- get_interactions_signor(head_only = TRUE)

expect_true(
  checkmate::testDataTable(signor_network),
  info = "signor interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(signor_network),
    must.include = c("from", "to")
  ),
  info = "correct colums for signor network"
)

pc_network <- get_interactions_pc(head_only = TRUE)

expect_true(
  checkmate::testDataTable(pc_network),
  info = "pathway commons interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(pc_network),
    must.include = c("from", "to")
  ),
  info = "correct colums for pathway commons network"
)

combined_network <- get_interactions_combined(head_only = TRUE)

expect_true(
  checkmate::testDataTable(combined_network),
  info = "combined interactions correctly returned"
)

expect_true(
  checkmate::testNames(
    names(combined_network),
    must.include = c("from", "to")
  ),
  info = "correct colums for combined network"
)
