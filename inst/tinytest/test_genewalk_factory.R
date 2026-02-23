# genewalk factory tests -------------------------------------------------------

## gw data ---------------------------------------------------------------------

### get the data ---------------------------------------------------------------

reactome_genes <- get_gene_to_reactome()
reactome_ppi <- get_interactions_reactome()
reactome_hierarchy <- get_reactome_hierarchy(relationship = "parent_of")
reactome_info <- get_reactome_info()

# for testing purposes
cell_cycle_reactome <- reactome_info[
  grepl("cell cycle", reactome_name),
  reactome_id
]

gois <- reactome_genes[to %in% cell_cycle_reactome, unique(from)]

pathways_with_genes <- reactome_genes[from %in% gois, unique(to)]

### generate the class ---------------------------------------------------------

gw_data <- new_data_builder(
  ppis = reactome_ppi,
  gene_to_pathways = reactome_genes,
  pathway_hierarchy = reactome_hierarchy
)

expect_true(
  current = checkmate::testClass(gw_data, "DataBuilder"),
  info = "gw_data class correctly generated"
)

expect_equal(
  current = nrow(gw_data),
  target = gw_data$n_edges,
  info = "nrow dispatch behaving"
)

expect_true(
  current = is.na(dim(gw_data)[2]),
  info = "second dimensions in dim(DataBuilder) = NA"
)

### methods --------------------------------------------------------------------

#### all data ------------------------------------------------------------------

internal_data <- get_gw_data(x = gw_data)

expect_true(
  current = checkmate::testDataTable(internal_data),
  info = "internal data correctly returned"
)

expect_true(
  current = checkmate::testNames(
    names(internal_data),
    must.include = c("from", "to", "type")
  ),
  info = "internal data has expected columns"
)

#### filtered data -------------------------------------------------------------

internal_data_filtered <- get_gw_data_filtered(x = gw_data, gene_ids = gois)

expect_true(
  current = checkmate::testDataTable(internal_data_filtered$gwn),
  info = "gwn data returned"
)

expect_true(
  current = checkmate::testDataTable(internal_data_filtered$genes_to_pathways),
  info = "gene to pathway data returned"
)

expect_equal(
  current = sort(gois),
  target = sort(internal_data_filtered$represented_genes),
  info = "genes of interest represented"
)

## gw factory ------------------------------------------------------------------

gw_factory <- GeneWalkGenerator$new()

## adding data -----------------------------------------------------------------

gw_factory$add_ppi(source = "reactome")

expect_equal(
  current = gw_factory$ppi_sources,
  target = "reactome",
  info = "ppi addition behaves - reactome"
)

gw_factory$add_ppi(source = "signor")

expect_equal(
  current = gw_factory$ppi_sources,
  target = c("reactome", "signor"),
  info = "ppi addition behaves - reactome and signor"
)

gw_factory$add_pathways(source = "reactome")

expect_equal(
  current = gw_factory$pathway_sources,
  target = c("reactome"),
  info = "pathway addition - reactome works"
)

## erasing data ----------------------------------------------------------------

gw_factory$reset_choices()

expect_equal(
  current = gw_factory$pathway_sources,
  target = character(0),
  info = "pathway deletion behaves"
)

expect_equal(
  current = gw_factory$ppi_sources,
  target = character(0),
  info = "ppi deletion behaves"
)

## building the network --------------------------------------------------------

gw_factory$add_ppi(source = "reactome")$add_pathways(
  source = "reactome"
)$build(.verbose = FALSE)

full_network <- gw_factory$return_full_network_dt()

expect_equal(
  current = full_network,
  target = internal_data,
  info = "data from the factory is the same as manually created one"
)

expect_true(
  current = checkmate::testDataTable(full_network),
  info = "returning the full graph dt is behaving"
)

expect_true(
  current = checkmate::testNames(
    names(full_network),
    must.include = c("from", "to", "type")
  ),
  info = "correct columns on the graph dt"
)

## checking sub genes ----------------------------------------------------------

### generate a gw object -------------------------------------------------------

gw_obj <- gw_factory$create_for_genes(genes = gois)

### test the logic -------------------------------------------------------------

generated_graph <- get_graph_dt(gw_obj)

expect_true(
  current = all(gois %in% unique(c(generated_graph$from, generated_graph$to))),
  info = "genes of interest are represented"
)

expect_true(
  current = all(
    pathways_with_genes %in% unique(c(generated_graph$from, generated_graph$to))
  ),
  info = "reactome pathways are represented"
)
