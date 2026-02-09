# genewalk factory tests -------------------------------------------------------

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

expect_true(
  current = checkmate::testDataTable(full_network),
  info = "returning the full graph dt is behaving"
)

expect_true(
  current = checkmate::testNames(
    names(full_network),
    must.include = c("from", "to", "edge_type", "source")
  ),
  info = "correct columns on the graph dt"
)

## checking sub genes ----------------------------------------------------------

### generate a gw object -------------------------------------------------------

reactome_info <- get_reactome_info()
reactome_genes <- get_gene_to_reactome()

cell_cycle_reactome <- reactome_info[
  grepl("cell cycle", reactome_name),
  reactome_id
]

gois <- reactome_genes[to %in% cell_cycle_reactome, unique(from)]

pathways_with_genes <- reactome_genes[to %in% cell_cycle_reactome, unique(to)]

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
