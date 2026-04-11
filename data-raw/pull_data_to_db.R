# libraries, paths, helpers ----------------------------------------------------

library(here)
library(data.table)
library(duckdb)
library(magrittr)

devtools::load_all()

# functions --------------------------------------------------------------------

## parameters ------------------------------------------------------------------

download_all_data <- TRUE
dir_data <- "~/Desktop/genewalk_data/"

## constants -------------------------------------------------------------------

relationship_translation <- setNames(
  c(
    "child_of",
    "parent_of",
    "ancestor_of",
    "descendant_of"
  ),
  c("parents", "children", "descendants", "ancestors")
)

# process data -----------------------------------------------------------------

## download the data -----------------------------------------------------------

options(timeout = 100000)

if (download_all_data) {
  download_data(target_path = dir_data)
  download_pathway_commons(target_path = dir_data)
}

list.files(dir_data)

## individual tables -----------------------------------------------------------

### genes ----------------------------------------------------------------------

gene_main <- ingest_gene_data(dir_data = dir_data)

### reactome -------------------------------------------------------------------

reactome_main <- ingest_reactome_data(dir_data = dir_data)
reactome_genes <- ingest_reactome_genes(dir_data = dir_data)
reactome_ontology <- ingest_reactome_ontology(dir_data = dir_data)

### gene ontology --------------------------------------------------------------

go_data_parsed <- parse_gene_ontology(dir_data = dir_data)

go_main <- ingest_go_data(obo_list = go_data_parsed)
go_genes <- ingest_gene_ontology_genes(dir_data = dir_data)
go_ontology <- ingest_gene_ontology_relations(obo_list = go_data_parsed)

### interaction ----------------------------------------------------------------

interactions_string <- ingest_string_interactions(dir_data = dir_data)
interactions_signor <- ingest_signor_interactions(dir_data = dir_data)
interactions_reactome <- ingest_reactome_interactions(dir_data = dir_data)
interactions_intact <- ingest_intact_interaction(dir_data = dir_data)
interactions_pc <- ingest_pc_interactions(dir_data = dir_data)
interactions_combined <- generate_combined_network(
  interactions_string = interactions_string,
  interactions_signor = interactions_signor,
  interactions_reactome = interactions_reactome,
  interactions_intact = interactions_intact,
  interactions_pc = interactions_pc
)

# generate the internal db -----------------------------------------------------

## provide a path and clean up the db ------------------------------------------

dir.create(
  file.path(dir_data, "duckdb"),
  recursive = TRUE,
  showWarnings = FALSE
)

db_path <- file.path(dir_data, "duckdb/genewalk.duckdb")

if (checkmate::testFileExists(db_path)) {
  unlink(db_path, force = TRUE)
}

## populate the db -------------------------------------------------------------

table_list <- list(
  # main tables
  gene_main = gene_main,
  reactome_main = reactome_main,
  gene_ontology_main = go_main,
  # pathways to genes
  reactome_genes = reactome_genes,
  gene_ontology_genes = go_genes,
  # ontology
  reactome_ontology = reactome_ontology,
  gene_ontology_dag = go_ontology,
  # interactions
  interactions_string = interactions_string,
  interactions_signor = interactions_signor,
  interactions_reactome = interactions_reactome,
  interactions_intact = interactions_intact,
  interactions_pc = interactions_pc,
  interactions_combined = interactions_combined
)

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = FALSE)

# should reduce file size
DBI::dbExecute(con, "SET default_null_order='nulls_first'")

for (i in seq_along(table_list)) {
  table_name <- names(table_list)[i]
  message(sprintf("Ingesting table %s into the internal DB.", table_name))
  DBI::dbWriteTable(conn = con, name = table_name, table_list[[i]])
  no_rows_written <- DBI::dbGetQuery(
    conn = con,
    statement = sprintf("SELECT COUNT(*) AS nrow FROM %s", table_name)
  )
  message(sprintf(
    " Wrote %s rows.",
    format(no_rows_written$nrow, big.mark = "_")
  ))
}

# vacuum helps with file size
DBI::dbExecute(con, "CHECKPOINT")
DBI::dbExecute(con, "VACUUM")

DBI::dbDisconnect(con, shutdown = TRUE)

db_path_zip = sprintf("%s.zip", db_path)

zip(path.expand(db_path_zip), path.expand(db_path))
