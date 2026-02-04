# libraries, paths, helpers ----------------------------------------------------

library(here)
library(data.table)
library(duckdb)
library(magrittr)

# functions --------------------------------------------------------------------

## parameters ------------------------------------------------------------------

download_all_data <- FALSE
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

## downloads -------------------------------------------------------------------

#' Helper function to download necessary data
#'
#' @param target_path String. Where do download the data
#'
#' @returns Downloads the necessary files to target_path.
download_data <- function(target_path) {
  checkmate::qassert(target_path, "S1")

  cmd <- sprintf(
    "bash %s/data-raw/scripts/downloads.sh -t %s",
    here::here(),
    target_path
  )

  system(cmd)
}

## parsing ---------------------------------------------------------------------

#' Parse OBO file to list of data.tables
#'
#' @description
#' This is a generalised helper function to read in .obo files (ontology files)
#' and generate some useful format out of it.
#'
#' @param file_path String. Path to the .obo file
#' @param verbose Boolean. Controls verbosity of the function.
#'
#' @return A list the following elements:
#' \itemize{
#'  \item names - data.table with name and ID of the term.
#'  \item alt_ids - data.table with alternative ids.
#'  \item namespace - data.table with the found name-spaces.
#'  \item term_subsets - data.table with the found subsets.
#'  \item term_obsolete - character vector of terms that were marked obsolete.
#'  \item term_relations - data.table with the ontological relationships found
#'  in the file.
#'  \item synonyms - data.table with synonyms provided in the .obo file.
#'  \item externalRef - data.table with external references
#' }
obo_parser <- function(file_path, verbose = TRUE) {
  # Checks
  checkmate::assertFileExists(file_path)
  checkmate::qassert(verbose, "B1")

  # Load in the file
  con <- file(file_path, "r")
  on.exit(close(con))

  # Raw
  raw <- readLines(con) %>%
    stringr::str_split_fixed(":", n = 2) %>%
    `colnames<-`(c("key", "value")) %>%
    data.table::as.data.table() %>%
    .[, `:=`(value = trimws(value), row = .I)]

  # Add a column with term IDs and remove unneeded stuff
  terms <- raw[key == "id"]
  raw[,
    term := cut(row, c(terms$row - 1, Inf), labels = terms$value, right = FALSE)
  ]
  raw[, row := NULL]
  externals <- raw[key == "namespace" & value == "external", term]
  raw <- raw[!is.na(term) & !term %in% externals]

  # Get the names, namespace, obsolete terms and subsets
  term_name <- raw[key == "name"][, .(ID = term, name = value)]
  term_namespace <- raw[key == "namespace"][,
    .(IDs = list(unique(as.character(term)))),
    .(name_space = value)
  ]
  obsolete_terms <- raw[
    key == "is_obsolete" & value == "true",
    as.character(term)
  ]
  term_subset <- raw[key == "subset"][,
    .(subset = list(unique(as.character(term)))),
    .(value)
  ]

  # Get the relationships
  if (verbose) {
    message("Extracting the ontology relationships.")
  }
  term_relations <- raw[key == "relationship", value] %>%
    stringr::str_split(., "!", n = 2) %>%
    purrr::map_chr(1) %>%
    trimws() %>%
    data.table::data.table(original = .) %>%
    .[,
      c("relationship", "ID") := data.table::tstrsplit(
        original,
        " ",
        fixed = T,
        keep = c(1:2)
      )
    ] %>%
    .[, `:=`(original = NULL, term = raw[key == "relationship", term])] %>%
    .[, .(from = term, to = ID, relationship)]
  is_a_relations <- raw[key == "is_a", value] %>%
    stringr::str_split(., "!", n = 2) %>%
    purrr::map_chr(1) %>%
    trimws() %>%
    data.table::data.table(to = .) %>%
    .[, `:=`(from = raw[key == "is_a", term], relationship = "is_a")] %>%
    .[, .(from, to, relationship)]
  all_relations <- rbind(term_relations, is_a_relations)

  # Synonyms
  if (verbose) {
    message(
      "Extracting the synonyms, alternative identifiers and external refs."
    )
  }
  synonyms <- raw[key == "synonym"][,
    c("synonym", "type") := list(
      stringr::str_extract(value, '"(.*?)"') %>% gsub('"', "", .),
      stringr::str_replace(value, '.*"', "") %>%
        trimws() %>%
        stringr::str_split(., " ", n = 2) %>%
        purrr::map_chr(1)
    )
  ] %>%
    .[, .(term, synonym, type)] %>%
    unique()

  # External references
  external_ref <- raw[key == "xref"] %>%
    .[,
      externalRef := data.table::tstrsplit(value, " ", fixed = T, keep = 1)
    ] %>%
    .[, .(term, externalRef)]

  # And alternative IDs
  alternative_ID = raw[key == "alt_id"] %>%
    .[, .(term, altID = value)]

  result <- list(
    names = term_name,
    alt_ids = alternative_ID,
    namespace = term_namespace,
    term_subsets = term_subset,
    term_obsolete = obsolete_terms,
    term_relations = all_relations,
    synonyms = synonyms,
    external_ref = external_ref
  )

  result
}

#' Parses the Gene Ontology
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A named list with the different sub tables and data extracted from
#' the .obo file
parse_gene_ontology <- function(dir_data) {
  obo_parser(file_path = file.path(dir_data, "obo_files/go.obo"))
}

## main tables -----------------------------------------------------------------

#' Generates a genes main table
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A data.table with the genes to gene information (Ensembl
#' identifiers, symbols and biotype)
get_gene_data <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  # connection
  con <- dbConnect(duckdb(), dbdir = ":memory:")

  gene_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT DISTINCT id AS ensembl_id, approvedSymbol AS symbol, biotype
      FROM read_parquet('%s')",
      file.path(dir_data, "targets/*.parquet")
    )
  ) %>%
    setDT()

  gene_data
}

#' Generates a Reactome main table
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A data.table with the Reactome main table
get_reactome_data <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  # connection
  con <- dbConnect(duckdb(), dbdir = ":memory:")

  reactome_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT id AS reactome_id, label AS reactome_name 
      FROM read_parquet('%s')",
      file.path(dir_data, "reactome/*.parquet")
    )
  ) %>%
    setDT()

  reactome_data
}

#' Generates the Gene Ontology table
#'
#' @param obo_list Named list. Output of `parse_gene_ontology()`
#'
#' @returns A data.table with the Gene ontology term ids, names and their
#' namespace
get_go_data <- function(obo_list) {
  # checks
  checkmate::assertList(obo_list)
  checkmate::checkNames(names(obo_list), must.include = c("names", "namespace"))

  go_main <- merge(
    obo_list$names,
    obo_list$namespace[, .(ID = unlist(IDs)), by = name_space],
    by = "ID"
  )
  setnames(go_main, old = c("ID", "name"), new = c("go_id", "go_name"))

  return(go_main)
}

## relationships ---------------------------------------------------------------

### pathway to genes -----------------------------------------------------------

#' Generates a Reactome to gene table
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A data.table with the genes to Reactome associations
get_reactome_genes <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  # connection
  con <- dbConnect(duckdb(), dbdir = ":memory:")

  # get data
  reactome_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT id AS from, pathway.pathwayId AS to
      FROM (
        SELECT id, approvedSymbol, UNNEST(pathways) AS pathway
        FROM read_parquet('%s')
        WHERE pathways IS NOT NULL
      )",
      file.path(dir_data, "targets/*.parquet")
    )
  ) %>%
    setDT()

  reactome_data
}

#' Generates a Gene Ontology to gene table
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A data.table with the Genes to Gene Ontology term associations
get_gene_ontology_genes <- function(dir_data) {
  checkmate::assertDirectoryExists(dir_data)

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  go_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT DISTINCT id AS from, go.id AS to
      FROM (
        SELECT id, approvedSymbol, UNNEST(go) AS go
        FROM read_parquet('%s')
        WHERE go IS NOT NULL
      )",
      file.path(dir_data, "targets/*.parquet")
    )
  )

  go_data
}

### ontologies -----------------------------------------------------------------

#' Generates a Reactome term ontology table
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns A data.table with the Reactome ontology (child of, etc.)
get_reactome_ontology <- function(dir_data) {
  checkmate::assertDirectoryExists(dir_data)

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  reactome_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT id AS from, UNNEST(value) AS to, variable AS relationship
      FROM (
        SELECT id, variable, value
        FROM read_parquet('%s')
        UNPIVOT (
          value FOR variable IN (parents, children, ancestors, descendants)
        )
        WHERE value IS NOT NULL
      )
      WHERE value IS NOT NULL",
      file.path(dir_data, "reactome/*.parquet")
    )
  ) %>%
    setDT() %>%
    .[, relationship := relationship_translation[relationship]]

  reactome_data
}

#' Generates the Gene Ontology table
#'
#' @param obo_list Named list. Output of `parse_gene_ontology()`
#'
#' @returns A data.table with the Gene Ontology relationships
get_gene_ontology_relations <- function(obo_list) {
  # checks
  checkmate::assertList(obo_list)
  checkmate::checkNames(names(obo_list), must.include = c("names", "namespace"))

  go_data_parsed$term_relations
}

### interactions ---------------------------------------------------------------

#' Pull out the string data from OpenTargets
#'
#' @param dir_data String. The directory in which you downloaded the data
#' @param threshold Numeric. The string threshold from which to consider a
#' STRING connection high quality.
#'
#' @returns The data.table with the STRING network data.
get_string_interactions <- function(dir_data, threshold = 0.85) {
  # checks
  checkmate::assertDirectoryExists(dir_data)
  checkmate::qassert(threshold, "N1[0, 1]")

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  # get the string stuff
  string_data <- DBI::dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT DISTINCT
      targetA AS from,
      targetB AS to,
      interactionScore AS string_score,
      (interactionScore >= %f) AS string_high_quality
      FROM (
        SELECT 
          targetA, 
          targetB, 
          interactionScore,
          UNNEST(interactionResources, recursive := true)
        FROM read_parquet('%s')
      )
      WHERE sourceDatabase = 'string'
        AND targetA IS NOT NULL
        AND targetB IS NOT NULL",
      threshold,
      file.path(dir_data, "interactions/*.parquet")
    )
  ) %>%
    setDT()

  # remove non-human stuff
  string_data <- string_data[
    (grepl("ENSG", to) & grepl("ENSG", from)) &
      !(grepl("ENSGGOG", to) | grepl("ENSGGOG", from))
  ]

  unique(string_data)
}

#' Pull out the Reactome interaction data from OpenTargets
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns The data.table with the Reactome interaction
get_reactome_interactions <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  reactome_interactions_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT DISTINCT
        targetA AS from,
        targetB AS to,
        interactionTypeShortName AS interaction_type
      FROM (
        SELECT 
          targetA, 
          targetB, 
          interactionTypeShortName,
          UNNEST(interactionResources, recursive := true)
        FROM read_parquet('%s')
      )
      WHERE sourceDatabase = 'reactome'
        AND targetA IS NOT NULL
        AND targetB IS NOT NULL",
      file.path(dir_data, "interactions/*.parquet")
    )
  ) %>%
    setDT()

  # remove non-human stuff
  reactome_interactions_data <- reactome_interactions_data[
    (grepl("ENSG", to) & grepl("ENSG", from)) &
      !(grepl("ENSGGOG", to) | grepl("ENSGGOG", from))
  ]

  unique(reactome_interactions_data)
}

#' Pull out the SIGNOR interaction data from OpenTargets
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns The data.table with the SIGNOR interaction
get_signor_interactions <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  signor_interactions_data <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT DISTINCT
        targetA AS from,
        targetB AS to,
        interactionTypeShortName AS interaction_type
      FROM (
        SELECT 
          targetA, 
          targetB, 
          interactionTypeShortName,
          UNNEST(interactionResources, recursive := true)
        FROM read_parquet('%s')
      )
      WHERE sourceDatabase = 'signor'
        AND targetA IS NOT NULL
        AND targetB IS NOT NULL",
      file.path(dir_data, "interactions/*.parquet")
    )
  ) %>%
    setDT()

  # remove non-human stuff
  signor_interactions_data <- signor_interactions_data[
    (grepl("ENSG", to) & grepl("ENSG", from)) &
      !(grepl("ENSGGOG", to) | grepl("ENSGGOG", from))
  ]

  unique(signor_interactions_data)
}

#' Pull out the Intact interaction data from OpenTargets
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns The data.table with the Intact interaction
get_intact_interaction <- function(dir_data) {
  # checks
  checkmate::assertDirectoryExists(dir_data)

  con <- dbConnect(duckdb(), dbdir = ":memory:")

  intact_interactions <- dbGetQuery(
    conn = con,
    statement = sprintf(
      "SELECT
        targetA AS from,
        targetB AS to,
        interactionScore AS intact_score,
        CONTAINS(STRING_AGG(DISTINCT interactionTypeShortName, '; '), 'physical association') AS physical_association
      FROM (
        SELECT 
          targetA, 
          targetB, 
          interactionScore,
          interactionTypeShortName,
          UNNEST(interactionResources, recursive := true)
        FROM read_parquet('%s')
      )
      WHERE sourceDatabase = 'intact'
        AND targetA IS NOT NULL
        AND targetB IS NOT NULL
      GROUP BY targetA, targetB, interactionScore",
      file.path(dir_data, "interactions/*.parquet")
    )
  ) %>%
    setDT()

  # remove non-human stuff
  intact_interactions <- intact_interactions[
    (grepl("ENSG", to) & grepl("ENSG", from)) &
      !(grepl("ENSGGOG", to) | grepl("ENSGGOG", from))
  ]

  unique(intact_interactions)
}

#' Combine
#'
#' @param dir_data String. The directory in which you downloaded the data
#'
#' @returns The data.table with the Intact interaction
generate_combined_network <- function(
  interactions_string,
  interactions_signor,
  interactions_reactome,
  interactions_intact
) {
  # checks
  checkmate::assertDataTable(interactions_string)
  checkmate::assertDataTable(interactions_signor)
  checkmate::assertDataTable(interactions_reactome)
  checkmate::assertDataTable(interactions_intact)

  # generate a combined network
  combined_network <- rbindlist(
    list(
      copy(interactions_string)[(string_high_quality), c("from", "to")][,
        `:=`(source = "string")
      ],
      copy(interactions_signor)[, c("from", "to")][, source := "signor"],
      copy(interactions_reactome)[, c("from", "to")][, source := "reactome"],
      copy(interactions_intact)[(physical_association), c("from", "to")][,
        `:=`(source = "intact")
      ]
    )
  )[, in_source := TRUE] %>%
    dcast(
      .,
      from + to ~ source,
      value.var = "in_source",
      fun.aggregate = sum
    )

  source_cols <- setdiff(names(combined_network), c("from", "to"))
  combined_network[,
    (source_cols) := lapply(.SD, \(x) x > 0),
    .SDcols = source_cols
  ]

  return(combined_network)
}

# process data -----------------------------------------------------------------

## download the data -----------------------------------------------------------

if (download_all_data) {
  download_data(target_path = dir_data)
}


## individual tables -----------------------------------------------------------

### genes ----------------------------------------------------------------------

gene_main <- get_gene_data(dir_data = dir_data)

### reactome -------------------------------------------------------------------

reactome_main <- get_reactome_data(dir_data = dir_data)
reactome_genes <- get_reactome_genes(dir_data = dir_data)
reactome_ontology <- get_reactome_ontology(dir_data = dir_data)

### gene ontology --------------------------------------------------------------

go_data_parsed <- parse_gene_ontology(dir_data = dir_data)

go_main <- get_go_data(obo_list = go_data_parsed)
go_genes <- get_gene_ontology_genes(dir_data = dir_data)
go_ontology <- get_gene_ontology_relations(obo_list = go_data_parsed)

### interaction ----------------------------------------------------------------

interactions_string <- get_string_interactions(dir_data = dir_data)
interactions_signor <- get_signor_interactions(dir_data = dir_data)
interactions_reactome <- get_reactome_interactions(dir_data = dir_data)
interactions_intact <- get_intact_interaction(dir_data = dir_data)
interactions_combined <- generate_combined_network(
  interactions_string = interactions_string,
  interactions_signor = interactions_signor,
  interactions_reactome = interactions_reactome,
  interactions_intact = interactions_intact
)

# generate the internal db -----------------------------------------------------

## provide a path and clean up the db ------------------------------------------

db_path <- file.path(here::here(), "inst/extdata/genewalk.duckdb")

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

DBI::dbGetQuery(
  con,
  "
  SELECT 
    table_name,
    column_name, 
    compression
  FROM duckdb_columns()
"
)

DBI::dbDisconnect(con, shutdown = TRUE)
