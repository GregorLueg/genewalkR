# exported data ----------------------------------------------------------------

## myc genes -------------------------------------------------------------------

#' MYC target genes
#'
#' The MYC target genes from the MSigDB (Hallmarks V1).
#'
#' @format ## `myc_genes`
#' A data.table with 200 rows and 3 columns:
#' \describe{
#'   \item{ensembl_gene}{Ensembl identifier}
#'   \item{gene_symbol}{Gene Symbol}
#'   \item{gene_set}{Name of the gene set}
#' }
#'
#' @references Liberzon, et al., Cell Syst, 2015
"myc_genes"

# db utils ---------------------------------------------------------------------

## utils -----------------------------------------------------------------------

#' Transform string to factor columns
#'
#' @param table data.table. The data.table for which to transform the factor
#' columns into strings.
#'
#' @returns The data.table with the factor columns transformed to strings.
#'
#' @keywords internal
.string_cols_to_factors <- function(table) {
  # checks
  checkmate::assertDataTable(table)

  table <- table[,
    (names(table)[sapply(table, is.factor)]) := lapply(.SD, as.character),
    .SDcols = is.factor
  ]

  table
}

## downloads -------------------------------------------------------------------

#' Check if database file is a Git LFS pointer
#'
#' @param db_path String. The path to the DB
#'
#' @keywords internal
is_lfs_pointer <- function(db_path) {
  # checks
  checkmate::qassert(db_path, "S1")

  if (!file.exists(db_path)) {
    return(FALSE)
  }

  if (file.size(db_path) > 200) {
    return(FALSE)
  }

  first_line <- readLines(db_path, n = 1, warn = FALSE)
  grepl("^version https://git-lfs\\.github\\.com", first_line)
}

#' Download actual database from Git LFS
#'
#' @param repo String. Repository in format "user/repo" (default:
#' `"GregorLueg/genewalkR"`).
#' @param branch String. Branch name (default: `"main"`)
#'
#' @export
download_database <- function(repo = "GregorLueg/genewalkR", branch = "main") {
  # checks
  checkmate::qassert(repo, "S1")
  checkmate::qassert(branch, "S1")

  db_path <- system.file("extdata", "genewalk.duckdb", package = "genewalkR")

  if (!is_lfs_pointer(db_path)) {
    message("Database already downloaded")
    return(invisible(db_path))
  }

  # Construct LFS URL
  lfs_url <- sprintf(
    "https://media.githubusercontent.com/media/%s/%s/inst/extdata/genewalk.duckdb",
    repo,
    branch
  )

  message("Downloading database (150 MB)...")
  tryCatch(
    {
      download.file(lfs_url, db_path, mode = "wb", quiet = FALSE)
      message("Download complete")
      invisible(db_path)
    },
    error = function(e) {
      stop("Failed to download database: ", e$message, call. = FALSE)
    }
  )
}

## connection ------------------------------------------------------------------

#' Internal connector to the DB
#'
#' @param repo String. Repository in format "user/repo" (default:
#' `"GregorLueg/genewalkR"`).
#' @param branch String. Branch name (default: `"main"`)
#'
#' @returns The DuckDB connection to the internal data
#'
#' @export
get_db_connection <- function(repo = "GregorLueg/genewalkR", branch = "main") {
  # checks
  checkmate::qassert(repo, "S1")
  checkmate::qassert(branch, "S1")

  db_path <- system.file("extdata", "genewalk.duckdb", package = "genewalkR")

  # check if the DB is only a pointer - if yes, download
  if (is_lfs_pointer(db_path)) {
    message("Downloading database (150 MB)...")
    lfs_url <- sprintf(
      "https://media.githubusercontent.com/media/%s/%s/inst/extdata/genewalk.duckdb",
      repo,
      branch
    )
    tryCatch(
      {
        download.file(lfs_url, db_path, mode = "wb", quiet = FALSE)
        message("Download complete")
        invisible(db_path)
      },
      error = function(e) {
        stop("Failed to download database: ", e$message, call. = FALSE)
      }
    )
  }

  DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE)
}

## getters ---------------------------------------------------------------------

### main tables ----------------------------------------------------------------

#' Get the gene main information
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the gene information
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_gene_info <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM gene_main"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

#' Get the Reactome main information
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Reactome information
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_reactome_info <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM reactome_main"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

#' Get the Gene Ontology main information
#'
#' @inherit gene_ontology_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Gene Ontology information
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_gene_ontology_info <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM gene_ontology_main"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

### pathway informations -------------------------------------------------------

#' Get the Gene Ontology to gene connections
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Gene Ontology to gene edges
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_gene_to_go <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM gene_ontology_genes"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}


#' Get the Reactome to gene connections
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Reactome to gene edges
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_gene_to_reactome <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM reactome_genes"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

### pathway hierarchy ----------------------------------------------------------

#' Get the Reactome pathway hierarchy
#'
#' @inherit opentargets_source details
#'
#' @param relationship String vector. Defines which relationships to return.
#' Multiple allowed and choices from
#' `c("child_of", "descendant_of", "parent_of", "ancestor_of")`
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Reactome to Reactome edges with edge type info.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_reactome_hierarchy <- function(
  relationship = c("child_of", "descendant_of", "parent_of", "ancestor_of"),
  head_only = FALSE
) {
  relationship <- match.arg(relationship, several.ok = TRUE)
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  relationship_str <- paste0("'", relationship, "'", collapse = ", ")
  query <- sprintf(
    "SELECT * FROM reactome_ontology WHERE relationship IN (%s)",
    relationship_str
  )

  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

#' Get the Gene Ontology
#'
#' @inherit gene_ontology_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the GO to GO edges with edge type info.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_gene_ontology_hierarchy <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM gene_ontology_dag"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

### interactions ---------------------------------------------------------------

#' Get the STRING network
#'
#' @inherit opentargets_source details
#'
#' @param threshold Optional numeric. If you wish to filter on the STRING score.
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the String network data.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_interactions_string <- function(threshold = NULL, head_only = FALSE) {
  checkmate::qassert(threshold, c("0", "N1[0, 1]"))
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- if (is.null(threshold)) {
    "SELECT * FROM interactions_string"
  } else {
    sprintf(
      "SELECT * FROM interactions_string WHERE string_score >= %f",
      threshold
    )
  }

  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

#' Get the SIGNOR network
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the SIGNOR network data.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_interactions_signor <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM interactions_signor"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}


#' Get the Reactome network
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Reactome network data.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_interactions_reactome <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM interactions_reactome"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}

#' Get the Intact network
#'
#' @inherit opentargets_source details
#'
#' @param threshold Optional numeric. If you wish to filter on the STRING score.
#' @param physical_interactions Boolean. Shall only physical interactions be
#' returned.
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Intact network data.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
get_interactions_intact <- function(
  threshold = NULL,
  physical_interactions = FALSE,
  head_only = FALSE
) {
  # checks
  checkmate::qassert(threshold, c("0", "N1[0, 1]"))
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- if (is.null(threshold) & !physical_interactions) {
    "SELECT * FROM interactions_intact"
  } else if (is.null(threshold) & physical_interactions) {
    "SELECT * FROM interactions_intact WHERE physical_association"
  } else if (!physical_interactions) {
    sprintf(
      "SELECT * FROM interactions_intact WHERE intact_score >= %f",
      threshold
    )
  } else {
    sprintf(
      "SELECT * FROM interactions_intact WHERE intact_score >= %f AND physical_association",
      threshold
    )
  }

  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}


#' Get the Combined network
#'
#' @description This returns a combined network akin to the workflow from
#' Barrio-Hernandez et al. based on high quality edges from STRING, physical
#' interactions from Intact and the full sets from SIGNOR and Reactome
#'
#' @inherit opentargets_source details
#'
#' @param head_only Logical. If TRUE, returns only the first 5 rows.
#'
#' @returns data.table with the Reactome network data.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @references Barrio-Hernandez et al., Nat Gen, 2023
get_interactions_combined <- function(head_only = FALSE) {
  checkmate::assert_flag(head_only)

  con <- get_db_connection()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- "SELECT * FROM interactions_combined"
  if (head_only) {
    query <- paste(query, "LIMIT 5")
  }

  table <- DBI::dbGetQuery(conn = con, statement = query) %>%
    data.table::setDT() %>%
    .string_cols_to_factors()

  return(table)
}
