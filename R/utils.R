#' Generate default graph_dt from the network datasets provided by the package
#'
#' @param resources Character vector of networks to include.
#'   Choices: "g2g_interactions", "reactome", "gene_ontology".
#' @param g2g_source Which G2G interaction sources to include.
#' @param go_type Gene Ontology namespace to include.
#'
#' @return A list with:
#' \describe{
#'   \item{graph_dt}{A `data.table` suitable for `genewalkR_class`}
#'   \item{graph_gene_params}{Metadata describing how the graph was generated}
#' }
#' @seealso [g2g_interactions], [reactome], [gene_ontology] for dataset details.
#'
#' @export
get_default_graph_dt <- function(
  resources = c("g2g_interactions", "reactome", "gene_ontology"),
  g2g_source = c('all', 'intact', 'string', 'signor', 'reactome'),
  go_type = c("all", "BP", "MF", "CC")
) {
  # Arguments checks
  checkmate::qassert(resources, 'S+')
  available_resources <- c("g2g_interactions", "reactome", "gene_ontology")
  checkmate::assert_subset(resources, available_resources)
  checkmate::assert_true('g2g_interactions' %in% resources)

  g2g_source <- match.arg(g2g_source)
  go_type <- match.arg(go_type)

  # Load networks
  networks <- lapply(resources, function(res) {
    get(res, envir = asNamespace("genewalkR"))
  })
  names(networks) <- resources

  # Filter G2G interactions by source
  if (g2g_source != 'all') {
    g2g <- networks[["g2g_interactions"]]
    source_cols = paste0('from_', g2g_source)
    checkmate::assert_subset(source_cols, colnames(g2g$graph))

    g2g$graph <- g2g$graph[
      g2g$graph[, Reduce(`|`, lapply(.SD, isTRUE)), .SDcols = source_cols]
    ]

    g2g$node_id <- g2g$node_id[
      id %in% unique(c(g2g$graph$from, g2g$graph$to))
    ]

    networks[["g2g_interactions"]] <- g2g
  }

  if ('gene_ontology' %in% resources && go_type != 'all') {
    go <- networks[["gene_ontology"]]

    namespace_map <- c(
      BP = "biological_process",
      MF = "molecular_function",
      CC = "cellular_component"
    )
    ns <- namespace_map[[go_type]]

    go$node_id <- go$node_id[namespace == ns]
    go_ids <- go$node_id$id

    go$graph <- go$graph[
      (source == "gene_ontology_relationships" &
        from %in% go_ids &
        to %in% go_ids) |
        (source == "gene_ontology_to_genes" & to %in% go_ids)
    ]

    networks[["gene_ontology"]] <- go
  }

  # Combine networks
  graph_dt <- data.table::rbindlist(
    lapply(networks, `[[`, "graph"),
    fill = TRUE
  )

  node_id <- data.table::rbindlist(
    lapply(networks, `[[`, "node_id"),
    fill = TRUE
  )

  graph_gene_params <- list(
    resources = resources,
    node_id = node_id,
    g2g_source = g2g_source
  )

  if ('gene_ontology' %in% resources) {
    graph_gene_params$go_type = go_type
  }

  list(
    graph_dt = graph_dt,
    graph_gene_params = graph_gene_params
  )
}

#' Trim a heterogeneous graph by gene set
#'
#' @param graph_dt data.table with at least `from` and `to`
#' @param node_id data.table with columns `id` and `type`
#' @param genes_to_keep character vector of gene IDs
#' @param gene_type value in `node_id$node_type` identifying genes
#'
#' @return list with trimmed `graph_dt` and `node_id`
#' @export
trim_graph <- function(
  graph_dt,
  node_id,
  genes_to_keep,
  gene_type = "gene"
) {
  checkmate::assertDataTable(graph_dt)
  checkmate::assertDataTable(node_id)
  checkmate::assertNames(names(graph_dt), must.include = c("from", "to"))
  checkmate::assertNames(names(node_id), must.include = c("id", "node_type"))
  checkmate::qassert(genes_to_keep, 'S+')

  # annotate node types onto edges
  graph_annotated <- merge(
    graph_dt,
    node_id[, .(from = id, from_type = node_type)],
    by = "from",
    all.x = TRUE
  )

  graph_annotated <- merge(
    graph_annotated,
    node_id[, .(to = id, to_type = node_type)],
    by = "to",
    all.x = TRUE
  )

  # Trim edges
  graph_trimmed <- graph_annotated[
    # gene–gene edges
    (from_type == gene_type &
      to_type == gene_type &
      from %in% genes_to_keep &
      to %in% genes_to_keep) |

      # gene–non-gene edges
      (from_type == gene_type &
        to_type != gene_type &
        from %in% genes_to_keep) |
      (to_type == gene_type & from_type != gene_type & to %in% genes_to_keep) |

      # non-gene–non-gene edges (ontology/pathway structure)
      (from_type != gene_type & to_type != gene_type)
  ][, c("from_type", "to_type") := NULL]

  # trim node table
  node_id_trimmed <- node_id[
    id %in% unique(c(graph_trimmed$from, graph_trimmed$to))
  ]

  list(
    graph_dt = graph_trimmed,
    node_id = node_id_trimmed
  )
}
