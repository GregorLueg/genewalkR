# Prepare internal network datasets
library(zeallot)

# Helper functions to fetch specific networks from DuckLake --------------------

fetch_tables_from_ducklake = function(ducklake_con, table_names) {
  lapply(table_names, function(tbl_name) {
    cosyne.ducklake::ducklake_fetch_table(
      ducklake_con = ducklake_con,
      table_name = tbl_name
    ) |>
      data.table::setDT()
  })
}

get_g2g_interactions_nw = function(ducklake_con) {
  `%<-%` <- zeallot::`%<-%`
  c(interactions, gene_annotation) %<-%
    fetch_tables_from_ducklake(
      ducklake_con,
      c(
        'opentargets_interactions_combined_network',
        'gex_feature_info'
      )
    )

  graph = interactions[,
    `:=`(
      edge_type = 'g2g_interacts',
      source = 'opentargets'
    )
  ]

  # Not all ensembl ids in the interactions are in the gene annotation
  gene_annotation = gene_annotation[
    .(gene_id = unique(c(graph$from, graph$to))),
    on = 'gene_id',
    nomatch = NA
  ]
  checkmate::assert_true(all(
    c(graph$from, graph$to) %in% gene_annotation$gene_id
  ))
  node_ids = gene_annotation[
    gene_id %in% c(graph$from, graph$to),
    .(id = gene_id, name = gene_name, node_type = 'gene')
  ][,
    unique(.SD)
  ]

  return(list(graph = graph, node_id = node_ids))
}

get_reactome_nw = function(ducklake_con) {
  `%<-%` <- zeallot::`%<-%`
  c(
    reactome_relationships,
    reactome_to_genes,
    reactome_main
  ) %<-%
    fetch_tables_from_ducklake(
      ducklake_con = ducklake_con,
      c('reactome_relationships', 'reactome_to_genes', 'reactome_main')
    )

  graph = data.table::rbindlist(list(
    reactome_relationships[,
      .(from, to, edge_type = relationship, source = 'reactome_relationships')
    ],
    reactome_to_genes[,
      .(from, to, edge_type = relationship, source = 'reactome_to_genes')
    ]
  ))

  reactome_main = reactome_main[,
    .(id = pathway_id, name = pathway_name, node_type = 'pathway')
  ]

  list(graph = graph, node_id = reactome_main)
}

get_gene_ontology_nw = function(
  ducklake_con
) {
  `%<-%` <- zeallot::`%<-%`
  c(
    go_relationships,
    go_to_genes,
    go_main
  ) %<-%
    fetch_tables_from_ducklake(
      ducklake_con = ducklake_con,
      c(
        'gene_ontology_relationships',
        'gene_ontology_to_genes',
        'gene_ontology_main'
      )
    )

  graph = data.table::rbindlist(list(
    go_relationships[,
      .(
        from,
        to,
        edge_type = relationship,
        source = 'gene_ontology_relationships'
      )
    ],
    go_to_genes[,
      .(from, to, edge_type = relationship, source = 'gene_ontology_to_genes')
    ]
  ))

  go_main = go_main[,
    .(id = go_id, name = go_name, namespace = namespace, node_type = 'go_term')
  ]

  list(graph = graph, node_id = go_main)
}

# Main function to generate exported dataset -----------------------------------
generate_exported_networks = function(
  resources = c(
    'g2g_interactions',
    'reactome',
    'gene_ontology'
  )
) {
  checkmate::assert_true(
    'g2g_interactions' %in% resources && length(resources) > 1
  )

  # Fetch data from DuckLake
  ducklake_con = cosyne.ducklake::ducklake_connect()
  if ('g2g_interactions' %in% resources) {
    g2g_interactions = get_g2g_interactions_nw(ducklake_con)
    usethis::use_data(g2g_interactions, overwrite = TRUE)
  }
  if ('reactome' %in% resources) {
    reactome = get_reactome_nw(ducklake_con)
    usethis::use_data(reactome, overwrite = TRUE)
  }
  if ('gene_ontology' %in% resources) {
    gene_ontology = get_gene_ontology_nw(ducklake_con)
    usethis::use_data(gene_ontology, overwrite = TRUE)
  }

  message("Exported datasets generated successfully.")
}
