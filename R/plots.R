# plotting functions -----------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Helper function to get the node degree from the PPI
#'
#' @param graph_dt data.table. The graph info for which to get the node degree
#' of the genes.
#'
#' @returns data.table with the gene and node degree in the network
#'
#' @keywords internal
.get_node_degree <- function(graph_dt) {
  # checks
  assertGeneWalkDataTable(graph_dt)

  node_degree <- igraph::graph_from_data_frame(
    graph_dt[type == "interaction"],
    directed = FALSE
  ) %>%
    igraph::degree()

  node_dt <- data.table::data.table(
    gene = names(node_degree),
    degree = node_degree
  )

  node_dt
}

# plot the results -------------------------------------------------------------

#' Plot the gene walk results
#'
#' @description
#' Plots the GeneWalk results. Creates an interactive plotly plot to explore
#' the results interactively.
#'
#' @param object The `GeneWalk` class. The stats need to be available in the
#' object.
#'
#' @returns ggplot2 object with the plot for the gene
#'
#' @export
plot_gw_results <- S7::new_generic(
  name = "plot_gw_results",
  dispatch_args = "object",
  fun = function(
    object,
    fdr_treshold = 0.25
  ) {
    S7::S7_dispatch()
  }
)

#' @method plot_gw_results GeneWalk
#'
#' @import ggplot2
S7::method(plot_gw_results, GeneWalk) <- function(
  object,
  fdr_treshold = 0.25
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  checkmate::qassert(fdr_treshold, "N1[0, 1]")

  # early return
  if (suppressWarnings(nrow(get_stats(object)) == 0)) {
    warning(paste(
      "No stats found in the object.",
      "Did you run calculate_genewalk_stats()?",
      "Returning NULL."
    ))
    return(NULL)
  }

  results <- get_stats(object)

  graph_dt <- get_graph_dt(object)

  node_dt <- .get_node_degree(graph_dt)

  plot_dt <- results[,
    .(
      ratio = sum(avg_gene_fdr <= fdr_treshold) / length(pathway),
      pathway_connections = length(pathway)
    ),
    .(gene)
  ]

  plot_dt <- node_dt[plot_dt, on = "gene"][,
    degree := fifelse(is.na(degree), 0, degree)
  ]

  p <- ggplot(data = plot_dt, mapping = aes(x = degree + 1, y = ratio)) +
    suppressWarnings(geom_point(
      mapping = aes(fill = ratio, text = gene, size = pathway_connections),
      shape = 21,
      alpha = 0.7
    )) +
    scale_x_log10() +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis") +
    ylim(0, 1) +
    xlab("log(degree + 1)") +
    ylab("Ratio enrichment") +
    labs(fill = "Ratio", size = "Pathway\nconnections") +
    theme_bw() +
    ggtitle("Gene Walk Results")

  plotly::ggplotly(p, tooltip = c("text", "x", "y", "size"))
}
