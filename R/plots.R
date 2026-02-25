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

## similarities ----------------------------------------------------------------

#' Plot similarity distributions
#'
#' @description
#' Plots the cosine similarity distributions of the permuted embeddings
#' alongside the actual values as density histograms.
#'
#' @param object The `GeneWalk` class. Permutations and stats need to be
#' available in the object.
#'
#' @returns A ggplot2 object with the similarity distributions.
#'
#' @export
plot_similarities <- S7::new_generic(
  name = "plot_similarities",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method plot_similarities GeneWalk
#'
#' @import ggplot2
S7::method(plot_similarities, GeneWalk) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  if (length(object@permuted_embd) == 0) {
    warning(paste(
      "No permutations found in the object.",
      "Did you run permute_graph()?",
      "Returning NULL."
    ))
    return(NULL)
  }
  if (suppressWarnings(nrow(get_stats(object)) == 0)) {
    warning(paste(
      "No stats found in the object.",
      "Did you run calculate_genewalk_stats()?",
      "Returning NULL."
    ))
    return(NULL)
  }
  plot_dt <- purrr::imap(object@permuted_embd, \(x, i) {
    data.table::data.table(similarity = x, sample = sprintf("perm_%i", i))
  }) %>%
    append(list(data.table::data.table(
      similarity = object@stats$similarity,
      sample = "actual"
    ))) %>%
    data.table::rbindlist()

  p <- ggplot(data = plot_dt, mapping = aes(x = similarity)) +
    facet_wrap(~sample, scales = "free") +
    geom_histogram(aes(y = after_stat(density)), bins = 50) +
    xlab("Cosine similarity") +
    ylab("Density") +
    theme_bw()

  p
}

## scatter plot ----------------------------------------------------------------

#' Plot the gene walk results
#'
#' @description
#' Plots the GeneWalk results. Creates an interactive plotly plot to explore
#' the results interactively.
#'
#' @param object The `GeneWalk` class. The stats need to be available in the
#' object.
#' @param fdr_treshold Numeric. The FDR approach you wish to apply. GeneWalk
#' is quite stringent, so, setting this to lower values can make sense.
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
