# plotting functions -----------------------------------------------------------

## individual genes ------------------------------------------------------------

#' Generate a plot for an individual gene
#'
#' @description
#' Plots the enrichment for individual genes within your GeneWalk results.
#' This will plot the top_n enriched pathways for the gene of interest.
#'
#' @param object The `GeneWalk` class. The stats need to be available in the
#' object.
#' @param gene_of_interest String. The gene of interest for which to plot
#' the enriched pathways.
#' @param top_n Integer. Number of genes to plot
#'
#' @returns ggplot2 object with the plot for the gene
#'
#' @export
plot_individual_gene <- S7::new_generic(
  name = "plot_individual_gene",
  dispatch_args = "object",
  fun = function(
    object,
    gene_of_interest,
    top_n = 10L
  ) {
    S7::S7_dispatch()
  }
)

#' @method plot_individual_gene GeneWalk
#'
#' @import ggplot2
S7::method(plot_individual_gene, GeneWalk) <- function(
  object,
  gene_of_interest,
  top_n = 10L
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  checkmate::qassert(gene_of_interest, "S1")
  checkmate::qassert(top_n, "I1")

  # early return
  if (suppressWarnings(nrow(get_stats(genewalk_obj)) == 0)) {
    warning(paste(
      "No stats found in the object.",
      "Did you run calculate_genewalk_stats()?",
      "Returning NULL."
    ))
    return(NULL)
  }

  plot_df <- data.table::copy(get_stats(genewalk_obj))

  # throw an error if a gene that is not in the data was used
  checkmate::assertTRUE(gene_of_interest %in% plot_df$gene)

  plot_df <- get_stats(obj)[gene == gene_of_interest]
  setorder(plot_df, avg_pval)

  plot_df <- plot_df[1:top_n, ][,
    pathway := factor(pathway, levels = rev(pathway))
  ]

  p <- ggplot(
    data = plot_df,
    mapping = aes(
      x = -log10(avg_gene_fdr + .Machine$double.eps),
      xend = 0,
      y = pathway,
      yend = pathway
    )
  ) +
    geom_col(
      mapping = aes(fill = pathway),
      colour = "black",
      show.legend = FALSE
    ) +
    geom_segment(aes(
      x = -log10(gene_fdr_ci_lower + .Machine$double.eps),
      xend = -log10(gene_fdr_ci_upper + .Machine$double.eps)
    )) +
    xlab("-log10(FDR)") +
    ylab("Pathway") +
    ggtitle(sprintf("Genewalk result - %s", gene_of_interest)) +
    scale_fill_viridis_d(option = "C") +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.5,
      color = "darkred"
    ) +
    theme_bw()

  return(p)
}
