# generate a nice plotting function --------------------------------------------

test <- myc_gwn

test

table(test@graph_dt$type)

object = myc_gwn

results <- get_stats(object)

graph_dt <- get_graph_dt(object)


node_dt <- .get_node_degree(graph_dt)

plot_dt <- results[,
  .(
    ratio = sum(avg_gene_fdr <= 0.1) / length(pathway),
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

ggplotly(p, tooltip = c("text", "x", "y", "size")) |>
  layout(legend = list(y = 0.5))
