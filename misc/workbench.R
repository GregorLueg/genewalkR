test_data <- differential_graph_test_data(
  n_stable = 50L,
  n_comm2 = 50L,
  n_comm3 = 50L,
  n_exclusive = 5L
)

# data_info holds the ground truth
head(test_data$data_info)

edges = test_data$g1_edges
nodes = test_data$g1_nodes
info = test_data$data_info

make_graph <- function(edges, nodes, info) {
  g <- igraph::graph_from_data_frame(d = edges, directed = FALSE)

  node_names <- igraph::V(g)$name
  status <- info[match(node_names, info$node), node_status]
  differential <- info[match(node_names, info$node), is_differential]

  layout <- igraph::layout_with_kk(g)

  igraph::V(g)$x <- layout[, 1]
  igraph::V(g)$y <- layout[, 2]

  igraph::V(g)$color <- data.table::fcase(
    grepl("only", status)             , "grey80"  ,
    differential & status == "shared" , "#E74C3C" ,
    default = "#3498DB"
  )
  g
}

g1 <- make_graph(
  edges = test_data$g1_edges,
  nodes = test_data$g1_nodes,
  info = test_data$data_info
)

g2 <- make_graph(
  edges = test_data$g2_edges,
  nodes = test_data$g2_nodes,
  info = test_data$data_info
)

par(mfrow = c(1, 2))
plot(g1, vertex.label = NA, main = "Graph 1", vertex.size = 6)
plot(g2, vertex.label = NA, main = "Graph 2", vertex.size = 6)
par(mfrow = c(1, 1))

obj <- EmbedDrift(
  graph_dt_1 = test_data$g1_edges,
  graph_dt_2 = test_data$g2_edges
)

print(obj)

obj <- generate_initial_embeddings(
  object = obj,
  embd_dim = 8L,
  node2vec_params = params_node2vec(window_size = 2L),
  .verbose = TRUE
)

print(obj)

obj <- calculate_drift(object = obj)

stats <- get_stats(obj)
stats


results <- merge(
  stats,
  test_data$data_info[, .(node, is_differential)],
  by = "node"
)

results[node_status == "shared"] |>
  ggplot(aes(
    x = is_differential,
    y = cosine_similarity,
    fill = is_differential
  )) +
  geom_boxplot(width = 0.4, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(width = 0.1, shape = 21, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("FALSE" = "#3498DB", "TRUE" = "#E74C3C")) +
  scale_x_discrete(labels = c("FALSE" = "Stable", "TRUE" = "Differential")) +
  labs(
    fill = "Pos Ctr"
  ) +
  theme_minimal() +
  xlab("Pos Ctr") +
  ylab("Cosine sim") +
  ggtitle("Context drift by ground truth label")
