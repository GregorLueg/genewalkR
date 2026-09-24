# metapath2vec tests -----------------------------------------------------------

## test data -------------------------------------------------------------------

# three gene communities, each with dense gene-gene edges and its own two
# pathways; one gene-gene edge bridges consecutive communities
n_comm <- 3L
n_genes <- 10L

genes <- data.table::data.table(
  id = sprintf("g%02i", seq_len(n_comm * n_genes)),
  type = "gene",
  cluster = rep(seq_len(n_comm), each = n_genes)
)
pathways <- data.table::data.table(
  id = sprintf("p%02i", seq_len(n_comm * 2L)),
  type = "pathway",
  cluster = rep(seq_len(n_comm), each = 2L)
)
node_dt <- rbind(genes, pathways)

gene_gene <- genes[,
  {
    pairs <- utils::combn(id, 2L)
    list(from = pairs[1, ], to = pairs[2, ])
  },
  by = cluster
]
gene_pathway <- merge(
  genes[, .(from = id, cluster)],
  pathways[, .(to = id, cluster)],
  by = "cluster",
  allow.cartesian = TRUE
)
bridges <- data.table::data.table(
  from = genes[cluster < n_comm, id[1], by = cluster]$V1,
  to = genes[cluster > 1L, id[n_genes], by = cluster]$V1
)
graph_dt <- rbind(
  gene_gene[, .(from, to)],
  gene_pathway[, .(from, to)],
  bridges
)

params <- params_metapath2vec(
  walks_per_node = 20L,
  walk_length = 20L,
  n_epochs = 10L
)

## tests -----------------------------------------------------------------------

### output shape ---------------------------------------------------------------

res <- metapath2vec(
  graph_dt = graph_dt,
  node_dt = node_dt[, .(id, type)],
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params,
  .verbose = FALSE
)

expect_equal(
  current = dim(res),
  target = c(nrow(node_dt), 8L),
  info = "metapath2vec - embedding dimensions"
)

expect_true(
  current = setequal(rownames(res), node_dt$id),
  info = "metapath2vec - rownames are the node set"
)

expect_equal(
  current = as.character(attr(res, "node_type")[node_dt$id]),
  target = node_dt$type,
  info = "metapath2vec - node_type attribute matches the node table"
)

walk_stats <- attr(res, "walk_stats")
expect_equal(
  current = walk_stats$start_nodes,
  target = nrow(genes),
  info = "metapath2vec - walks start on every gene"
)
expect_equal(
  current = walk_stats$dropped,
  target = 0L,
  info = "metapath2vec - no dropped walks on a fully connected schema"
)

expect_equal(
  current = attr(res, "unvisited_nodes"),
  target = character(0),
  info = "metapath2vec - every node visited on a fully connected schema"
)

### unvisited nodes ------------------------------------------------------------

# a pathway hanging off another pathway has no genes, so gene-pathway-gene
# never reaches it
node_dt_orphan <- rbind(
  node_dt[, .(id, type)],
  data.table::data.table(id = "p_orphan", type = "pathway")
)
graph_dt_orphan <- rbind(
  graph_dt,
  data.table::data.table(from = "p_orphan", to = "p01")
)

res_filtered <- metapath2vec(
  graph_dt = graph_dt_orphan,
  node_dt = node_dt_orphan,
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params,
  .verbose = FALSE
)

expect_equal(
  current = attr(res_filtered, "unvisited_nodes"),
  target = "p_orphan",
  info = "metapath2vec - orphan pathway is reported as unvisited"
)

expect_true(
  current = setequal(
    c(rownames(res_filtered), attr(res_filtered, "unvisited_nodes")),
    node_dt_orphan$id
  ),
  info = "metapath2vec - rows plus unvisited nodes give the node set"
)

expect_equal(
  current = names(attr(res_filtered, "node_type")),
  target = rownames(res_filtered),
  info = "metapath2vec - node_type attribute follows the filtered rows"
)

res_unfiltered <- metapath2vec(
  graph_dt = graph_dt_orphan,
  node_dt = node_dt_orphan,
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params,
  filter_unvisited = FALSE,
  .verbose = FALSE
)

expect_true(
  current = "p_orphan" %in% rownames(res_unfiltered) &&
    identical(attr(res_unfiltered, "unvisited_nodes"), "p_orphan"),
  info = "metapath2vec - filter_unvisited = FALSE keeps the row, still reports it"
)

### cluster recovery -----------------------------------------------------------

gene_labels <- genes[, .(node = id, cluster)]

metrics <- evaluate_node2vec_test(
  embeddings = res[gene_labels$node, ],
  node_labels = gene_labels
)

expect_true(
  current = metrics$between_cluster_sim < metrics$within_cluster_sim,
  info = "metapath2vec - genes cluster by community"
)

res_plus <- metapath2vec(
  graph_dt = graph_dt,
  node_dt = node_dt[, .(id, type)],
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params,
  metapath_plus = TRUE,
  .verbose = FALSE
)

metrics_plus <- evaluate_node2vec_test(
  embeddings = res_plus[gene_labels$node, ],
  node_labels = gene_labels
)

expect_true(
  current = metrics_plus$between_cluster_sim <
    metrics_plus$within_cluster_sim,
  info = "metapath2vec++ - genes cluster by community"
)

### seed reproducibility -------------------------------------------------------

# single thread, the hogwild gradient updates race otherwise
params_single <- params_metapath2vec(
  walks_per_node = 20L,
  walk_length = 20L,
  n_epochs = 10L,
  num_workers = 1L
)

res_1 <- metapath2vec(
  graph_dt = graph_dt,
  node_dt = node_dt[, .(id, type)],
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params_single,
  .verbose = FALSE
)

res_2 <- metapath2vec(
  graph_dt = graph_dt,
  node_dt = node_dt[, .(id, type)],
  metapath = c("gene", "pathway", "gene"),
  metapath2vec_params = params_single,
  .verbose = FALSE
)

expect_equal(
  current = res_1,
  target = res_2,
  info = "metapath2vec - reproducibility with num_workers set to 1"
)

### input validation -----------------------------------------------------------

expect_error(
  current = metapath2vec(
    graph_dt = graph_dt,
    node_dt = node_dt[, .(id, type)],
    metapath = c("gene", "term", "gene"),
    .verbose = FALSE
  ),
  info = "metapath2vec - unknown type in metapath"
)

expect_error(
  current = metapath2vec(
    graph_dt = graph_dt,
    node_dt = node_dt[, .(id, type)],
    metapath = c("gene", "pathway", "pathway"),
    .verbose = FALSE
  ),
  pattern = "close",
  info = "metapath2vec - open metapath"
)

expect_error(
  current = metapath2vec(
    graph_dt = rbind(graph_dt, data.table::data.table(from = "g01", to = "x")),
    node_dt = node_dt[, .(id, type)],
    metapath = c("gene", "pathway", "gene"),
    .verbose = FALSE
  ),
  pattern = "endpoints",
  info = "metapath2vec - edge to unknown node"
)

expect_error(
  current = metapath2vec(
    graph_dt = graph_dt,
    node_dt = node_dt[, .(id, type)],
    metapath = c("gene", "pathway", "gene"),
    metapath2vec_params = list(walks_per_node = 10),
    .verbose = FALSE
  ),
  info = "metapath2vec - malformed params"
)

expect_error(
  current = params_metapath2vec(walk_length = 1L),
  info = "params_metapath2vec - walk_length below 2"
)
