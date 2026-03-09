# synthetic data ---------------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Generate Simple Node2Vec Test Data
#'
#' @description Creates graphs with known structure to verify node2vec is
#' working correctly. The following options exists:
#' \itemize{
#'  \item `"barbell"` - Two dense clusters connected by single edge.
#'  \item `"cavemen"` - Multiple isolated cliques with sparse inter-connections.
#'  \item `"stochastic_block"` - Communities with controlled edge probabilities.
#' }
#'
#' @param test_type String. One of
#' `c("barbell", "cavemen", "stochastic_block")`
#' @param n_nodes_per_cluster Integer. Nodes per cluster/community.
#' @param n_clusters Integer. Number of clusters.
#' @param p_within Numeric. Probability of edges within cluster (0-1).
#' @param p_between Numeric. Probability of edges between clusters (0-1).
#' @param seed Integer. Random seed.
#'
#' @return A list containing:
#' \item{edges}{data.table with from, to, edge_type}
#' \item{node_labels}{data.table with node, cluster}
#'
#' @export
#'
#' @import data.table
node2vec_test_data <- function(
  test_type = c("barbell", "cavemen", "stochastic_block"),
  n_nodes_per_cluster = 10L,
  n_clusters = 3L,
  p_within = 0.8,
  p_between = 0.05,
  seed = 42L
) {
  test_type <- match.arg(test_type)

  # checks
  checkmate::assertChoice(
    test_type,
    c("barbell", "cavemen", "stochastic_block")
  )
  checkmate::qassert(n_nodes_per_cluster, "I1")
  checkmate::qassert(n_clusters, "I1")
  checkmate::qassert(p_within, "N1[0, 1]")
  checkmate::qassert(p_between, "N1[0, 1]")
  checkmate::qassert(seed, "I1")

  zeallot::`%<-%`(
    c(edges_list, nodes_list),
    rs_node2vec_synthetic_data(
      test_data = test_type,
      n_nodes_per_cluster = n_nodes_per_cluster,
      n_clusters = n_clusters,
      p_within = p_within,
      p_between = p_between,
      seed = seed
    )
  )

  edges_dt <- data.table::setDT(edges_list)
  node_labels <- data.table::setDT(nodes_list)

  list(
    edges = edges_dt,
    node_labels = node_labels
  )
}

#' Evaluate node2vec test results
#'
#' Compares embedding similarities against ground truth clusters.
#'
#' @param embeddings Matrix of node embeddings (nodes x dimensions).
#' @param node_labels data.table with columns: node, cluster.
#'
#' @return A list with:
#' \itemize{
#'   \item within_cluster_sim - Mean similarity within clusters.
#'   \item between_cluster_sim - Mean similarity between clusters.
#'   \item within_dist - The within cluster similarities as vector.
#'   \item between_dist - The between cluster similarities as vector.
#' }
#'
#' @export
#'
#' @import data.table
evaluate_node2vec_test <- function(embeddings, node_labels) {
  # checks
  checkmate::assertMatrix(embeddings, mode = "numeric")
  checkmate::assertDataTable(node_labels)
  checkmate::assertNames(
    names(node_labels),
    must.include = c("node", "cluster")
  )

  n_nodes <- nrow(embeddings)

  # compute cosine similarities for all pairs
  within_sims <- c()
  between_sims <- c()

  for (i in 1:(n_nodes - 1)) {
    for (j in (i + 1):n_nodes) {
      sim <- rs_cosine_sim(embeddings[i, ], embeddings[j, ])
      cluster_i <- node_labels$cluster[i]
      cluster_j <- node_labels$cluster[j]

      if (cluster_i == cluster_j) {
        within_sims <- c(within_sims, sim)
      } else {
        between_sims <- c(between_sims, sim)
      }
    }
  }

  mean_within <- mean(within_sims, na.rm = TRUE)
  mean_between <- mean(between_sims, na.rm = TRUE)

  list(
    within_cluster_sim = mean_within,
    between_cluster_sim = mean_between,
    within_dist = within_sims,
    between_dist = between_sims
  )
}

## gene walk -------------------------------------------------------------------

#' Generate synthetic GeneWalk data
#'
#' @description
#' Generates a synthetic network that mimics real GeneWalk input data for
#' testing and benchmarking. The network consists of three components: a
#' hierarchical ontology (tree structure with lateral edges), gene-to-ontology
#' annotations, and protein-protein interactions (PPI).
#'
#' Signal genes are annotated exclusively within a single ontology subtree and
#' form PPI edges among themselves, creating coherent graph neighbourhoods.
#' Noise genes receive the same number of annotations but are forced to span
#' multiple subtrees, and their PPI edges are random across the noise group.
#' Both groups are degree-matched so that significance is driven by annotation
#' coherence rather than connectivity alone.
#'
#' @param data_params Named list. Return from
#' [genewalkR::params_genewalk_data()].
#' @param seed Integer. Random seed. Defaults to `42L`.
#'
#' @returns A list with the following items:
#' \describe{
#'   \item{full_data}{A data.table with columns `from` and `to` containing
#'   all edges (ontology hierarchy, gene-to-pathway, and PPI).}
#'   \item{gene_to_pathways}{A data.table with columns `from` (gene) and
#'   `to` (pathway) for gene-to-ontology annotation edges only.}
#'   \item{gene_ids}{Character vector of all gene identifiers (signal and
#'   noise).}
#'   \item{pathway_ids}{Character vector of all ontology term identifiers.}
#' }
#'
#' @export
synthetic_genewalk_data <- function(
  data_params = params_genewalk_data(),
  seed = 42L
) {
  assertGeneWalkDataParams(data_params)
  checkmate::qassert(seed, "I1")

  data <- with(
    data_params,
    rs_build_synthetic_genewalk(
      n_signal_genes = n_signal_genes,
      n_noise_genes = n_noise_genes,
      n_roots = n_roots,
      depth = depth,
      branching = branching,
      p_lateral = p_lateral,
      p_ppi = p_ppi,
      min_annotations = min_annotations,
      max_annotations = max_annotations,
      min_noise_subtrees = min_noise_subtrees,
      seed = seed
    )
  )

  ontology <- data.table::as.data.table(data[c("ont_from", "ont_to")]) %>%
    `colnames<-`(c("from", "to")) %>%
    .[, type := "hierarchy"]
  pathway_genes <- data.table::as.data.table(data[c(
    "gene_ont_from",
    "gene_ont_to"
  )]) %>%
    `colnames<-`(c("from", "to")) %>%
    .[, type := "part_of"]
  ppi <- data.table::as.data.table(data[c(
    "ppi_from",
    "ppi_to"
  )]) %>%
    `colnames<-`(c("from", "to")) %>%
    .[, type := "interaction"]

  full_dt <- data.table::rbindlist(list(ontology, pathway_genes, ppi))

  pathways <- (ontology %$% unique(c(from, to)))
  genes <- unique(c(pathway_genes$from, (ppi %$% c(from, to))))

  res <- list(
    full_data = full_dt,
    gene_to_pathways = pathway_genes,
    gene_ids = genes,
    pathway_ids = pathways
  )
}

## diff cor graphs -------------------------------------------------------------

#' Generate Synthetic Differential Graph Data
#'
#' @description
#' Creates a pair of graphs with known topological differences to verify
#' context-aware differential embedding analysis. The two graphs share a common
#' backbone but differ in defined regions:
#'
#' \itemize{
#'   \item Community 1 is a stable negative control, identical across both
#'   graphs.
#'   \item Community 2 contains a hub node demoted to a peripheral node in
#'   graph 2.
#'   \item Two bridge nodes span communities 2 and 3 in graph 1 but are fully
#'   embedded within community 3 in graph 2.
#'   \item A set of exclusive nodes appears in only one graph.
#' }
#'
#' @param n_stable Integer. Number of nodes in the stable negative-control
#' community.
#' @param n_comm2 Integer. Number of regular nodes in community 2, excluding
#' the hub.
#' @param n_comm3 Integer. Number of nodes in community 3, excluding bridge
#' nodes.
#' @param n_exclusive Integer. Number of exclusive nodes per graph.
#'
#' @return A named list containing:
#' \itemize{
#'   \item g1_edges - data.table with columns `from`, `to` for graph 1.
#'   \item g1_nodes - data.table with columns `node`, `cluster` for graph 1.
#'   \item g2_edges - data.table with columns `from`, `to` for graph 2.
#'   \item g2_nodes - data.table with columns `node`, `cluster` for graph 2.
#'   \item data_info - data.table with columns `node`, `is_differential`,
#'   `node_status` (`"shared"`, `"g1_only"`, `"g2_only"`).
#' }
#'
#' @export
differential_graph_test_data <- function(
  n_stable = 50L,
  n_comm2 = 50L,
  n_comm3 = 50L,
  n_exclusive = 5L
) {
  # checks
  checkmate::qassert(n_stable, "I1[2,)")
  checkmate::qassert(n_comm2, "I1[2,)")
  checkmate::qassert(n_comm3, "I1[2,)")
  checkmate::qassert(n_exclusive, "I1[0,)")

  raw <- rs_differential_graph_data(
    n_stable = n_stable,
    n_comm2 = n_comm2,
    n_comm3 = n_comm3,
    n_exclusive = n_exclusive
  )

  info <- data.table::data.table(
    node = raw$data_info$shared_nodes,
    is_differential = raw$data_info$is_differential,
    node_status = "shared"
  )
  info <- data.table::rbindlist(list(
    info,
    data.table::data.table(
      node = raw$data_info$g1_only,
      is_differential = TRUE,
      node_status = "g1_only"
    ),
    data.table::data.table(
      node = raw$data_info$g2_only,
      is_differential = TRUE,
      node_status = "g2_only"
    )
  ))

  list(
    g1_edges = data.table::setDT(raw$g1_edges),
    g1_nodes = data.table::setDT(raw$g1_nodes),
    g2_edges = data.table::setDT(raw$g2_edges),
    g2_nodes = data.table::setDT(raw$g2_nodes),
    data_info = info
  )
}
