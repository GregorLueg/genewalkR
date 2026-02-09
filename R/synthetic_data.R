# synthetic data ---------------------------------------------------------------

## node2vec --------------------------------------------------------------------

#' Generate Simple Node2Vec Test Data
#'
#' @description Creates graphs with known structure to verify node2vec is
#' working correctly. The following options exists:
#' \itemize{
#'  \item `"barbell"` - Two dense clusters connected by single edge.
#'  \item `"caveman"` - Multiple isolated cliques with sparse inter-connections.
#'  \item `"stochastic_block"` - Communities with controlled edge probabilities.
#' }
#'
#' @param test_type String. One of
#' `c("barbell", "caveman", "stochastic_block")`
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
  test_type = c("barbell", "caveman", "stochastic_block"),
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
    c("barbell", "caveman", "stochastic_block")
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
#'   \item{within_cluster_sim}{Mean similarity within clusters}
#'   \item{between_cluster_sim}{Mean similarity between clusters}
#'   \item{within_dist}{The within cluster similarities as vector}
#'   \item{between_dist}{The between cluster similarities as vector}
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
#' @param ppi_params List. PPI parameters from [params_ppi()].
#' @param pathway_params List. Pathway parameters from [params_pathway()].
#' @param n_communities Integer. Number of gene communities. Defaults to `3L`.
#' @param signal_strength Numeric. Probability of connecting to focal pathways
#'   (0-1). Defaults to `0.8`.
#' @param seed Integer. Random seed. Defaults to `42L`.
#'
#' @returns A list with edges, ground truth, and pathway metadata.
#'
#' @export
synthetic_genewalk_data <- function(
  ppi_params = params_ppi(),
  pathway_params = params_pathway(),
  n_communities = 3L,
  signal_strength = 0.8,
  seed = 42L
) {
  assertPPIParams(ppi_params)
  assertPathwayParams(pathway_params)
  checkmate::qassert(n_communities, "I1")
  checkmate::qassert(signal_strength, "N1[0, 1]")
  checkmate::qassert(seed, "I1")

  if (n_communities * ppi_params$min_community_size > ppi_params$n_genes) {
    stop(
      "Cannot create ",
      n_communities,
      " communities with minimum ",
      ppi_params$min_community_size,
      " genes each from ",
      ppi_params$n_genes,
      " total genes"
    )
  }

  set.seed(seed)

  gene_ids <- sprintf("gene_%04d", seq_len(ppi_params$n_genes))
  ppi_graph <- igraph::sample_pa(
    ppi_params$n_genes,
    power = 1,
    m = ppi_params$ppi_m,
    directed = FALSE
  )
  igraph::V(ppi_graph)$name <- gene_ids

  communities <- igraph::cluster_louvain(ppi_graph)
  gene_communities <- igraph::membership(communities)

  unique_comms <- unique(gene_communities)
  if (length(unique_comms) > n_communities) {
    comm_sizes <- table(gene_communities)
    keep_comms <- names(sort(comm_sizes, decreasing = TRUE)[seq_len(
      n_communities
    )])
    reassign_to <- as.integer(keep_comms[1])
    gene_communities[
      !gene_communities %in% as.integer(keep_comms)
    ] <- reassign_to
  }

  comm_sizes <- table(gene_communities)
  while (any(comm_sizes < ppi_params$min_community_size)) {
    small_comms <- as.integer(names(comm_sizes[
      comm_sizes < ppi_params$min_community_size
    ]))
    large_comm <- as.integer(names(which.max(comm_sizes)))

    for (small_comm in small_comms) {
      gene_communities[gene_communities == small_comm] <- large_comm
    }

    comm_sizes <- table(gene_communities)
    unique_comms <- unique(gene_communities)

    if (length(unique_comms) < n_communities) {
      stop(
        "Cannot maintain ",
        n_communities,
        " communities with minimum size ",
        ppi_params$min_community_size,
        ". Reduce n_communities or min_community_size."
      )
    }
  }

  ppi_edges <- igraph::as_edgelist(ppi_graph)
  ppi_dt <- data.table(
    from = ppi_edges[, 1],
    to = ppi_edges[, 2],
    edge_type = "interacts"
  )

  # call rust helper for pathway generation
  pathway_result <- rs_generate_pathway_structure(
    n_pathways = pathway_params$n_pathways,
    pathway_depth = pathway_params$pathway_depth,
    pathway_branching = pathway_params$pathway_branching,
    n_communities = n_communities,
    gene_ids = gene_ids,
    gene_communities = as.integer(gene_communities),
    n_focal_pathways = pathway_params$n_focal_pathways,
    signal_strength = signal_strength,
    connections_per_gene = pathway_params$connections_per_gene,
    seed = seed
  )

  pathway_dt <- data.table(
    from = pathway_result$pathway_edges_from,
    to = pathway_result$pathway_edges_to,
    edge_type = "parent_of"
  )

  pathway_metadata <- data.table(
    pathway = pathway_result$pathway_ids,
    subtree = pathway_result$pathway_subtrees,
    depth = pathway_result$pathway_depths,
    is_hub = pathway_result$is_hub
  )

  gene_pathway_dt <- data.table(
    from = pathway_result$gene_pathway_edges_from,
    to = pathway_result$gene_pathway_edges_to,
    edge_type = "part_of"
  )

  # reconstruct community_focal_pathways as named list
  cfp_dt <- data.table(
    community = pathway_result$cfp_communities,
    pathway = pathway_result$cfp_pathways
  )
  community_focal_pathways <- split(cfp_dt$pathway, cfp_dt$community)

  # reconstruct comm_to_subtree as named vector
  comm_to_subtree <- setNames(
    pathway_result$cts_subtrees,
    as.character(pathway_result$cts_communities)
  )

  all_edges <- rbindlist(list(ppi_dt, pathway_dt, gene_pathway_dt))

  ground_truth <- data.table(
    gene = gene_ids,
    community = gene_communities
  )

  list(
    edges = all_edges,
    ground_truth = ground_truth,
    pathway_metadata = pathway_metadata,
    community_focal_pathways = community_focal_pathways,
    comm_to_subtree = comm_to_subtree
  )
}


#' Get Expected Gene-Pathway Associations
#'
#' Returns gene-pathway pairs that should show high cosine similarity based on
#' the underlying community-subtree structure in synthetic data.
#'
#' @param synthetic_data List. Output from
#' [genewalkR::synthetic_genewalk_data()].
#' @param return_negatives Logical. If TRUE, also returns gene-pathway pairs
#' from different communities/subtrees as negative controls.
#' @param exclude_hubs Logical. If TRUE, excludes hub pathways (depth <= 1)
#' from evaluation as they have high similarity to everything.
#'
#' @return data.table with columns: gene, pathway, expected_signal
#'
#' @export
#'
#' @import data.table
get_expected_associations <- function(
  synthetic_data,
  return_negatives = FALSE,
  exclude_hubs = TRUE
) {
  checkmate::assertList(synthetic_data)
  checkmate::assertNames(
    names(synthetic_data),
    must.include = c(
      "ground_truth",
      "pathway_metadata",
      "community_focal_pathways",
      "comm_to_subtree"
    )
  )

  ground_truth <- copy(synthetic_data$ground_truth)
  pathway_metadata <- copy(synthetic_data$pathway_metadata)
  community_focal_pathways <- synthetic_data$community_focal_pathways
  comm_to_subtree <- synthetic_data$comm_to_subtree

  # filter out hub pathways if requested
  if (exclude_hubs) {
    pathway_metadata <- pathway_metadata[is_hub == FALSE]
  }

  positive_pairs_list <- list()

  for (comm_id in names(community_focal_pathways)) {
    focal_pathways <- community_focal_pathways[[comm_id]]

    if (exclude_hubs) {
      focal_pathways <- intersect(
        focal_pathways,
        pathway_metadata$pathway
      )
    }

    if (length(focal_pathways) == 0) {
      next
    }

    genes_in_comm <- ground_truth[community == as.integer(comm_id), gene]

    positive_pairs_list[[comm_id]] <- CJ(
      gene = genes_in_comm,
      pathway = focal_pathways
    )
  }

  positive_pairs <- rbindlist(positive_pairs_list)
  positive_pairs[, expected_signal := TRUE]

  if (!return_negatives) {
    return(positive_pairs)
  }

  negative_pairs_list <- list()

  for (comm_id in names(comm_to_subtree)) {
    genes_in_comm <- ground_truth[community == as.integer(comm_id), gene]
    target_subtree <- comm_to_subtree[comm_id]

    other_pathways <- pathway_metadata[
      subtree != target_subtree,
      pathway
    ]

    if (length(other_pathways) == 0) {
      next
    }

    n_negatives_per_gene <- ceiling(
      nrow(positive_pairs) / length(genes_in_comm)
    )

    for (gene in genes_in_comm) {
      n_sample <- min(n_negatives_per_gene, length(other_pathways))
      sampled_pathways <- sample(other_pathways, n_sample)

      negative_pairs_list[[paste0(comm_id, "_", gene)]] <- data.table(
        gene = gene,
        pathway = sampled_pathways
      )
    }
  }

  negative_pairs <- rbindlist(negative_pairs_list)
  negative_pairs[, expected_signal := FALSE]

  if (nrow(negative_pairs) > nrow(positive_pairs)) {
    negative_pairs <- negative_pairs[sample(.N, nrow(positive_pairs))]
  }

  rbindlist(list(positive_pairs, negative_pairs))
}
