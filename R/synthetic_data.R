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

  if (test_type == "barbell") {
    n_total <- n_nodes_per_cluster * 2
    node_ids <- sprintf("node_%03d", 1:n_total)

    edges_list <- list()

    for (i in 1:(n_nodes_per_cluster - 1)) {
      for (j in (i + 1):n_nodes_per_cluster) {
        edges_list[[length(edges_list) + 1]] <- c(node_ids[i], node_ids[j])
      }
    }

    offset <- n_nodes_per_cluster
    for (i in 1:(n_nodes_per_cluster - 1)) {
      for (j in (i + 1):n_nodes_per_cluster) {
        edges_list[[length(edges_list) + 1]] <- c(
          node_ids[offset + i],
          node_ids[offset + j]
        )
      }
    }

    edges_list[[length(edges_list) + 1]] <- c(
      node_ids[n_nodes_per_cluster],
      node_ids[n_nodes_per_cluster + 1]
    )

    edges_dt <- data.table(
      from = sapply(edges_list, `[`, 1),
      to = sapply(edges_list, `[`, 2),
      edge_type = "connects"
    )

    node_labels <- data.table(
      node = node_ids,
      cluster = c(rep(1, n_nodes_per_cluster), rep(2, n_nodes_per_cluster))
    )
  } else if (test_type == "caveman") {
    n_total <- n_nodes_per_cluster * n_clusters
    node_ids <- sprintf("node_%03d", 1:n_total)

    edges_list <- list()

    for (cluster in 1:n_clusters) {
      offset <- (cluster - 1) * n_nodes_per_cluster
      for (i in 1:(n_nodes_per_cluster - 1)) {
        for (j in (i + 1):n_nodes_per_cluster) {
          edges_list[[length(edges_list) + 1]] <- c(
            node_ids[offset + i],
            node_ids[offset + j]
          )
        }
      }
    }

    for (cluster1 in 1:(n_clusters - 1)) {
      for (cluster2 in (cluster1 + 1):n_clusters) {
        if (runif(1) < 0.3) {
          node1 <- sample(
            (cluster1 - 1) * n_nodes_per_cluster + (1:n_nodes_per_cluster),
            1
          )
          node2 <- sample(
            (cluster2 - 1) * n_nodes_per_cluster + (1:n_nodes_per_cluster),
            1
          )
          edges_list[[length(edges_list) + 1]] <- c(
            node_ids[node1],
            node_ids[node2]
          )
        }
      }
    }

    edges_dt <- data.table(
      from = sapply(edges_list, `[`, 1),
      to = sapply(edges_list, `[`, 2),
      edge_type = "connects"
    )

    node_labels <- data.table(
      node = node_ids,
      cluster = rep(1:n_clusters, each = n_nodes_per_cluster)
    )
  } else if (test_type == "stochastic_block") {
    n_total <- n_nodes_per_cluster * n_clusters
    node_ids <- sprintf("node_%03d", 1:n_total)

    edges_list <- list()

    for (i in 1:(n_total - 1)) {
      for (j in (i + 1):n_total) {
        cluster_i <- ceiling(i / n_nodes_per_cluster)
        cluster_j <- ceiling(j / n_nodes_per_cluster)

        prob <- if (cluster_i == cluster_j) p_within else p_between

        if (runif(1) < prob) {
          edges_list[[length(edges_list) + 1]] <- c(node_ids[i], node_ids[j])
        }
      }
    }

    edges_dt <- data.table(
      from = sapply(edges_list, `[`, 1),
      to = sapply(edges_list, `[`, 2),
      edge_type = "connects"
    )

    node_labels <- data.table(
      node = node_ids,
      cluster = rep(1:n_clusters, each = n_nodes_per_cluster)
    )
  }

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
#' @param verbose Boolean. Print diagnostic information.
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

#' Generate Synthetic GeneWalk Data
#'
#' Creates a synthetic dataset combining a scale-free PPI network with a
#' hierarchical pathway DAG, where gene communities are preferentially
#' associated with specific pathway subtrees to enable signal recovery testing.
#'
#' @param n_genes Integer. Number of genes to generate.
#' @param n_pathways Integer. Number of pathways to generate.
#' @param n_communities Integer. Target number of gene communities in PPI
#' network.
#' @param signal_strength Numeric. Probability (0-1) that genes connect to
#' their associated pathway subtree.
#' @param noise_level Numeric. Baseline probability (0-1) for random
#' gene-pathway associations.
#' @param seed Integer. Random seed for reproducibility.
#' @param ppi_m Integer. Number of edges attached per new node in
#' Barabási-Albert model.
#' @param pathway_depth Integer. Depth of pathway hierarchy.
#' @param pathway_branching Integer. Average number of children per parent
#' pathway.
#'
#' @return List containing:
#'   \item{edges}{data.table with columns: from, to, edge_type}
#'   \item{ground_truth}{data.table mapping genes to communities}
#'   \item{pathway_groups}{data.table mapping pathways to subtrees}
#'   \item{comm_to_subtree}{Named vector mapping community IDs to subtree IDs}
#'
#' @export
#'
#' @import data.table
synthetic_genewalk_data <- function(
  n_genes = 500L,
  n_pathways = 250L,
  n_communities = 3L,
  signal_strength = 0.8,
  noise_level = 0.05,
  seed = 42L,
  ppi_m = 3L,
  pathway_depth = 4L,
  pathway_branching = 3L
) {
  checkmate::qassert(n_genes, "I1")
  checkmate::qassert(n_pathways, "I1")
  checkmate::qassert(n_communities, "I1")
  checkmate::qassert(signal_strength, "N1[0, 1]")
  checkmate::qassert(noise_level, "N1[0, 1]")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(ppi_m, "I1")
  checkmate::qassert(pathway_depth, "I1")
  checkmate::qassert(pathway_branching, "I1")

  set.seed(seed)

  # Generate PPI network
  gene_ids <- sprintf("gene_%04d", seq_len(n_genes))
  ppi_graph <- igraph::sample_pa(
    n_genes,
    power = 1,
    m = ppi_m,
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
    gene_communities[!gene_communities %in% keep_comms] <- keep_comms[1]
  }

  ppi_edges <- igraph::as_edgelist(ppi_graph)
  ppi_dt <- data.table(
    from = ppi_edges[, 1],
    to = ppi_edges[, 2],
    edge_type = "interacts"
  )

  # Generate pathway DAG - build dynamically to avoid orphans
  pathway_edges <- list()
  pathway_subtree <- integer(0)
  pathway_ancestors <- list()
  pathway_ids_built <- character(0)

  n_roots <- max(2L, n_communities)
  current_id <- 1L
  subtree_id <- 1L

  for (root_idx in seq_len(n_roots)) {
    if (current_id > n_pathways) {
      break
    }

    root_name <- sprintf("pathway_%04d", current_id)
    pathway_ids_built <- c(pathway_ids_built, root_name)
    pathway_subtree <- c(pathway_subtree, subtree_id)
    pathway_ancestors[[root_name]] <- character(0)

    current_level <- root_name
    current_id <- current_id + 1L

    for (depth in seq_len(pathway_depth - 1L)) {
      if (current_id > n_pathways) {
        break
      }
      next_level <- character(0)

      for (parent in current_level) {
        n_children <- sample(
          max(1L, pathway_branching - 1L):(pathway_branching + 1L),
          1L
        )

        for (child_idx in seq_len(n_children)) {
          if (current_id > n_pathways) {
            break
          }

          child_name <- sprintf("pathway_%04d", current_id)
          pathway_ids_built <- c(pathway_ids_built, child_name)
          pathway_subtree <- c(pathway_subtree, subtree_id)

          pathway_edges[[length(pathway_edges) + 1L]] <- c(parent, child_name)
          pathway_ancestors[[child_name]] <- c(
            pathway_ancestors[[parent]],
            parent
          )

          next_level <- c(next_level, child_name)
          current_id <- current_id + 1L

          # Add second parent occasionally (DAG structure)
          eligible_alt_parents <- setdiff(current_level, parent)
          if (runif(1) < 0.2 && length(eligible_alt_parents) > 0) {
            alt_parent <- sample(eligible_alt_parents, 1L)
            pathway_edges[[length(pathway_edges) + 1L]] <- c(
              alt_parent,
              child_name
            )
            pathway_ancestors[[child_name]] <- unique(c(
              pathway_ancestors[[child_name]],
              pathway_ancestors[[alt_parent]],
              alt_parent
            ))
          }
        }
      }
      current_level <- next_level
      if (length(current_level) == 0) break
    }
    subtree_id <- subtree_id + 1L
  }

  # Use only the pathways we actually built
  pathway_ids <- pathway_ids_built
  n_pathways_actual <- length(pathway_ids)

  pathway_dt <- data.table(
    from = vapply(pathway_edges, `[`, character(1), 1L),
    to = vapply(pathway_edges, `[`, character(1), 2L),
    edge_type = "parent_of"
  )

  # Gene-pathway associations with GO propagation
  gene_pathway_edges <- list()

  actual_subtrees <- unique(pathway_subtree)
  actual_communities <- unique(gene_communities)
  n_mappable <- min(length(actual_communities), length(actual_subtrees))

  comm_to_subtree <- setNames(
    actual_subtrees[seq_len(n_mappable)],
    as.character(actual_communities[seq_len(n_mappable)])
  )

  for (gene_idx in seq_len(n_genes)) {
    gene <- gene_ids[gene_idx]
    gene_comm <- as.character(gene_communities[gene_idx])
    target_subtree <- comm_to_subtree[gene_comm]

    n_connections <- rpois(1L, lambda = 3) + 1L
    connected_pathways <- character(0)

    for (i in seq_len(n_connections)) {
      use_signal <- runif(1) < signal_strength &&
        !is.na(target_subtree) &&
        length(target_subtree) > 0

      if (use_signal) {
        target_pathways <- pathway_ids[pathway_subtree == target_subtree]
        if (length(target_pathways) > 0) {
          connected_pathways <- c(
            connected_pathways,
            sample(target_pathways, 1L)
          )
        }
      } else if (runif(1) < noise_level) {
        connected_pathways <- c(connected_pathways, sample(pathway_ids, 1L))
      }
    }

    # Ensure every gene connects to at least one pathway
    if (length(connected_pathways) == 0) {
      if (!is.na(target_subtree) && length(target_subtree) > 0) {
        target_pathways <- pathway_ids[pathway_subtree == target_subtree]
        if (length(target_pathways) > 0) {
          connected_pathways <- sample(target_pathways, 1L)
        }
      }
      if (length(connected_pathways) == 0) {
        connected_pathways <- sample(pathway_ids, 1L)
      }
    }

    # For each connected pathway, also connect to ALL ancestors (GO propagation)
    for (pathway in connected_pathways) {
      gene_pathway_edges[[length(gene_pathway_edges) + 1L]] <- c(gene, pathway)

      ancestors <- pathway_ancestors[[pathway]]
      for (ancestor in ancestors) {
        gene_pathway_edges[[length(gene_pathway_edges) + 1L]] <- c(
          gene,
          ancestor
        )
      }
    }
  }

  gene_pathway_edges <- unique(gene_pathway_edges)

  gene_pathway_dt <- data.table(
    from = vapply(gene_pathway_edges, `[`, character(1), 1L),
    to = vapply(gene_pathway_edges, `[`, character(1), 2L),
    edge_type = "part_of"
  )

  all_edges <- rbindlist(list(ppi_dt, pathway_dt, gene_pathway_dt))

  ground_truth <- data.table(
    gene = gene_ids,
    community = gene_communities
  )

  pathway_groups <- data.table(
    pathway = pathway_ids,
    subtree = pathway_subtree
  )

  list(
    edges = all_edges,
    ground_truth = ground_truth,
    pathway_groups = pathway_groups,
    comm_to_subtree = comm_to_subtree,
    pathway_hierarchy = pathway_ancestors
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
#'
#' @return data.table with columns: gene, pathway, expected_signal (TRUE for
#' same community/subtree, FALSE for different)
#'
#' @export
#'
#' @import data.table
get_expected_associations <- function(
  synthetic_data,
  return_negatives = FALSE
) {
  checkmate::assertList(synthetic_data)
  checkmate::assertNames(
    names(synthetic_data),
    must.include = c("ground_truth", "pathway_groups", "comm_to_subtree")
  )

  ground_truth <- copy(synthetic_data$ground_truth)
  pathway_groups <- copy(synthetic_data$pathway_groups)
  comm_to_subtree <- synthetic_data$comm_to_subtree

  # Map genes to their target subtree
  ground_truth[, target_subtree := comm_to_subtree[as.character(community)]]

  # Positive pairs: inner join to only keep pathways with mapped genes
  positive_pairs <- pathway_groups[
    ground_truth,
    on = .(subtree = target_subtree),
    allow.cartesian = TRUE,
    nomatch = NULL
  ][, .(gene, pathway, expected_signal = TRUE)]

  if (!return_negatives) {
    return(positive_pairs)
  }

  # Negative pairs: cross join then filter out matches
  all_pairs <- CJ(
    gene = ground_truth$gene,
    pathway = pathway_groups$pathway
  )

  # Add subtree info
  all_pairs[ground_truth, target_subtree := i.target_subtree, on = .(gene)]
  all_pairs[pathway_groups, pathway_subtree := i.subtree, on = .(pathway)]

  # Keep only mismatches where both have valid subtrees
  negative_pairs <- all_pairs[
    !is.na(target_subtree) &
      !is.na(pathway_subtree) &
      target_subtree != pathway_subtree
  ]

  # Sample if too many
  n_negatives <- nrow(positive_pairs)
  if (nrow(negative_pairs) > n_negatives) {
    negative_pairs <- negative_pairs[sample(.N, n_negatives)]
  }

  negative_pairs <- negative_pairs[, .(gene, pathway, expected_signal = FALSE)]

  rbindlist(list(positive_pairs, negative_pairs))
}
