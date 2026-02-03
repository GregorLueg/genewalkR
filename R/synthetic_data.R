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

#' Generate Synthetic GeneWalk Data
#'
#' Creates a synthetic dataset combining a scale-free PPI network with a
#' hierarchical pathway DAG. Genes within communities share specific focal
#' pathways to enable signal recovery testing.
#'
#' @param n_genes Integer. Number of genes to generate.
#' @param n_pathways Integer. Number of pathways to generate.
#' @param n_communities Integer. Number of gene communities in PPI network.
#' @param signal_strength Numeric. Probability (0-1) that genes connect to
#' their community's focal pathways vs random pathways.
#' @param seed Integer. Random seed for reproducibility.
#' @param ppi_m Integer. Number of edges attached per new node in
#' Barabási-Albert model.
#' @param pathway_depth Integer. Depth of pathway hierarchy.
#' @param pathway_branching Integer. Average number of children per parent
#' pathway.
#' @param n_focal_pathways Integer. Number of focal pathways assigned to each
#' community.
#' @param connections_per_gene Integer. Number of pathway connections per gene.
#' @param min_community_size Integer. Minimum genes per community.
#'
#' @return List containing:
#'   \item{edges}{data.table with columns: from, to, edge_type}
#'   \item{ground_truth}{data.table mapping genes to communities}
#'   \item{pathway_metadata}{data.table with pathway properties (depth, subtree,
#'     is_hub)}
#'   \item{community_focal_pathways}{List mapping community IDs to focal
#'     pathway sets}
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
  seed = 42L,
  ppi_m = 3L,
  pathway_depth = 4L,
  pathway_branching = 3L,
  n_focal_pathways = 15L,
  connections_per_gene = 3L,
  min_community_size = 20L
) {
  checkmate::qassert(n_genes, "I1")
  checkmate::qassert(n_pathways, "I1")
  checkmate::qassert(n_communities, "I1")
  checkmate::qassert(signal_strength, "N1[0, 1]")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(ppi_m, "I1")
  checkmate::qassert(pathway_depth, "I1")
  checkmate::qassert(pathway_branching, "I1")
  checkmate::qassert(n_focal_pathways, "I1")
  checkmate::qassert(connections_per_gene, "I1")
  checkmate::qassert(min_community_size, "I1")

  if (n_communities * min_community_size > n_genes) {
    stop(
      "Cannot create ",
      n_communities,
      " communities with minimum ",
      min_community_size,
      " genes each from ",
      n_genes,
      " total genes"
    )
  }

  set.seed(seed)

  # generate PPI network with controlled communities
  gene_ids <- sprintf("gene_%04d", seq_len(n_genes))
  ppi_graph <- igraph::sample_pa(
    n_genes,
    power = 1,
    m = ppi_m,
    directed = FALSE
  )
  igraph::V(ppi_graph)$name <- gene_ids

  # detect communities and consolidate to target number
  communities <- igraph::cluster_louvain(ppi_graph)
  gene_communities <- igraph::membership(communities)

  unique_comms <- unique(gene_communities)
  if (length(unique_comms) > n_communities) {
    comm_sizes <- table(gene_communities)
    keep_comms <- names(sort(comm_sizes, decreasing = TRUE)[seq_len(
      n_communities
    )])

    # reassign small communities to largest community
    reassign_to <- as.integer(keep_comms[1])
    gene_communities[!gene_communities %in% as.integer(keep_comms)] <-
      reassign_to
  }

  # ensure minimum community sizes by reassigning if necessary
  comm_sizes <- table(gene_communities)
  while (any(comm_sizes < min_community_size)) {
    small_comms <- as.integer(names(comm_sizes[
      comm_sizes < min_community_size
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
        min_community_size,
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

  # generate pathway DAG with depth tracking
  pathway_edges <- list()
  pathway_depth_map <- integer()
  pathway_subtree <- integer()
  pathway_ancestors <- list()
  pathway_ids_built <- character()

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
    pathway_depth_map <- c(pathway_depth_map, 0L)
    pathway_ancestors[[root_name]] <- character(0)

    current_level <- root_name
    current_id <- current_id + 1L

    for (depth in seq_len(pathway_depth - 1L)) {
      if (current_id > n_pathways) {
        break
      }
      next_level <- character()

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
          pathway_depth_map <- c(pathway_depth_map, depth)

          pathway_edges[[length(pathway_edges) + 1L]] <- c(parent, child_name)
          pathway_ancestors[[child_name]] <- c(
            pathway_ancestors[[parent]],
            parent
          )

          next_level <- c(next_level, child_name)
          current_id <- current_id + 1L

          # occasional second parent for DAG structure
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

  pathway_ids <- pathway_ids_built
  n_pathways_actual <- length(pathway_ids)

  pathway_dt <- data.table(
    from = vapply(pathway_edges, `[`, character(1), 1L),
    to = vapply(pathway_edges, `[`, character(1), 2L),
    edge_type = "parent_of"
  )

  # build pathway metadata
  pathway_metadata <- data.table(
    pathway = pathway_ids,
    subtree = pathway_subtree,
    depth = pathway_depth_map
  )

  # mark hub pathways (roots and first level)
  pathway_metadata[, is_hub := depth <= 1L]

  # assign focal pathways to each community
  actual_communities <- sort(unique(gene_communities))
  actual_subtrees <- unique(pathway_subtree)
  n_mappable <- min(length(actual_communities), length(actual_subtrees))

  comm_to_subtree <- setNames(
    actual_subtrees[seq_len(n_mappable)],
    as.character(actual_communities[seq_len(n_mappable)])
  )

  community_focal_pathways <- list()

  for (comm_id in names(comm_to_subtree)) {
    target_subtree <- comm_to_subtree[comm_id]
    target_pathways <- pathway_metadata[subtree == target_subtree]

    # prefer mid-to-deep level pathways (not roots, not necessarily leaves)
    candidate_pathways <- target_pathways[depth >= 2L, pathway]

    if (length(candidate_pathways) == 0) {
      candidate_pathways <- target_pathways[depth >= 1L, pathway]
    }

    if (length(candidate_pathways) == 0) {
      candidate_pathways <- target_pathways$pathway
    }

    n_to_pick <- min(n_focal_pathways, length(candidate_pathways))
    community_focal_pathways[[comm_id]] <- sample(candidate_pathways, n_to_pick)
  }

  # gene-pathway associations with GO propagation
  gene_pathway_edges <- list()

  for (gene_idx in seq_len(n_genes)) {
    gene <- gene_ids[gene_idx]
    gene_comm <- as.character(gene_communities[gene_idx])
    focal_pathways <- community_focal_pathways[[gene_comm]]

    # fallback if community has no focal pathways
    if (is.null(focal_pathways) || length(focal_pathways) == 0) {
      focal_pathways <- sample(
        pathway_ids,
        min(n_focal_pathways, n_pathways_actual)
      )
    }

    connected_pathways <- character()

    for (i in seq_len(connections_per_gene)) {
      if (runif(1) < signal_strength) {
        pathway <- sample(focal_pathways, 1L)
      } else {
        pathway <- sample(pathway_ids, 1L)
      }
      connected_pathways <- c(connected_pathways, pathway)
    }

    connected_pathways <- unique(connected_pathways)

    # add connections with GO propagation to all ancestors
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

  list(
    edges = all_edges,
    ground_truth = ground_truth,
    pathway_metadata = pathway_metadata,
    community_focal_pathways = community_focal_pathways,
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
