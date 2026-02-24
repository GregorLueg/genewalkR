# genewalkR main functions -----------------------------------------------------

## initial embedding generation ------------------------------------------------

#' Generate the initial graph embedding
#'
#' @description
#' This function will generate the initial node2vec embedding with the provided
#' parameters. In brief, it will generate random walks over the graph and fit
#' a Skip Gram model for the number of defined epochs to generate on a per node
#' basis an embedding vector. These will be subsequently used for the
#' calculation of test statistics. For more details, please refer to Ietswaart,
#' et al.
#'
#' @param object The `GeneWalk` class, please see [genewalkR::GeneWalk()].
#' @param embd_dim Integer. Size of the embedding dimensions to create. Defaults
#' to `8L` in line with the authors recommendations.
#' @param node2vec_params Named list. Contains the node2vec parameters. The
#' list has the following elements:
#' \itemize{
#'  \item p - Numeric. Return parameter for biased random walks. Defaults to
#'  `1.0`.
#'  \item q - Numeric. In-out parameter for biased random walks. Defaults to
#'  `1.0.`.
#'  \item walks_per_node - Integer. Number of random walks per node. Defaults to
#'  ` 40L.`.
#'  \item walk_length - Integer. Length of each random walk. Defaults to `40L`.
#'  \item num_workers - Number of worker threads during batching. Defaults to
#'  `4L`. To note: during graph generation, Rust will use via the Rayon backend
#'  all available threads.
#'  \item batch_size - Integer. Batch size for training. Defaults to `256L`.
#'  \item n_epochs - Integer. Number of training epochs. Defaults to `15L`.
#'  \item n_negatives - Integer. Number of negative samples. Defaults to `5L`.
#'  \item window_size - Integer.  Context window size. Defaults to `2L`.
#'  \item lr - Numeric. Learning rate. Defaults to `1e-2`.
#' }
#' @param directed Boolean. Indicates if this is a directed or undirected
#' network. Defaults to `FALSE`.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function
#'
#' @returns The object with the added initial embeddings.
#'
#' @export
#'
#' @references Ietswaart, et al., Genome Biol, 2021
generate_initial_emb <- S7::new_generic(
  name = "generate_initial_emb",
  dispatch_args = "object",
  fun = function(
    object,
    embd_dim = 8L,
    node2vec_params = params_node2vec(),
    directed = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_initial_emb GeneWalk
#'
#' @export
S7::method(generate_initial_emb, GeneWalk) <- function(
  object,
  embd_dim = 8L,
  node2vec_params = params_node2vec(),
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  assertNode2VecParam(node2vec_params)
  checkmate::qassert(directed, "B1")
  checkmate::qassert(.verbose, "B1")

  graph_dt <- S7::prop(object, "graph_dt")

  # function body
  nodes <- sort(unique(c(graph_dt$from, graph_dt$to)))
  from_idx <- match(graph_dt$from, nodes)
  to_idx <- match(graph_dt$to, nodes)

  weights <- if ("weights" %in% names(graph_dt)) {
    graph_dt$weights
  } else {
    NULL
  }

  # generate the embedding
  embd <- rs_gene_walk(
    from = from_idx,
    to = to_idx,
    weights = weights,
    gene_walk_params = node2vec_params,
    embd_dim = embd_dim,
    directed = directed,
    seed = seed,
    verbose = .verbose
  )

  rownames(embd) <- nodes
  colnames(embd) <- sprintf("emb_dim_%i", 1:embd_dim)

  full_params <- node2vec_params
  full_params[["embd_size"]] <- embd_dim
  full_params[["directed"]] <- directed

  S7::prop(object, "embd") <- embd
  S7::prop(object, "params")[["node2vec"]] <- full_params

  return(object)
}

## generate random embeddings --------------------------------------------------

#' Generate the random, permuted embedding
#'
#' @description
#' This function will generate generate permuted subgraphs that show the same
#' node degree and run the node2vec algorithm on top of these subgraphs. The
#' cosine distance between the real edges will be calculated and stored in
#' the object. These will be used subsequently for statistical significance
#' testing. For more details, please refer to Ietswaart, et al. The function
#' will return a warning and the class as is if you have not run
#' [genewalkR::generate_initial_emb()].
#'
#' @param object The `GeneWalk` class, please see [genewalkR::GeneWalk()].
#' @param n_perm Inter. Number of permutations to generate. Defaults to `3L`,
#' the recommendation of the authors.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function
#'
#' @returns The object with the added random embeddings for statistical testing.
#'
#' @details
#' The function will be using the same parameters for the node2vec training
#' as were used for [genewalkR::generate_initial_emb()].
#'
#' @export
#'
#' @references Ietswaart, et al., Genome Biol, 2021
generate_permuted_emb <- S7::new_generic(
  name = "generate_permuted_emb",
  dispatch_args = "object",
  fun = function(
    object,
    n_perm = 3L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_permuted_emb GeneWalk
#'
#' @export
S7::method(generate_permuted_emb, GeneWalk) <- function(
  object,
  n_perm = 3L,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  checkmate::qassert(n_perm, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # early return
  if (is.null(S7::prop(object, "params")[["node2vec"]])) {
    warning(paste(
      "No initial embedding was generated. Returning object as is.",
      "Please run generate_initial_emb()."
    ))

    return(object)
  }

  graph_dt <- S7::prop(object, "graph_dt")

  n_gene_pathway_edges <- nrow(graph_dt[type == "part_of"])

  nodes <- sort(unique(c(graph_dt$from, graph_dt$to)))
  from_idx <- match(graph_dt$from, nodes)
  to_idx <- match(graph_dt$to, nodes)

  weights <- if ("weights" %in% names(graph_dt)) {
    graph_dt$weights
  } else {
    NULL
  }

  # generate random embeddings — now returns list of matrices
  rnd_embd <- rs_gene_walk_perm(
    from = from_idx,
    to = to_idx,
    weights = weights,
    gene_walk_params = S7::prop(object, "params")[["node2vec"]],
    n_gene_pathway_edges = n_gene_pathway_edges,
    n_perm = n_perm,
    embd_dim = S7::prop(object, "params")[["node2vec"]][["embd_size"]],
    directed = S7::prop(object, "params")[["node2vec"]][["directed"]],
    seed = seed,
    verbose = .verbose
  )

  S7::prop(object, "permuted_embd") <- rnd_embd

  return(object)
}

## calculate statistics --------------------------------------------------------

#' Calculate the pathway and gene statistics
#'
#' @description
#' Calculates statistical significance for gene-pathway associations using
#' the GeneWalk algorithm. For each gene-pathway pair, this function computes
#' p-values and FDR corrections (both globally across all pairs and per-gene)
#' across multiple network permutations. The final output includes mean values
#' and 95% confidence intervals for all statistics.
#'
#' The function processes all permutations in Rust for optimal performance,
#' calculating cosine similarities between gene and pathway embeddings,
#' determining p-values based on null distributions from randomised networks,
#' and applying Benjamini-Hochberg FDR corrections.
#'
#' @param object The `GeneWalk` object containing actual and permuted
#' embeddings. See [genewalkR::GeneWalk()].
#' @param .verbose Logical. If `TRUE`, prints progress messages during
#' calculation. Default is `TRUE`.
#'
#' @returns The input object with an added `stats` property containing a
#' data.table with the following columns for each gene-pathway pair:
#' \itemize{
#'   \item gene - Gene identifier
#'   \item pathway - Pathway identifier
#'   \item avg_pval - Mean p-value across permutations
#'   \item pval_ci_lower, pval_ci_upper - 95% CI for p-values
#'   \item avg_global_fdr - Mean global FDR across permutations
#'   \item global_fdr_ci_lower,global_fdr_ci_upper - 95% CI for global FDR
#'   \item avg_gene_fdr - Mean gene-wise FDR across permutations
#'   \item gene_fdr_ci_lower,gene_fdr_ci_upper - 95% CI for gene-wise FDR
#' }
#'
#' @export
#'
#' @references Ietswaart, et al., Genome Biol, 2021
calculate_genewalk_stats <- S7::new_generic(
  name = "calculate_genewalk_stats",
  dispatch_args = "object",
  fun = function(
    object,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_genewalk_stats GeneWalk
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @import data.table
S7::method(calculate_genewalk_stats, GeneWalk) <- function(
  object,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, GeneWalk))
  checkmate::qassert(.verbose, "B1")

  # early return
  if (length(S7::prop(object, "permuted_embd")) == 0) {
    warning(paste(
      "No permuted embeddings found. Returning object as is.",
      "Please run generate_permuted_emb()."
    ))

    return(object)
  }

  actual_embd <- S7::prop(object, "embd")
  permuted_embds <- S7::prop(object, "permuted_embd")
  gene_ids <- S7::prop(object, "gene_ids")
  pathway_ids <- S7::prop(object, "pathway_ids")
  gene_to_pathways <- S7::prop(object, "gene_to_pathway_dt")

  connected <- split(
    match(gene_to_pathways$to, pathway_ids),
    factor(match(gene_to_pathways$from, gene_ids), levels = seq_along(gene_ids))
  )

  # Row indices into the full embedding (1-based for R/Rust bridge)
  all_nodes <- rownames(actual_embd)
  gene_indices <- match(gene_ids, all_nodes)
  pathway_indices <- match(pathway_ids, all_nodes)

  stats <- rs_gene_walk_test(
    gene_embds = actual_embd[gene_ids, ],
    pathway_embds = actual_embd[pathway_ids, ],
    null_similarities = permuted_embds,
    connected_pathways = connected,
    verbose = .verbose
  )

  stats <- data.table::as.data.table(stats) %>%
    .[, `:=`(
      gene = gene_ids[gene],
      pathway = pathway_ids[pathway]
    )]

  S7::prop(object, "stats") <- stats

  return(object)
}
