pub mod data;
pub mod embedding;
pub mod graph;
pub mod utils;

use extendr_api::prelude::*;
use faer::Mat;
use std::time::Instant;

use crate::data::*;
use crate::embedding::*;
use crate::graph::*;
use crate::utils::*;

extendr_module! {
    mod genewalkR;
    fn rs_gene_walk;
    fn rs_gene_walk_perm;
    fn rs_gene_walk_test;
    fn rs_cosine_sim;
    fn rs_node2vec_synthetic_data;
    fn rs_build_synthetic_genewalk;
}

/////////////////////////
// Gene Walk functions //
/////////////////////////

/// Generate GeneWalk node embeddings
///
/// @description Uses a SIMD-accelerated CPU implementation of word2vec
/// with negative sampling to learn node representations from biased
/// random walks (node2vec).
///
/// @param from Integer vector. Node indices for edge origins.
/// @param to Integer vector. Node indices for edge destinations.
/// @param weights Optional numeric vector. Edge weights, defaults to 1.0.
/// @param gene_walk_params Named list. Training parameters (p, q,
///   walks_per_node, walk_length, num_workers, n_epochs, num_negatives,
///   window_size, lr, dim).
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Treat graph as directed.
/// @param seed Integer. Random seed.
/// @param verbose Boolean. Controls verbosity.
///
/// @return A numeric matrix of n_nodes x embedding dimensions.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> RMatrix<f64> {
    let mut config = GeneWalkConfig::from_r_list(gene_walk_params, seed);
    config.train_args.dim = embd_dim;

    if verbose {
        println!("Preparing the network edges.");
    }
    let start_edge_prep = Instant::now();

    let from = from.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let edges = prepare_edges(from, to, weights);

    if verbose {
        println!(
            "Processed the edges in {:.2}s",
            start_edge_prep.elapsed().as_secs_f64()
        );
        println!("Running the random walks.");
    }

    let start_rw = Instant::now();

    let node2vec_graph = create_graph(&edges, config.p, config.q, directed);

    let vocab_size = edges
        .iter()
        .flat_map(|(from, to, _)| [*from, *to])
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0);

    let random_walks =
        node2vec_graph.generate_walks(config.walks_per_node, config.walk_length, config.seed);

    if verbose {
        println!(
            "Generated the random walks in {:.2}s",
            start_rw.elapsed().as_secs_f64()
        );
        println!(
            "Training the word2vec model for {} epochs",
            config.train_args.epochs
        );
    }

    let start_training = Instant::now();

    let embedding = train_node2vec(random_walks, vocab_size, &config, verbose);

    if verbose {
        println!(
            "Training completed in {:.2}s.",
            start_training.elapsed().as_secs_f64()
        );
    }

    let nrows = embedding.len();
    let ncols = embedding[0].len();

    RMatrix::new_matrix(nrows, ncols, |r, c| embedding[r][c] as f64)
}

//////////////////
// Permutations //
//////////////////

/// Generate permuted embeddings for null distribution
///
/// @description Generates degree-preserving random networks and trains node2vec
/// on each, returning the raw embedding matrices for downstream statistical
/// testing.
///
/// @param from Integer vector. Node indices for edge origins.
/// @param to Integer vector. Node indices for edge destinations.
/// @param weights Optional numeric vector. Edge weights, defaults to 1.0.
/// @param gene_walk_params Named list. Training parameters.
/// @param n_perm Integer. Number of permutations.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Treat graph as directed.
/// @param seed Integer. Random seed.
/// @param verbose Boolean. Controls verbosity.
///
/// @returns A list of n_perm embedding matrices (each n_nodes x embd_dim).
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk_perm(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    n_perm: usize,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let mut config = GeneWalkConfig::from_r_list(gene_walk_params, seed);
    config.train_args.dim = embd_dim;

    if verbose {
        println!("Preparing the network edges.");
    }
    let start_edge_prep = Instant::now();

    let from = from.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let edges = prepare_edges(from, to, weights);

    if verbose {
        println!(
            "Processed the edges in {:.2}s",
            start_edge_prep.elapsed().as_secs_f64()
        );
    }

    let vocab_size = edges
        .iter()
        .flat_map(|(from, to, _)| [*from, *to])
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0);

    let mut res_ls = List::new(n_perm);

    for perm in 0..n_perm {
        let start_perm = Instant::now();

        if verbose {
            println!("Starting permutation: {}", perm + 1);
            println!("  Generating randomised edges and random walks.");
        }

        let start_rw = Instant::now();

        let perm_edges = generate_random_network(&edges, directed, (seed + perm) as u64);
        let node2vec_graph = create_graph(&perm_edges, config.p, config.q, directed);
        let random_walks_perm =
            node2vec_graph.generate_walks(config.walks_per_node, config.walk_length, config.seed);

        if verbose {
            println!(
                "  Random walks on permuted data in {:.2}s",
                start_rw.elapsed().as_secs_f64()
            );
            println!(
                "  Training word2vec model on permuted network for {} epochs",
                config.train_args.epochs
            );
        }

        let embd_i = train_node2vec(random_walks_perm, vocab_size, &config, verbose);

        // return as R matrix — same shape as the real embedding
        let nrows = embd_i.len();
        let ncols = embd_i[0].len();
        let mat = RMatrix::new_matrix(nrows, ncols, |r, c| embd_i[r][c] as f64);

        res_ls.set_elt(perm, Robj::from(mat))?;

        if verbose {
            println!(
                "  Permutation {} completed in {:.2}s",
                perm + 1,
                start_perm.elapsed().as_secs_f64()
            );
        }
    }

    Ok(res_ls)
}

/////////////////////
// Test statistics //
/////////////////////

/// Calculate the test statistics
///
/// @description Calculates test statistics for gene-pathway pairs. The null
/// distribution is derived from gene-pathway cosine similarities in permuted
/// embeddings (full cross-cosine, not just connected pairs).
///
/// @param gene_embds Matrix of n_genes x embedding dimensions.
/// @param pathway_embds Matrix of n_pathways x embedding dimensions.
/// @param permuted_embds List of permuted embedding matrices (n_nodes x dim).
/// @param gene_indices Integer vector. 1-based row indices for genes in
/// permuted embeddings.
/// @param pathway_indices Integer vector. 1-based row indices for pathways in
/// permuted embeddings.
/// @param connected_pathways List. Gene to pathway connections (1-indexed).
/// @param verbose Controls verbosity.
///
/// @returns A list with per-pair statistics (see original docs).
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk_test(
    gene_embds: RMatrix<f64>,
    pathway_embds: RMatrix<f64>,
    permuted_embds: List,
    gene_indices: Vec<i32>,
    pathway_indices: Vec<i32>,
    connected_pathways: List,
    verbose: bool,
) -> List {
    let gene_embds = r_matrix_to_vec(gene_embds);
    let pathway_embds = r_matrix_to_vec(pathway_embds);

    let n_genes = gene_embds.len();
    let n_perms = permuted_embds.len();

    // Convert 1-based R indices to 0-based
    let gene_idx: Vec<usize> = gene_indices.iter().map(|&x| (x - 1) as usize).collect();
    let pathway_idx: Vec<usize> = pathway_indices.iter().map(|&x| (x - 1) as usize).collect();

    let connectivity: Vec<Vec<usize>> = connected_pathways
        .iter()
        .map(|(_, v)| {
            v.as_integer_vector()
                .unwrap_or_default()
                .iter()
                .map(|&x| (x - 1) as usize)
                .collect()
        })
        .collect();

    // observed cosine similarities (gene embeddings vs pathway embeddings)
    let cosine_sim: Mat<f64> = compute_cross_cosine(&gene_embds, &pathway_embds);

    if verbose {
        let total_pairs: usize = connectivity.iter().map(|c| c.len()).sum();
        println!("Testing {} connected gene-pathway pairs", total_pairs);
    }

    // build null distribution and compute p-values per permutation
    let mut raw_pvals_per_perm: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_perms);
    let mut gene_fdr_per_perm: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_perms);
    let mut global_fdr_per_perm: Vec<Vec<f64>> = Vec::with_capacity(n_perms);

    for (perm_i, perm_item) in permuted_embds.iter().enumerate() {
        if verbose {
            println!("Processing permutation {}/{}", perm_i + 1, n_perms);
        }

        let perm_mat = r_matrix_to_vec(RMatrix::try_from(&perm_item.1).unwrap());

        // extract gene and pathway rows from permuted embedding
        let perm_genes: Vec<Vec<f32>> = gene_idx.iter().map(|&i| perm_mat[i].clone()).collect();
        let perm_pathways: Vec<Vec<f32>> =
            pathway_idx.iter().map(|&i| perm_mat[i].clone()).collect();

        // full gene x pathway cross-cosine on the permuted embedding = null
        let null_cross = compute_cross_cosine(&perm_genes, &perm_pathways);
        let mut null_vec: Vec<f64> = Vec::with_capacity(null_cross.nrows() * null_cross.ncols());
        for r in 0..null_cross.nrows() {
            for c in 0..null_cross.ncols() {
                null_vec.push(null_cross[(r, c)]);
            }
        }

        // compute p-vals observed similarities against this permutation's null
        let pvals_mat: Mat<f64> = calculate_p_vals(&cosine_sim.as_ref(), null_vec);

        // extract connected p-values per gene
        let gene_pvals: Vec<Vec<f64>> = (0..n_genes)
            .map(|gi| {
                connectivity[gi]
                    .iter()
                    .map(|&pi| pvals_mat[(gi, pi)])
                    .collect()
            })
            .collect();

        // global FDR over all connected pairs
        let all_connected_pvals: Vec<f64> = gene_pvals.iter().flatten().cloned().collect();
        let global_fdr_flat = calc_fdr(&all_connected_pvals);

        // per-gene FDR
        let gene_fdr: Vec<Vec<f64>> = gene_pvals.iter().map(|p| calc_fdr(p)).collect();

        raw_pvals_per_perm.push(gene_pvals);
        gene_fdr_per_perm.push(gene_fdr);
        global_fdr_per_perm.push(global_fdr_flat);
    }

    // aggregate across permutations
    let mut out_gene: Vec<i32> = Vec::new();
    let mut out_pathway: Vec<i32> = Vec::new();
    let mut out_similarity: Vec<f64> = Vec::new();
    let mut out_avg_pval: Vec<f64> = Vec::new();
    let mut out_pval_ci_lower: Vec<f64> = Vec::new();
    let mut out_pval_ci_upper: Vec<f64> = Vec::new();
    let mut out_avg_global_fdr: Vec<f64> = Vec::new();
    let mut out_global_fdr_ci_lower: Vec<f64> = Vec::new();
    let mut out_global_fdr_ci_upper: Vec<f64> = Vec::new();
    let mut out_avg_gene_fdr: Vec<f64> = Vec::new();
    let mut out_gene_fdr_ci_lower: Vec<f64> = Vec::new();
    let mut out_gene_fdr_ci_upper: Vec<f64> = Vec::new();

    let mut global_flat_idx = 0usize;

    for gene_i in 0..n_genes {
        for (local_idx, &pathway_i) in connectivity[gene_i].iter().enumerate() {
            let raw_pvals: Vec<f64> = (0..n_perms)
                .map(|p| raw_pvals_per_perm[p][gene_i][local_idx])
                .collect();
            let gene_fdrs: Vec<f64> = (0..n_perms)
                .map(|p| gene_fdr_per_perm[p][gene_i][local_idx])
                .collect();
            let global_fdrs: Vec<f64> = (0..n_perms)
                .map(|p| global_fdr_per_perm[p][global_flat_idx])
                .collect();

            let (avg_pval, pval_lo, pval_hi) = log_stats(&raw_pvals);
            let (avg_gfdr, gfdr_lo, gfdr_hi) = log_stats(&global_fdrs);
            let (avg_gene_fdr, gene_lo, gene_hi) = log_stats(&gene_fdrs);

            out_gene.push((gene_i + 1) as i32);
            out_pathway.push((pathway_i + 1) as i32);
            out_similarity.push(cosine_sim[(gene_i, pathway_i)]);
            out_avg_pval.push(avg_pval);
            out_pval_ci_lower.push(pval_lo);
            out_pval_ci_upper.push(pval_hi);
            out_avg_global_fdr.push(avg_gfdr);
            out_global_fdr_ci_lower.push(gfdr_lo);
            out_global_fdr_ci_upper.push(gfdr_hi);
            out_avg_gene_fdr.push(avg_gene_fdr);
            out_gene_fdr_ci_lower.push(gene_lo);
            out_gene_fdr_ci_upper.push(gene_hi);

            global_flat_idx += 1;
        }
    }

    list![
        gene = out_gene,
        pathway = out_pathway,
        similarity = out_similarity,
        avg_pval = out_avg_pval,
        pval_ci_lower = out_pval_ci_lower,
        pval_ci_upper = out_pval_ci_upper,
        avg_global_fdr = out_avg_global_fdr,
        global_fdr_ci_lower = out_global_fdr_ci_lower,
        global_fdr_ci_upper = out_global_fdr_ci_upper,
        avg_gene_fdr = out_avg_gene_fdr,
        gene_fdr_ci_lower = out_gene_fdr_ci_lower,
        gene_fdr_ci_upper = out_gene_fdr_ci_upper
    ]
}

///////////
// Utils //
///////////

/// Cosine similarity between two vectors
///
/// @description Computes cosine similarity between two numeric vectors.
///
/// @param a Numeric vector.
/// @param b Numeric vector.
///
/// @returns Cosine similarity between the two vectors.
///
/// @export
#[extendr]
fn rs_cosine_sim(a: &[f64], b: &[f64]) -> f64 {
    cosine_similarity(a, b)
}

/// Generate synthetic data for node2vec
///
/// @param test_data String. One of
///   c("barbell", "caveman", "stochastic_block").
/// @param n_nodes_per_cluster Integer. Nodes per cluster.
/// @param n_clusters Integer. Number of clusters.
/// @param p_within Numeric. Within-cluster edge probability (0-1).
/// @param p_between Numeric. Between-cluster edge probability (0-1).
/// @param seed Integer. Random seed.
///
/// @returns A list with edges and nodes.
///
/// @export
#[extendr]
fn rs_node2vec_synthetic_data(
    test_data: String,
    n_nodes_per_cluster: usize,
    n_clusters: usize,
    p_within: f64,
    p_between: f64,
    seed: usize,
) -> List {
    let data_type = parse_node2vec_data(&test_data).unwrap_or_default();

    let synthetic_data = match data_type {
        Node2VecDataType::Barbell => node2vec_barbell(n_nodes_per_cluster),
        Node2VecDataType::Cavemen => {
            node2vec_caveman(n_nodes_per_cluster, n_clusters, p_between, seed)
        }
        Node2VecDataType::Stochastic => {
            node2vec_stochastic_block(n_nodes_per_cluster, n_clusters, p_within, p_between, seed)
        }
    };

    let (edge_data, node_data) = synthetic_data.generate_lists();

    list!(edges = edge_data, nodes = node_data)
}

/// Generate synthetic GeneWalk data with controlled signal structure
///
/// @param n_signal_genes Integer. Genes annotated to a single ontology subtree.
/// @param n_noise_genes Integer. Genes with annotations scattered across
///   subtrees.
/// @param n_roots Integer. Number of ontology root terms.
/// @param depth Integer. Depth of each ontology subtree.
/// @param branching Integer. Average branching factor per node.
/// @param p_lateral Numeric. Probability of lateral edges within each ontology
///   level.
/// @param p_ppi Numeric. PPI connection probability within gene groups.
/// @param min_annotations Integer. Minimum annotations per gene.
/// @param max_annotations Integer. Maximum annotations per gene.
/// @param min_noise_subtrees Integer. Minimum number of different subtrees
///   each noise gene must span.
/// @param seed Integer. Random seed.
///
/// @returns A list with: ontology_edges, gene_ont_edges, ppi_edges,
///   signal_genes, noise_genes.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_build_synthetic_genewalk(
    n_signal_genes: usize,
    n_noise_genes: usize,
    n_roots: usize,
    depth: usize,
    branching: usize,
    p_lateral: f64,
    p_ppi: f64,
    min_annotations: usize,
    max_annotations: usize,
    min_noise_subtrees: usize,
    seed: usize,
) -> List {
    let data = build_synthetic_genewalk(
        n_signal_genes,
        n_noise_genes,
        n_roots,
        depth,
        branching,
        p_lateral,
        p_ppi,
        min_annotations,
        max_annotations,
        min_noise_subtrees,
        seed as u64,
    );

    list!(
        ont_from = data
            .ontology_edges
            .iter()
            .map(|(f, _)| f.clone())
            .collect::<Vec<_>>(),
        ont_to = data
            .ontology_edges
            .iter()
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>(),
        gene_ont_from = data
            .gene_ont_edges
            .iter()
            .map(|(f, _)| f.clone())
            .collect::<Vec<_>>(),
        gene_ont_to = data
            .gene_ont_edges
            .iter()
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>(),
        ppi_from = data
            .ppi_edges
            .iter()
            .map(|(f, _)| f.clone())
            .collect::<Vec<_>>(),
        ppi_to = data
            .ppi_edges
            .iter()
            .map(|(_, t)| t.clone())
            .collect::<Vec<_>>(),
        signal_genes = data.signal_genes,
        noise_genes = data.noise_genes
    )
}
