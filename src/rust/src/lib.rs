pub mod data;
pub mod embedding;
pub mod graph;
pub mod utils;

use extendr_api::prelude::*;
use faer::Mat;
use rand::{rngs::StdRng, Rng, SeedableRng};
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

/// Generate GeneWalk node embeddings (multiple reps)
///
/// @description Trains node2vec on the original GeneWalk network `n_graph`
/// times with different random seeds, returning a list of embedding matrices.
/// This provides the variance across embedding reps needed for the final
/// log_stats aggregation.
///
/// @param from Integer vector. Node indices for edge origins.
/// @param to Integer vector. Node indices for edge destinations.
/// @param weights Optional numeric vector. Edge weights, defaults to 1.0.
/// @param gene_walk_params Named list. Training parameters (p, q,
///   walks_per_node, walk_length, num_workers, n_epochs, num_negatives,
///   window_size, lr, dim).
/// @param n_graph Integer. Number of independent embedding reps.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Treat graph as directed.
/// @param seed Integer. Random seed (incremented per rep).
/// @param verbose Boolean. Controls verbosity.
///
/// @return A list of n_graph numeric matrices, each n_nodes x embedding dim.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    n_graph: usize,
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
    let edges = prepare_edges(from, to, weights, directed);

    let vocab_size = edges
        .iter()
        .flat_map(|(from, to, _)| [*from, *to])
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0);

    if verbose {
        println!(
            "Processed {} edges ({} nodes) in {:.2}s",
            edges.len(),
            vocab_size,
            start_edge_prep.elapsed().as_secs_f64()
        );
    }

    // build the node2vec graph once (transition probs are deterministic)
    let node2vec_graph = create_graph(&edges, config.p, config.q, directed);

    let mut res_ls = List::new(n_graph);

    let mut master_rng = StdRng::seed_from_u64(seed as u64);

    for rep in 0..n_graph {
        let rep_seed = master_rng.random::<u64>() as usize;
        config.seed = rep_seed;

        if verbose {
            println!("Rep {}/{}: generating walks and training", rep + 1, n_graph);
        }

        let start_rep = Instant::now();

        let random_walks =
            node2vec_graph.generate_walks(config.walks_per_node, config.walk_length, rep_seed);
        let embedding = train_node2vec(random_walks, vocab_size, &config, verbose);

        let nrows = embedding.len();
        let ncols = embedding[0].len();
        let mat = RMatrix::new_matrix(nrows, ncols, |r, c| embedding[r][c] as f64);

        if verbose {
            println!(
                "Rep {}/{} completed in {:.2}s",
                rep + 1,
                n_graph,
                start_rep.elapsed().as_secs_f64()
            );
        }

        res_ls.set_elt(rep, Robj::from(mat))?;
    }

    Ok(res_ls)
}

//////////////////
// Permutations //
//////////////////

/// Generate null distributions from degree-matched random networks
///
/// @description For each permutation, generates a random graph via the
/// configuration model (matching the original degree distribution), trains
/// node2vec on it, then collects cosine similarities between each node and
/// its unique neighbours in the random graph. This matches the original Python
/// GeneWalk procedure.
///
/// @param from Integer vector. Node indices for edge origins.
/// @param to Integer vector. Node indices for edge destinations.
/// @param weights Optional numeric vector. Edge weights, defaults to 1.0.
/// @param gene_walk_params Named list. Training parameters.
/// @param n_perm Integer. Number of null permutations.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Treat graph as directed.
/// @param seed Integer. Random seed.
/// @param verbose Boolean. Controls verbosity.
///
/// @returns A list of n_perm numeric vectors of null cosine similarities.
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

    let from = from.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let edges = prepare_edges(from, to, weights, directed);

    let mut res_ls = List::new(n_perm);

    let mut master_rng = StdRng::seed_from_u64(seed as u64);

    for perm in 0..n_perm {
        let perm_seed = master_rng.random::<u64>();

        if verbose {
            println!("Null permutation {}/{}", perm + 1, n_perm);
        }

        let start_perm = Instant::now();

        // configuration model - degree-matched random graph with new node IDs
        // this caused me such a headache...
        let (rand_edges, rand_vocab_size) =
            generate_configuration_model(&edges, directed, perm_seed);

        // reduce to largest connected component to avoid weirdness with
        // small isolated communities (also important and headache inducing)
        let (rand_edges, rand_vocab_size) =
            largest_connected_component(&rand_edges, rand_vocab_size);

        // train node2vec on the random graph
        let rand_graph = create_graph(&rand_edges, config.p, config.q, directed);
        let rand_walks =
            rand_graph.generate_walks(config.walks_per_node, config.walk_length, config.seed);
        let rand_embedding = train_node2vec(rand_walks, rand_vocab_size, &config, verbose);

        // sample null: for each node, cosine similarity with its unique(!)
        // neighbours in the random graph (matching Python's get_null_distributions)
        let adj = build_adjacency(&rand_edges);
        let mut null_sims: Vec<f64> = Vec::new();

        for (node, neighbours) in &adj {
            let node_idx = *node as usize;
            if node_idx >= rand_embedding.len() {
                continue;
            }
            for &nbr in neighbours {
                let nbr_idx = nbr as usize;
                if nbr_idx >= rand_embedding.len() {
                    continue;
                }
                null_sims.push(cosine_similarity(
                    &rand_embedding[node_idx],
                    &rand_embedding[nbr_idx],
                ) as f64);
            }
        }

        if verbose {
            println!(
                "  {} null similarities from {} random edges in {:.2}s",
                null_sims.len(),
                rand_edges.len(),
                start_perm.elapsed().as_secs_f64()
            );
        }

        res_ls.set_elt(perm, Robj::from(null_sims))?;
    }

    Ok(res_ls)
}

/////////////////////
// Test statistics //
/////////////////////

/// Calculate GeneWalk test statistics
///
/// @description Pools all null distributions into one, then for each observed
/// embedding rep: computes p-values, per-gene FDR, and global FDR. Finally
/// aggregates across reps with geometric mean and 95\% CI (log_stats).
///
/// Matches the paper: "we repeat the above-described network representation
/// learning and significance testing procedures of the GWN nreps_graph times
/// and provide the mean and 95\% confidence intervals".
///
/// @param gene_embds_list List of n_graph gene embedding matrices
///   (n_genes x dim).
/// @param pathway_embds_list List of n_graph pathway embedding matrices
///   (n_pathways x dim).
/// @param null_similarities List of n_perm numeric vectors (null cosine
///   similarities to be pooled).
/// @param connected_pathways List. Per-gene integer vectors of connected
///   pathway indices (1-indexed).
/// @param verbose Controls verbosity.
///
/// @returns A list with per-pair statistics.
///
/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk_test(
    gene_embds_list: List,
    pathway_embds_list: List,
    null_similarities: List,
    connected_pathways: List,
    verbose: bool,
) -> List {
    let n_reps = gene_embds_list.len();
    let n_genes = {
        let first: RMatrix<f64> = gene_embds_list.elt(0).unwrap().try_into().unwrap();
        first.nrows()
    };

    // parse connectivity (1-indexed -> 0-indexed)
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

    // pool all null similarities into one sorted distribution
    let mut pooled_null: Vec<f64> = Vec::new();
    for (_, perm_item) in null_similarities.iter() {
        if let Some(v) = perm_item.as_real_vector() {
            pooled_null.extend(v);
        }
    }
    pooled_null.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if verbose {
        let total_pairs: usize = connectivity.iter().map(|c| c.len()).sum();
        println!(
            "Pooled null: {} similarities, testing {} connected pairs across {} reps",
            pooled_null.len(),
            total_pairs,
            n_reps
        );
    }

    // compute cosine similarities, p-values, and FDR per rep
    let mut raw_pvals_per_rep: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_reps);
    let mut gene_fdr_per_rep: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_reps);
    let mut global_fdr_per_rep: Vec<Vec<f64>> = Vec::with_capacity(n_reps);
    let mut cosine_sims_per_rep: Vec<Mat<f64>> = Vec::with_capacity(n_reps);

    for rep in 0..n_reps {
        if verbose {
            println!("Processing rep {}/{}", rep + 1, n_reps);
        }

        let gene_mat: RMatrix<f64> = gene_embds_list.elt(rep).unwrap().try_into().unwrap();
        let pathway_mat: RMatrix<f64> = pathway_embds_list.elt(rep).unwrap().try_into().unwrap();

        let gene_embds = r_matrix_to_vec(gene_mat);
        let pathway_embds = r_matrix_to_vec(pathway_mat);

        let cosine_sim: Mat<f64> = compute_cross_cosine(&gene_embds, &pathway_embds);

        // P-values for connected pairs against pooled null
        let gene_pvals: Vec<Vec<f64>> = (0..n_genes)
            .map(|gi| {
                connectivity[gi]
                    .iter()
                    .map(|&pi| {
                        let obs = cosine_sim[(gi, pi)];
                        let rank = pooled_null.partition_point(|&x| x < obs);
                        let pval = 1.0 - (rank as f64 / pooled_null.len() as f64);
                        pval.max(1e-16)
                    })
                    .collect()
            })
            .collect();

        // global FDR: across all connected pairs
        let all_connected_pvals: Vec<f64> = gene_pvals.iter().flatten().cloned().collect();
        let global_fdr_flat = calc_fdr(&all_connected_pvals);

        // per-gene FDR
        let gene_fdr: Vec<Vec<f64>> = gene_pvals.iter().map(|p| calc_fdr(p)).collect();

        raw_pvals_per_rep.push(gene_pvals);
        gene_fdr_per_rep.push(gene_fdr);
        global_fdr_per_rep.push(global_fdr_flat);
        cosine_sims_per_rep.push(cosine_sim);
    }

    // aggregate across reps

    let mut out_gene: Vec<i32> = Vec::new();
    let mut out_pathway: Vec<i32> = Vec::new();
    let mut out_similarity: Vec<f64> = Vec::new();
    let mut out_sem_sim: Vec<f64> = Vec::new();
    let mut out_avg_pval: Vec<f64> = Vec::new();
    let mut out_pval_ci_lower: Vec<f64> = Vec::new();
    let mut out_pval_ci_upper: Vec<f64> = Vec::new();
    let mut out_avg_global_fdr: Vec<f64> = Vec::new();
    let mut out_global_fdr_ci_lower: Vec<f64> = Vec::new();
    let mut out_global_fdr_ci_upper: Vec<f64> = Vec::new();
    let mut out_avg_gene_fdr: Vec<f64> = Vec::new();
    let mut out_gene_fdr_ci_lower: Vec<f64> = Vec::new();
    let mut out_gene_fdr_ci_upper: Vec<f64> = Vec::new();

    let mut global_idx = 0usize;

    #[allow(clippy::needless_range_loop)]
    for gene_i in 0..n_genes {
        for (local_idx, &pathway_i) in connectivity[gene_i].iter().enumerate() {
            // aimilarity: mean and SEM across reps
            let sims: Vec<f64> = (0..n_reps)
                .map(|r| cosine_sims_per_rep[r][(gene_i, pathway_i)])
                .collect();
            let mean_sim = sims.iter().sum::<f64>() / n_reps as f64;
            let sem_sim = if n_reps > 1 {
                let var =
                    sims.iter().map(|&s| (s - mean_sim).powi(2)).sum::<f64>() / (n_reps - 1) as f64;
                var.sqrt() / (n_reps as f64).sqrt()
            } else {
                0.0
            };

            // pvals across reps
            let raw_pvals: Vec<f64> = (0..n_reps)
                .map(|r| raw_pvals_per_rep[r][gene_i][local_idx])
                .collect();
            let gene_fdrs: Vec<f64> = (0..n_reps)
                .map(|r| gene_fdr_per_rep[r][gene_i][local_idx])
                .collect();
            let global_fdrs: Vec<f64> = (0..n_reps)
                .map(|r| global_fdr_per_rep[r][global_idx])
                .collect();

            let (avg_pval, pval_lo, pval_hi) = log_stats(&raw_pvals);
            let (avg_gfdr, gfdr_lo, gfdr_hi) = log_stats(&global_fdrs);
            let (avg_gene_fdr, gene_lo, gene_hi) = log_stats(&gene_fdrs);

            out_gene.push((gene_i + 1) as i32);
            out_pathway.push((pathway_i + 1) as i32);
            out_similarity.push(mean_sim);
            out_sem_sim.push(sem_sim);
            out_avg_pval.push(avg_pval);
            out_pval_ci_lower.push(pval_lo);
            out_pval_ci_upper.push(pval_hi);
            out_avg_global_fdr.push(avg_gfdr);
            out_global_fdr_ci_lower.push(gfdr_lo);
            out_global_fdr_ci_upper.push(gfdr_hi);
            out_avg_gene_fdr.push(avg_gene_fdr);
            out_gene_fdr_ci_lower.push(gene_lo);
            out_gene_fdr_ci_upper.push(gene_hi);

            global_idx += 1;
        }
    }

    list![
        gene = out_gene,
        pathway = out_pathway,
        similarity = out_similarity,
        sem_sim = out_sem_sim,
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

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::graph::*;
    use super::utils::*;
    use rustc_hash::FxHashMap;

    /// Helper: compute degree sequence from an edge list (undirected)
    fn degree_sequence(edges: &[(u32, u32, f32)]) -> Vec<usize> {
        let mut deg: FxHashMap<u32, usize> = FxHashMap::default();
        for (u, v, _) in edges {
            *deg.entry(*u).or_insert(0) += 1;
            *deg.entry(*v).or_insert(0) += 1;
        }
        let mut seq: Vec<usize> = deg.values().cloned().collect();
        seq.sort_unstable_by(|a, b| b.cmp(a));
        seq
    }

    /// Helper: check no self-loops exist
    fn has_self_loops(edges: &[(u32, u32, f32)]) -> bool {
        edges.iter().any(|(u, v, _)| u == v)
    }

    /// Helper: check no duplicate edges exist (undirected)
    fn has_duplicate_edges(edges: &[(u32, u32, f32)]) -> bool {
        let mut seen = std::collections::HashSet::new();
        for (u, v, _) in edges {
            let key = if u <= v { (*u, *v) } else { (*v, *u) };
            if !seen.insert(key) {
                return true;
            }
        }
        false
    }

    /// Helper: count connected components via union-find
    fn count_components(edges: &[(u32, u32, f32)], n_nodes: usize) -> usize {
        let mut parent: Vec<usize> = (0..n_nodes).collect();
        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }
        for (u, v, _) in edges {
            let pu = find(&mut parent, *u as usize);
            let pv = find(&mut parent, *v as usize);
            if pu != pv {
                parent[pu] = pv;
            }
        }
        let mut roots = std::collections::HashSet::new();
        for i in 0..n_nodes {
            roots.insert(find(&mut parent, i));
        }
        roots.len()
    }

    // -------------------------------------------------------
    // Configuration model: structural invariants
    // -------------------------------------------------------

    #[test]
    fn config_model_no_self_loops() {
        let edges = vec![
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
            (0, 2, 1.0),
            (1, 3, 1.0),
        ];
        for seed in 0..20 {
            let (rand_edges, _) = generate_configuration_model(&edges, false, seed);
            assert!(
                !has_self_loops(&rand_edges),
                "Self-loop found with seed {seed}"
            );
        }
    }

    #[test]
    fn config_model_no_multi_edges() {
        let edges = vec![
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
            (0, 2, 1.0),
            (1, 3, 1.0),
        ];
        for seed in 0..20 {
            let (rand_edges, _) = generate_configuration_model(&edges, false, seed);
            assert!(
                !has_duplicate_edges(&rand_edges),
                "Duplicate edge found with seed {seed}"
            );
        }
    }

    #[test]
    fn config_model_preserves_degree_distribution_small() {
        // Small dense graphs lose more edges to collisions -- that's expected
        let edges = vec![
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 4, 1.0),
            (4, 0, 1.0),
            (0, 3, 1.0),
            (1, 4, 1.0),
            (2, 4, 1.0),
        ];
        let orig_deg = degree_sequence(&edges);
        let orig_sum: usize = orig_deg.iter().sum();

        for seed in 0..10 {
            let (rand_edges, _) = generate_configuration_model(&edges, false, seed);
            let rand_deg = degree_sequence(&rand_edges);
            let rand_sum: usize = rand_deg.iter().sum();

            let loss = 1.0 - (rand_sum as f64 / orig_sum as f64);
            assert!(
                loss <= 0.5,
                "seed {seed}: lost {:.0}% of degree sum ({rand_sum} vs {orig_sum})",
                loss * 100.0
            );
        }
    }

    #[test]
    fn config_model_preserves_degree_distribution_realistic() {
        // Larger sparser graph closer to real GeneWalk networks
        // 50 nodes, ~100 edges -- collision rate should be low
        let mut edges = Vec::new();

        // Hub with 20 connections
        for i in 1..21u32 {
            edges.push((0, i, 1.0));
        }
        // Chain through the rest
        for i in 1..49u32 {
            edges.push((i, i + 1, 1.0));
        }
        // A few cross-links
        for i in (5..45).step_by(10) {
            edges.push((i, i + 3, 1.0));
        }

        let orig_deg = degree_sequence(&edges);
        let orig_sum: usize = orig_deg.iter().sum();

        for seed in 0..20 {
            let (rand_edges, _) = generate_configuration_model(&edges, false, seed);
            let rand_deg = degree_sequence(&rand_edges);
            let rand_sum: usize = rand_deg.iter().sum();

            let loss = 1.0 - (rand_sum as f64 / orig_sum as f64);
            assert!(
                loss < 0.15,
                "seed {seed}: lost {:.0}% of degree sum ({rand_sum} vs {orig_sum})",
                loss * 100.0
            );
        }
    }

    #[test]
    fn config_model_node_count_matches() {
        let edges = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0)];
        for seed in 0..10 {
            let (_, vocab_size) = generate_configuration_model(&edges, false, seed);
            assert_eq!(vocab_size, 5, "Expected 5 nodes, got {vocab_size}");
        }
    }

    #[test]
    fn config_model_different_seeds_give_different_graphs() {
        let edges: Vec<(u32, u32, f32)> = (0..20)
            .flat_map(|i| (i + 1..20).map(move |j| (i, j, 1.0)))
            .take(50)
            .collect();

        let (e1, _) = generate_configuration_model(&edges, false, 42);
        let (e2, _) = generate_configuration_model(&edges, false, 99);

        // Not a hard guarantee but with 50 edges, extremely unlikely to match
        assert_ne!(e1.len(), 0);
        assert_ne!(e1, e2, "Different seeds produced identical graphs");
    }

    // -------------------------------------------------------
    // Configuration model: fragmentation check
    // -------------------------------------------------------

    #[test]
    fn config_model_fragmentation_report() {
        // Use a graph that resembles a real GeneWalk network: bipartite-ish
        // with hubs. This test doesn't assert a threshold but prints stats
        // so you can see how bad fragmentation gets.
        let mut edges = Vec::new();
        // Hub node 0 connects to 1..20
        for i in 1..20u32 {
            edges.push((0, i, 1.0));
        }
        // Chain among 1..20
        for i in 1..19u32 {
            edges.push((i, i + 1, 1.0));
        }
        // Some cross-links
        edges.push((5, 15, 1.0));
        edges.push((3, 17, 1.0));
        edges.push((7, 12, 1.0));

        let n_nodes = 20;
        let orig_components = count_components(&edges, n_nodes);

        let mut max_components = 0usize;
        let mut total_components = 0usize;
        let n_trials = 50;

        for seed in 0..n_trials {
            let (rand_edges, vocab) = generate_configuration_model(&edges, false, seed as u64);
            let nc = count_components(&rand_edges, vocab);
            total_components += nc;
            if nc > max_components {
                max_components = nc;
            }
        }

        let avg = total_components as f64 / n_trials as f64;
        println!("Original: {orig_components} component(s)");
        println!("Config model over {n_trials} seeds: avg={avg:.1}, max={max_components}");

        // Soft assertion: average shouldn't be absurd
        assert!(
            avg < n_nodes as f64 / 2.0,
            "Average component count {avg:.1} is too high for {n_nodes} nodes"
        );
    }

    // -------------------------------------------------------
    // Adjacency builder
    // -------------------------------------------------------

    #[test]
    fn adjacency_excludes_self_loops() {
        let edges = vec![(0, 0, 1.0), (0, 1, 1.0), (1, 2, 1.0)];
        let adj = build_adjacency(&edges);
        // Node 0 should only have neighbour 1, not itself
        assert!(!adj.get(&0).unwrap().contains(&0));
        assert!(adj.get(&0).unwrap().contains(&1));
    }

    #[test]
    fn adjacency_is_symmetric() {
        let edges = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)];
        let adj = build_adjacency(&edges);
        for (node, neighbours) in &adj {
            for nbr in neighbours {
                assert!(
                    adj.get(nbr).unwrap().contains(node),
                    "Edge ({node}, {nbr}) not symmetric"
                );
            }
        }
    }

    #[test]
    fn adjacency_neighbour_count_matches_degree() {
        let edges = vec![(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0)];
        let adj = build_adjacency(&edges);
        assert_eq!(adj.get(&0).unwrap().len(), 3); // connected to 1, 2, 3
        assert_eq!(adj.get(&1).unwrap().len(), 2); // connected to 0, 2
        assert_eq!(adj.get(&3).unwrap().len(), 1); // connected to 0 only
    }

    // -------------------------------------------------------
    // Null distribution: sanity checks
    // -------------------------------------------------------

    #[test]
    fn null_similarities_are_bounded() {
        // Cosine similarity must be in [-1, 1]
        // Simulate what rs_gene_walk_perm does: fake embeddings, compute sims
        let embeddings: Vec<Vec<f32>> = (0..10)
            .map(|i| (0..8).map(|j| ((i * 7 + j) as f32).sin()).collect())
            .collect();

        let edges = vec![
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 4, 1.0),
            (4, 5, 1.0),
        ];
        let adj = build_adjacency(&edges);

        for (node, neighbours) in &adj {
            for nbr in neighbours {
                let sim = cosine_similarity(&embeddings[*node as usize], &embeddings[*nbr as usize])
                    as f64;
                assert!(
                    (-1.0..=1.0).contains(&sim),
                    "Cosine similarity {sim} out of bounds for ({node}, {nbr})"
                );
            }
        }
    }

    #[test]
    fn null_from_adjacency_has_expected_count() {
        // Each undirected edge contributes two similarities (once from each
        // endpoint's perspective), so total null samples = 2 * n_edges
        let edges = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0)];
        let adj = build_adjacency(&edges);

        let mut count = 0usize;
        for neighbours in adj.values() {
            count += neighbours.len();
        }

        assert_eq!(
            count,
            2 * edges.len(),
            "Expected {} null samples, got {count}",
            2 * edges.len()
        );
    }

    // -------------------------------------------------------
    // P-value calculation sanity
    // -------------------------------------------------------

    #[test]
    fn pval_of_max_null_is_near_zero() {
        let null: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        let obs = 0.999; // higher than almost everything
        let rank = null.partition_point(|&x| x < obs);
        let pval = 1.0 - (rank as f64 / null.len() as f64);
        assert!(pval < 0.01, "p-value {pval} should be near zero");
    }

    #[test]
    fn pval_of_min_null_is_near_one() {
        let null: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        let obs = 0.001;
        let rank = null.partition_point(|&x| x < obs);
        let pval = 1.0 - (rank as f64 / null.len() as f64);
        assert!(pval > 0.99, "p-value {pval} should be near one");
    }

    #[test]
    fn pval_of_median_is_near_half() {
        let null: Vec<f64> = (0..1000).map(|i| i as f64 / 1000.0).collect();
        let obs = 0.5;
        let rank = null.partition_point(|&x| x < obs);
        let pval = 1.0 - (rank as f64 / null.len() as f64);
        assert!((pval - 0.5).abs() < 0.05, "p-value {pval} should be ~0.5");
    }
}
