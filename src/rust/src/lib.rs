pub mod data;
pub mod embedding;
pub mod graph;
pub mod utils;

use burn::prelude::Backend;
use extendr_api::prelude::*;
use faer::Mat;
use itertools::Itertools;
use node2vec_rs::model::SkipGramConfig;
use rayon::prelude::*;
use std::time::Instant;

// torch is only support on NoneWindows
#[cfg(not(target_os = "windows"))]
use burn::backend::libtorch::{LibTorch, LibTorchDevice};

use burn::backend::{
    ndarray::{NdArray, NdArrayDevice},
    Autodiff,
};

use crate::data::*;
use crate::embedding::*;
use crate::graph::*;
use crate::utils::*;

// export the function to R
extendr_module! {
    mod genewalkR;
    fn rs_gene_walk;
    fn rs_gene_walk_perm;
    fn rs_gene_walk_test;
    fn rs_cosine_sim;
    fn rs_node2vec_synthetic_data;
}

/////////////////////////
// Gene Walk functions //
/////////////////////////

/// Function that generates the GeneWalk embedding
///
/// @description Wrapper function that leverages Rust and the Burn Tensor
/// framework to generate node representations for subsequent usage in the
/// GeneWalk approach. The default back end is the torch CPU backend. In the
/// future there might be other backends added, but testing revealed this one
/// to be very performant.
///
/// @param from Integer vector. The node indices indicating where the edge
/// originates from.
/// @param to Integer vector. The node indices indicating where the edge
/// goes to.
/// @param weights Optional numeric vector. If not supplied, defaults to `1.0`
/// as edge weight.
/// @param gene_walk_params Named list. Contains the parameters for running
/// gene walk.
/// @param backend String. One of `c("tch-cpu", "ndarray")`. Torch is not
/// supported on Windows and will panic.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @return A numerical matrix of number nodes x desired embedding dimensios.
///
/// @export
#[cfg(not(target_os = "windows"))]
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    backend: String,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> RMatrix<f64> {
    let gene_walk_config = GeneWalkConfig::from_r_list(gene_walk_params, seed);

    if verbose {
        println!("Preparing the network edges.")
    }
    let start_edge_prep = Instant::now();

    let from = from.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();

    let edges = prepare_edges(from, to, weights);

    let end_edge_prep = start_edge_prep.elapsed();

    if verbose {
        println!("Processed the edges in {:.2}s", end_edge_prep.as_secs_f64());
        println!("Running the random walks.")
    }

    let start_rw = Instant::now();

    let node2vec_graph = create_graph(&edges, gene_walk_config.p, gene_walk_config.q, directed);

    let vocab_size = edges
        .iter()
        .flat_map(|(from, to, _)| [*from, *to])
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0);

    let random_walks = node2vec_graph.generate_walks(
        gene_walk_config.walks_per_node,
        gene_walk_config.walk_length,
        gene_walk_config.seed,
    );

    let end_rw = start_rw.elapsed();

    if verbose {
        println!("Generated the random walks in {:.2}s", end_rw.as_secs_f64());
        println!(
            "Training the SkipGram model now for {:?} epochs",
            gene_walk_config.num_epochs
        );
    }

    let backend = parse_backend(&backend).unwrap_or_default();
    let start_training = Instant::now();

    let embedding: Vec<Vec<f32>> = match backend {
        Node2VecBackEnd::TorchCpu => {
            let device = LibTorchDevice::Cpu;
            LibTorch::<f32>::seed(&device, gene_walk_config.seed);

            let skipgram_config = SkipGramConfig {
                vocab_size,
                embedding_dim: embd_dim,
            };

            train_node2vec::<Autodiff<LibTorch>>(
                skipgram_config,
                gene_walk_config,
                random_walks,
                device,
                verbose,
            )
        }
        Node2VecBackEnd::Ndarray => {
            let device = NdArrayDevice::Cpu;
            NdArray::<f32>::seed(&device, gene_walk_config.seed);

            let skipgram_config = SkipGramConfig {
                vocab_size,
                embedding_dim: embd_dim,
            };

            train_node2vec::<Autodiff<NdArray>>(
                skipgram_config,
                gene_walk_config,
                random_walks,
                device,
                verbose,
            )
        }
    };

    let end_training = start_training.elapsed();

    if verbose {
        println!("Training completed in {:.2}s.", end_training.as_secs_f64());
    }

    let nrows = embedding.len();
    let ncols = embedding[0].len();

    RMatrix::new_matrix(nrows, ncols, |r, c| embedding[r][c] as f64)
}

/// Function that generates the GeneWalk embedding
///
/// @description Wrapper function that leverages Rust and the Burn Tensor
/// framework to generate node representations for subsequent usage in the
/// GeneWalk approach. The default back end is the torch CPU backend. In the
/// future there might be other backends added, but testing revealed this one
/// to be very performant.
///
/// @param from Integer vector. The node indices indicating where the edge
/// originates from.
/// @param to Integer vector. The node indices indicating where the edge
/// goes to.
/// @param weights Optional numeric vector. If not supplied, defaults to `1.0`
/// as edge weight.
/// @param gene_walk_params Named list. Contains the parameters for running
/// gene walk.
/// @param backend String. One of `c("tch-cpu", "ndarray")`. Torch is not
/// supported on Windows and will panic.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @return A numerical matrix of number nodes x desired embedding dimensios.
///
/// @export
#[cfg(target_os = "windows")]
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    backend: String,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> RMatrix<f64> {
    let gene_walk_config = GeneWalkConfig::from_r_list(gene_walk_params, seed);

    if verbose {
        println!("Preparing the network edges.")
    }
    let start_edge_prep = Instant::now();

    let from = from.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| (*x - 1) as u32).collect::<Vec<u32>>();

    let edges = prepare_edges(from, to, weights);

    let end_edge_prep = start_edge_prep.elapsed();

    if verbose {
        println!("Processed the edges in {:.2}s", end_edge_prep.as_secs_f64());
        println!("Running the random walks.")
    }

    let start_rw = Instant::now();

    let node2vec_graph = create_graph(&edges, gene_walk_config.p, gene_walk_config.q, directed);

    let vocab_size = edges
        .iter()
        .flat_map(|(from, to, _)| [*from, *to])
        .max()
        .map(|max_id| max_id as usize + 1)
        .unwrap_or(0);

    let random_walks = node2vec_graph.generate_walks(
        gene_walk_config.walks_per_node,
        gene_walk_config.walk_length,
        gene_walk_config.seed,
    );

    let end_rw = start_rw.elapsed();

    if verbose {
        println!("Generated the random walks in {:.2}s", end_rw.as_secs_f64());
        println!(
            "Training the SkipGram model now for {:?} epochs",
            gene_walk_config.num_epochs
        );
    }

    let backend = parse_backend(&backend).unwrap_or_default();
    let start_training = Instant::now();

    let embedding: Vec<Vec<f32>> = match backend {
        Node2VecBackEnd::TorchCpu => {
            panic!("Torch is not supported on Windows!")
        }
        Node2VecBackEnd::Ndarray => {
            let device = NdArrayDevice::Cpu;
            NdArray::<f32>::seed(&device, gene_walk_config.seed);

            let skipgram_config = SkipGramConfig {
                vocab_size,
                embedding_dim: embd_dim,
            };

            train_node2vec::<Autodiff<NdArray>>(
                skipgram_config,
                gene_walk_config,
                random_walks,
                device,
                verbose,
            )
        }
    };

    let end_training = start_training.elapsed();

    if verbose {
        println!("Training completed in {:.2}s.", end_training.as_secs_f64());
    }

    let nrows = embedding.len();
    let ncols = embedding[0].len();

    RMatrix::new_matrix(nrows, ncols, |r, c| embedding[r][c] as f64)
}

//////////////////
// Permutations //
//////////////////

/// Function to generate permuted embeddings
///
/// @description
/// Wrapper function to generate permuted embeddings and subsequently the
/// null distribution of cosine similarities for statistical testing.
///
/// @param from Integer vector. The node indices indicating where the edge
/// originates from.
/// @param to Integer vector. The node indices indicating where the edge
/// goes to.
/// @param weights Optional numeric vector. If not supplied, defaults to `1.0`
/// as edge weight.
/// @param gene_walk_params Named list. Contains the parameters for running
/// gene walk.
/// @param backend String. One of `c("tch-cpu", "ndarray")`. Torch is not
/// supported on Windows and will panic.
/// @param n_perm Integer.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @returns A list with the null distribution of the cosine similarities per
/// given permutation.
///
/// @export
#[extendr]
#[cfg(not(target_os = "windows"))]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk_perm(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    backend: String,
    n_perm: usize,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let gene_walk_config = GeneWalkConfig::from_r_list(gene_walk_params, seed);

    if verbose {
        println!("Preparing the network edges.")
    }
    let start_edge_prep = Instant::now();

    // go from R i32 to u32
    let from = from.iter().map(|x| *x as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| *x as u32).collect::<Vec<u32>>();

    let edges = prepare_edges(from, to, weights);

    let end_edge_prep = start_edge_prep.elapsed();
    if verbose {
        println!("Processed the edges in {:.2}s", end_edge_prep.as_secs_f64());
    }

    let mut null_similarities: Vec<Vec<f64>> = Vec::with_capacity(n_perm);

    let backend = parse_backend(&backend).unwrap_or_default();

    for perm in 0..n_perm {
        let start_perm = Instant::now();

        let gene_walk_config_i = gene_walk_config.clone();

        if verbose {
            println!("Starting permutation: {:?}", perm + 1);
            println!(" Generating randomised edges and random walks.")
        }

        let start_rw = Instant::now();

        let perm_edges = generate_random_network(&edges, directed, (seed + perm) as u64);

        let node2vec_graph = create_graph(
            &perm_edges,
            gene_walk_config_i.p,
            gene_walk_config_i.q,
            directed,
        );

        let vocab_size = edges
            .iter()
            .flat_map(|(from, to, _)| [*from, *to])
            .max()
            .map(|max_id| max_id as usize + 1)
            .unwrap_or(0);

        let random_walks_perm = node2vec_graph.generate_walks(
            gene_walk_config_i.walks_per_node,
            gene_walk_config_i.walk_length,
            gene_walk_config_i.seed,
        );

        let end_rw = start_rw.elapsed();

        if verbose {
            println!(
                " Finished the generation of the random walks on permutated data in {:.2}s",
                end_rw.as_secs_f64()
            );
            println!(
                " Training the SkipGram model on permuted network for {:?} epochs",
                gene_walk_config.num_epochs
            );
        }

        let skipgram_config = SkipGramConfig {
            vocab_size,
            embedding_dim: embd_dim,
        };

        let embd_i: Vec<Vec<f32>> = match backend {
            Node2VecBackEnd::TorchCpu => {
                let device = LibTorchDevice::Cpu;
                LibTorch::<f32>::seed(&device, gene_walk_config.seed);

                train_node2vec::<Autodiff<LibTorch>>(
                    skipgram_config,
                    gene_walk_config_i,
                    random_walks_perm,
                    device,
                    verbose,
                )
            }
            Node2VecBackEnd::Ndarray => {
                let device = NdArrayDevice::Cpu;
                NdArray::<f32>::seed(&device, gene_walk_config.seed);

                let skipgram_config = SkipGramConfig {
                    vocab_size,
                    embedding_dim: embd_dim,
                };

                train_node2vec::<Autodiff<NdArray>>(
                    skipgram_config,
                    gene_walk_config_i,
                    random_walks_perm,
                    device,
                    verbose,
                )
            }
        };

        let end_perm = start_perm.elapsed();

        let null_vals_i = perm_edges
            .par_iter()
            .map(|(node_i, node_j, _)| {
                cosine_similarity(&embd_i[*node_i as usize], &embd_i[*node_j as usize]) as f64
            })
            .collect::<Vec<f64>>();

        null_similarities.push(null_vals_i);

        if verbose {
            println!(
                " Finished the generation of the embedding on permutated data in {:.2}s",
                end_perm.as_secs_f64()
            );
        }
    }

    let mut res_ls = List::new(n_perm);

    #[allow(clippy::needless_range_loop)]
    for i in 0..n_perm {
        res_ls.set_elt(i, Robj::from(&null_similarities[i]))?;
    }

    Ok(res_ls)
}

/// Function to generate permuted embeddings
///
/// @description
/// Wrapper function to generate permuted embeddings and subsequently the
/// null distribution of cosine similarities for statistical testing.
///
/// @param from Integer vector. The node indices indicating where the edge
/// originates from.
/// @param to Integer vector. The node indices indicating where the edge
/// goes to.
/// @param weights Optional numeric vector. If not supplied, defaults to `1.0`
/// as edge weight.
/// @param gene_walk_params Named list. Contains the parameters for running
/// gene walk.
/// @param backend String. One of `c("tch-cpu", "ndarray")`. Torch is not
/// supported on Windows and will panic.
/// @param n_perm Integer.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @returns A list with the null distribution of the cosine similarities per
/// given permutation.
///
/// @export
#[cfg(target_os = "windows")]
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_gene_walk_perm(
    from: Vec<i32>,
    to: Vec<i32>,
    weights: Option<Vec<f64>>,
    gene_walk_params: List,
    backend: String,
    n_perm: usize,
    embd_dim: usize,
    directed: bool,
    seed: usize,
    verbose: bool,
) -> extendr_api::Result<List> {
    let gene_walk_config = GeneWalkConfig::from_r_list(gene_walk_params, seed);

    if verbose {
        println!("Preparing the network edges.")
    }
    let start_edge_prep = Instant::now();

    // go from R i32 to u32
    let from = from.iter().map(|x| *x as u32).collect::<Vec<u32>>();
    let to = to.iter().map(|x| *x as u32).collect::<Vec<u32>>();

    let edges = prepare_edges(from, to, weights);

    let end_edge_prep = start_edge_prep.elapsed();
    if verbose {
        println!("Processed the edges in {:.2}s", end_edge_prep.as_secs_f64());
    }

    let mut null_similarities: Vec<Vec<f64>> = Vec::with_capacity(n_perm);

    let backend = parse_backend(&backend).unwrap_or_default();

    for perm in 0..n_perm {
        let start_perm = Instant::now();

        let gene_walk_config_i = gene_walk_config.clone();

        if verbose {
            println!("Starting permutation: {:?}", perm + 1);
            println!(" Generating randomised edges and random walks.")
        }

        let start_rw = Instant::now();

        let perm_edges = generate_random_network(&edges, directed, (seed + perm) as u64);

        let node2vec_graph = create_graph(
            &perm_edges,
            gene_walk_config_i.p,
            gene_walk_config_i.q,
            directed,
        );

        let vocab_size = edges
            .iter()
            .flat_map(|(from, to, _)| [*from, *to])
            .max()
            .map(|max_id| max_id as usize + 1)
            .unwrap_or(0);

        let random_walks_perm = node2vec_graph.generate_walks(
            gene_walk_config_i.walks_per_node,
            gene_walk_config_i.walk_length,
            gene_walk_config_i.seed,
        );

        let end_rw = start_rw.elapsed();

        if verbose {
            println!(
                " Finished the generation of the random walks on permutated data in {:.2}s",
                end_rw.as_secs_f64()
            );
            println!(
                " Training the SkipGram model on permuted network for {:?} epochs",
                gene_walk_config.num_epochs
            );
        }

        let skipgram_config = SkipGramConfig {
            vocab_size,
            embedding_dim: embd_dim,
        };

        let embd_i: Vec<Vec<f32>> = match backend {
            Node2VecBackEnd::TorchCpu => {
                panic!("Windows does not support Torch CPU!")
            }
            Node2VecBackEnd::Ndarray => {
                let device = NdArrayDevice::Cpu;
                NdArray::<f32>::seed(&device, gene_walk_config.seed);

                let skipgram_config = SkipGramConfig {
                    vocab_size,
                    embedding_dim: embd_dim,
                };

                train_node2vec::<Autodiff<NdArray>>(
                    skipgram_config,
                    gene_walk_config_i,
                    random_walks_perm,
                    device,
                    verbose,
                )
            }
        };

        let end_perm = start_perm.elapsed();

        let null_vals_i = perm_edges
            .par_iter()
            .map(|(node_i, node_j, _)| {
                cosine_similarity(&embd_i[*node_i as usize], &embd_i[*node_j as usize]) as f64
            })
            .collect::<Vec<f64>>();

        null_similarities.push(null_vals_i);

        if verbose {
            println!(
                " Finished the generation of the embedding on permutated data in {:.2}s",
                end_perm.as_secs_f64()
            );
        }
    }

    let mut res_ls = List::new(n_perm);

    #[allow(clippy::needless_range_loop)]
    for i in 0..n_perm {
        res_ls.set_elt(i, Robj::from(&null_similarities[i]))?;
    }

    Ok(res_ls)
}

/////////////////////
// Test statistics //
/////////////////////

/// Calculate the test statistics
///
/// @description Calculates the test statistic for the gene / pathway pairs.
///
/// @param gene_embds Matrix of n_genes x their graph embeddings
/// @param pathway_embds Matrix of n_pathways x their graph embeddings
/// @param null_distributions List of null distribution vectors (one per
/// permutation).
/// @param verbose Controls verbosity of the function.
///
/// @returns A list with vectors:
/// \itemize{
///   \item gene - Gene indices (1-based for R)
///   \item pathway - Pathway indices (1-based for R)
///   \item similarity - Cosine similarity between gene and pathway
///   \item avg_pval - Mean p-value across permutations
///   \item pval_ci_lower - Lower 95% CI for p-value
///   \item pval_ci_upper - Upper 95% CI for p-value
///   \item avg_global_fdr - Mean global FDR across permutations
///   \item global_fdr_ci_lower - Lower 95% CI for global FDR
///   \item global_fdr_ci_upper - Upper 95% CI for global FDR
///   \item avg_gene_fdr - Mean gene-specific FDR across permutations
///   \item gene_fdr_ci_lower - Lower 95% CI for gene FDR
///   \item gene_fdr_ci_upper - Upper 95% CI for gene FDR
/// }
///
/// @export
#[extendr]
fn rs_gene_walk_test(
    gene_embds: RMatrix<f64>,
    pathway_embds: RMatrix<f64>,
    null_distributions: List,
    verbose: bool,
) -> List {
    let gene_embds = r_matrix_to_vec(gene_embds);
    let pathway_embds = r_matrix_to_vec(pathway_embds);

    if verbose {
        println!("Calculating Cosine similarities");
    }
    let cosine_sim: Mat<f64> = compute_cross_cosine(&gene_embds, &pathway_embds);

    let n_genes = cosine_sim.nrows();
    let n_pathways = cosine_sim.ncols();
    let total_size = n_genes * n_pathways;
    let n_perms = null_distributions.len();

    if verbose {
        println!("Processing {} permutations", n_perms);
    }

    // Store results from all permutations
    let mut all_pvals = Vec::with_capacity(n_perms);
    let mut all_global_fdr = Vec::with_capacity(n_perms);
    let mut all_gene_fdr = Vec::with_capacity(n_perms);

    for (perm_idx, null_dist) in null_distributions.iter().enumerate() {
        if verbose {
            println!("Processing permutation {}/{}", perm_idx + 1, n_perms);
        }

        let null_vec: Vec<f64> = null_dist.1.as_real_vector().unwrap();
        let pvals: Mat<f64> = calculate_p_vals(&cosine_sim.as_ref(), null_vec);

        // Flatten to vector
        let mut pval_vec = Vec::with_capacity(total_size);
        for i in 0..n_genes {
            for j in 0..n_pathways {
                pval_vec.push(pvals[(i, j)]);
            }
        }

        // Calculate FDRs
        let global_fdr = calc_fdr(&pval_vec);
        let gene_fdr: Vec<f64> = (0..n_genes)
            .into_par_iter()
            .flat_map(|gene_idx| {
                let start_idx = gene_idx * n_pathways;
                let end_idx = start_idx + n_pathways;
                calc_fdr(&pval_vec[start_idx..end_idx])
            })
            .collect();

        all_pvals.push(pval_vec);
        all_global_fdr.push(global_fdr);
        all_gene_fdr.push(gene_fdr);
    }

    if verbose {
        println!("Aggregating results across permutations");
    }

    let results: Vec<_> = (0..total_size)
        .into_par_iter()
        .map(|idx| {
            let gene_idx = idx / n_pathways;
            let pathway_idx = idx % n_pathways;

            // Get the actual cosine similarity
            let similarity = cosine_sim[(gene_idx, pathway_idx)];

            // Collect values across permutations
            let mut pvals_for_pair: Vec<f64> = all_pvals.iter().map(|v| v[idx]).collect();
            let mut global_fdr_for_pair: Vec<f64> = all_global_fdr.iter().map(|v| v[idx]).collect();
            let mut gene_fdr_for_pair: Vec<f64> = all_gene_fdr.iter().map(|v| v[idx]).collect();

            // Unstable sort (faster, order of equal elements doesn't matter)
            pvals_for_pair.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            global_fdr_for_pair.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            gene_fdr_for_pair.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            let n_perms_f64 = n_perms as f64;

            (
                (gene_idx + 1) as i32,
                (pathway_idx + 1) as i32,
                similarity,
                pvals_for_pair.iter().sum::<f64>() / n_perms_f64,
                quantile(&pvals_for_pair, 0.025),
                quantile(&pvals_for_pair, 0.975),
                global_fdr_for_pair.iter().sum::<f64>() / n_perms_f64,
                quantile(&global_fdr_for_pair, 0.025),
                quantile(&global_fdr_for_pair, 0.975),
                gene_fdr_for_pair.iter().sum::<f64>() / n_perms_f64,
                quantile(&gene_fdr_for_pair, 0.025),
                quantile(&gene_fdr_for_pair, 0.975),
            )
        })
        .collect();

    // Unpack parallel results
    #[allow(clippy::complexity)]
    let (
        gene_indices,
        pathway_indices,
        similarities,
        avg_pval,
        pval_ci_lower,
        pval_ci_upper,
        avg_global_fdr,
        global_fdr_ci_lower,
        global_fdr_ci_upper,
        avg_gene_fdr,
        gene_fdr_ci_lower,
        gene_fdr_ci_upper,
    ): (
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
        Vec<_>,
    ) = results.into_iter().multiunzip();

    list![
        gene = gene_indices,
        pathway = pathway_indices,
        similarity = similarities,
        avg_pval = avg_pval,
        pval_ci_lower = pval_ci_lower,
        pval_ci_upper = pval_ci_upper,
        avg_global_fdr = avg_global_fdr,
        global_fdr_ci_lower = global_fdr_ci_lower,
        global_fdr_ci_upper = global_fdr_ci_upper,
        avg_gene_fdr = avg_gene_fdr,
        gene_fdr_ci_lower = gene_fdr_ci_lower,
        gene_fdr_ci_upper = gene_fdr_ci_upper
    ]
}

///////////
// Utils //
///////////

/// Cosine similarity between two vectors
///
/// @description
/// Rust function to quickly calculate Cosine distances.
///
/// @param a Numeric vector. Vector a for which to calculate the similarity.
/// @param b Numeric vector. Vector b for which to calculate the similarity.
///
/// @returns Cosine similarity between the two vectors
///
/// @export
#[extendr]
fn rs_cosine_sim(a: &[f64], b: &[f64]) -> f64 {
    cosine_similarity(a, b)
}

/// Generate synthetic data for node2vec
///
/// @param test_data String. One of
/// `c("barbell", "caveman", "stochastic_block")`. Weird strings will default to
/// "barbell" data.
/// @param n_nodes_per_cluster Integer. Number of nodes in the test graph.
/// @param n_clusters Integer. Number of nodes per cluster.
/// @param p_within Numeric. Probability of edges within cluster (0-1).
/// @param p_between Numeric. Probability of edges between clusters (0-1).
/// @param seed Integer. Random seed.
///
/// @returns A list with the following elements:
/// \itemize{
///  \item edges - List with the edge data.
///  \item nodes - List with the node information.
/// }
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
        Node2VecDataType::Cavemen => node2vec_caveman(n_nodes_per_cluster, n_clusters, seed),
        Node2VecDataType::Stochastic => {
            node2vec_stochastic_block(n_nodes_per_cluster, n_clusters, p_within, p_between, seed)
        }
    };

    let (edge_data, node_data) = synthetic_data.generate_lists();

    list!(edges = edge_data, nodes = node_data)
}
