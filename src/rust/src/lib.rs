pub mod embedding;
pub mod graph;
pub mod utils;

use burn::backend::{
    libtorch::{LibTorch, LibTorchDevice},
    Autodiff,
};
use burn::prelude::Backend;
use extendr_api::prelude::*;
use node2vec_rs::model::SkipGramConfig;
use std::time::Instant;

use crate::embedding::*;
use crate::graph::*;

// export the function to R
extendr_module! {
    mod genewalkR;
    fn rs_gene_walk;
    fn rs_gene_walk_perm;
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
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @return A numerical matrix of number nodes x desired embedding dimensios.
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
        println!("Running the random walks.")
    }

    let start_rw = Instant::now();

    let node2vec_graph = create_graph(&edges, gene_walk_config.p, gene_walk_config.q, directed);

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

    let device = LibTorchDevice::Cpu;
    LibTorch::<f32>::seed(&device, gene_walk_config.seed);

    let skipgram_config = SkipGramConfig {
        vocab_size: random_walks.len(),
        embedding_dim: embd_dim,
    };

    let embedding: Vec<Vec<f32>> = train_node2vec::<Autodiff<LibTorch>>(
        skipgram_config,
        gene_walk_config,
        random_walks,
        device,
        verbose,
    );

    let nrows = embedding.len();
    let ncols = embedding[0].len();

    RMatrix::new_matrix(nrows, ncols, |r, c| embedding[r][c] as f64)
}

/// Function to generate permuted embeddings
///
/// @description
/// Wrapper function to generate permuted embeddings.
///
/// @param from Integer vector. The node indices indicating where the edge
/// originates from.
/// @param to Integer vector. The node indices indicating where the edge
/// goes to.
/// @param weights Optional numeric vector. If not supplied, defaults to `1.0`
/// as edge weight.
/// @param gene_walk_params Named list. Contains the parameters for running
/// gene walk.
/// @param n_perm Integer.
/// @param embd_dim Integer. Embedding dimension.
/// @param directed Boolean. Is the graph directed. If set to `FALSE` reverse
/// edges will be added.
/// @param seed Integer. For reproducibility of the algorithm.
/// @param verbose Boolean. Controls verbosity of the algorithm.
///
/// @returns A list of matrices of permuted embeddings.
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
    let gene_walk_config = GeneWalkConfig::from_r_list(gene_walk_params, seed);

    let device = LibTorchDevice::Cpu;
    LibTorch::<f32>::seed(&device, gene_walk_config.seed);

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

    let mut perm_res: Vec<Vec<Vec<f32>>> = Vec::with_capacity(n_perm);

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
            vocab_size: random_walks_perm.len(),
            embedding_dim: embd_dim,
        };

        let embd_i: Vec<Vec<f32>> = train_node2vec::<Autodiff<LibTorch>>(
            skipgram_config,
            gene_walk_config_i,
            random_walks_perm,
            device,
            verbose,
        );

        perm_res.push(embd_i);

        let end_perm = start_perm.elapsed();

        if verbose {
            println!(
                " Finished the generation of the embedding on permutated data in {:.2}s",
                end_perm.as_secs_f64()
            );
        }
    }
    let nrows = perm_res[0].len();
    let ncols = perm_res[0][0].len();

    let mut res_ls = List::new(n_perm);

    #[allow(clippy::needless_range_loop)]
    for i in 0..n_perm {
        let r_mat_i = RMatrix::new_matrix(nrows, ncols, |r, c| perm_res[i][r][c] as f64);
        res_ls.set_elt(i, Robj::from(r_mat_i))?;
    }

    Ok(res_ls)
}
