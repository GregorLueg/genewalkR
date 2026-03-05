//! Generating a node2vec embedding module

#![warn(missing_docs)]

use node2vec_rs::prelude::*;

////////////
// Config //
////////////

/// Configuration for GeneWalk training
///
/// Holds both random walk parameters and CPU training parameters.
pub struct GeneWalkConfig {
    /// Number of random walks per node
    pub walks_per_node: usize,
    /// Length of each random walk
    pub walk_length: usize,
    /// Return parameter for node2vec
    pub p: f32,
    /// In-out parameter for node2vec
    pub q: f32,
    /// Random seed for reproducibility
    pub seed: usize,
    /// CPU training arguments
    pub train_args: CpuTrainArgs,
}

impl GeneWalkConfig {
    /// Create a GeneWalkConfig from an R parameter list
    ///
    /// ### Params
    ///
    /// * `param_list` - Parameter list from R
    /// * `seed` - Random seed
    ///
    /// ### Returns
    ///
    /// Initialised GeneWalkConfig
    pub fn from_r_list(param_list: extendr_api::List, seed: usize) -> GeneWalkConfig {
        let param_list = param_list.into_hashmap();

        let p = param_list.get("p").and_then(|v| v.as_real()).unwrap_or(1.0) as f32;
        let q = param_list.get("q").and_then(|v| v.as_real()).unwrap_or(1.0) as f32;

        let walks_per_node = param_list
            .get("walks_per_node")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let walk_length = param_list
            .get("walk_length")
            .and_then(|v| v.as_integer())
            .unwrap_or(80) as usize;

        let num_workers = param_list
            .get("num_workers")
            .and_then(|v| v.as_integer())
            .unwrap_or(4) as usize;

        let n_epochs = param_list
            .get("n_epochs")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;

        let num_negatives = param_list
            .get("num_negatives")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;

        let window_size = param_list
            .get("window_size")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let lr = param_list
            .get("lr")
            .and_then(|v| v.as_real())
            .unwrap_or(0.025) as f32;

        let dim = param_list
            .get("dim")
            .and_then(|v| v.as_integer())
            .unwrap_or(128) as usize;

        let verbose = param_list
            .get("verbose")
            .and_then(|v| v.as_logical())
            .map(|v| v.is_true())
            .unwrap_or(false);

        let sample = param_list
            .get("sample")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-3) as f32;

        GeneWalkConfig {
            walks_per_node,
            walk_length,
            p,
            q,
            seed,
            train_args: CpuTrainArgs {
                dim,
                lr,
                epochs: n_epochs,
                neg: num_negatives,
                window: window_size,
                lr_update_rate: 10_000,
                n_threads: num_workers,
                verbose,
                sample,
            },
        }
    }
}

/// Train node2vec on the CPU and return combined embeddings
///
/// ### Params
///
/// * `walks` - The generated random walks
/// * `vocab_size` - Number of unique nodes
/// * `config` - The GeneWalk configuration
/// * `verbose` - Controls verbosity
///
/// ### Returns
///
/// Combined embeddings as Vec<Vec<f32>>, indexed by node ID
pub fn train_node2vec(
    walks: Vec<Vec<u32>>,
    vocab_size: usize,
    config: &GeneWalkConfig,
    verbose: bool,
) -> Vec<Vec<f32>> {
    let neg_table = create_negative_table(vocab_size, &walks, NEGATIVE_TABLE_SIZE, config.seed);

    let mut args = config.train_args.clone();
    args.verbose = verbose;

    let (mut input_mat, _) = train_node2vec_cpu(walks, vocab_size, args, neg_table, config.seed);

    input_mat.norm_self();

    (0..input_mat.n_rows())
        .map(|i| input_mat.row_as_slice(i).to_vec())
        .collect()
}
