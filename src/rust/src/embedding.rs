use burn::data::dataloader::DataLoaderBuilder;
use burn::optim::{AdamConfig, GradientsParams, Optimizer};
use burn::prelude::ElementConversion;
use burn::tensor::backend::AutodiffBackend;
use extendr_api::List;
use node2vec_rs::prelude::*;
use std::ops::Deref;
use std::time::Instant;

////////////
// Params //
////////////

#[derive(Clone, Debug)]
pub struct GeneWalkConfig(TrainingConfig);

impl GeneWalkConfig {
    /// Helper function to generate the Trainings configuration
    ///
    /// ### Params
    ///
    /// * `param_list` - Parameter list for training from R
    /// * `seed` - Random seed
    ///
    /// ### Returns
    ///
    /// `GeneWalkConfig` that stores the TrainingConfig internally
    pub fn from_r_list(param_list: List, seed: usize) -> GeneWalkConfig {
        let param_list = param_list.into_hashmap();

        // random walk params
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

        // training params
        let num_workers = param_list
            .get("num_workers")
            .and_then(|v| v.as_integer())
            .unwrap_or(4) as usize;

        let batch_size = param_list
            .get("batch_size")
            .and_then(|v| v.as_integer())
            .unwrap_or(256) as usize;

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
            .unwrap_or(2) as usize;

        let lr = param_list
            .get("lr")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-3);

        Self(TrainingConfig {
            walks_per_node,
            walk_length,
            window_size,
            batch_size,
            num_workers,
            num_epochs: n_epochs,
            num_negatives,
            p,
            q,
            learning_rate: lr,
            seed: seed as u64,
        })
    }

    /// Get a reference to the TrainingConfig
    ///
    /// ### Returns
    ///
    /// Reference to the TrainingConfig
    pub fn inner(&self) -> &TrainingConfig {
        &self.0
    }

    /// Get the TrainingConfig
    ///
    /// ### Returns
    ///
    /// Returns the TrainingConfig
    pub fn into_innter(self) -> TrainingConfig {
        self.0
    }
}

/// Derefencing trait for the GeneWalkConfig
impl Deref for GeneWalkConfig {
    type Target = TrainingConfig;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

/// Train and return the `node2vec` model
///
/// ### Params
///
/// * `model_config` - The model configuration. Contains the embedding
///   dimensions and vocabulary size.
/// * `genewalk_config` - The GeneWalk configuration that contains all of the
///   parameters.
/// * `random_walks` - The generated random walks to train the model.
/// * `device` - The device on which to run the training.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// The embedding of the model
pub fn train_node2vec<B: AutodiffBackend>(
    model_config: SkipGramConfig,
    genewalk_config: GeneWalkConfig,
    random_walks: Vec<Vec<u32>>,
    device: B::Device,
    verbose: bool,
) -> Vec<Vec<f32>> {
    let mut model = model_config.init::<B>(&device);
    let mut optim = AdamConfig::new().init();
    let batcher = SkipGramBatcher::new(genewalk_config.window_size);

    let data_loader = DataLoaderBuilder::new(batcher.clone())
        .batch_size(genewalk_config.batch_size)
        .shuffle(genewalk_config.seed)
        .num_workers(genewalk_config.num_workers)
        .build(WalkDataset::new(random_walks));

    let total_batches = data_loader.num_items() / genewalk_config.batch_size;

    for epoch in 1..=genewalk_config.num_epochs {
        let epoch_start = Instant::now();
        let mut total_loss = 0.0;
        let mut num_batches = 0;

        for batch in data_loader.iter() {
            let batch_size = batch.contexts.dims()[0];
            let negatives = sample_negatives(
                batch_size,
                model.vocab_size,
                genewalk_config.num_negatives,
                &batch.contexts.device(),
            );

            let loss = model
                .forward(batch.centers, batch.contexts, negatives)
                .mean();

            let loss_scalar: f64 = loss.clone().into_scalar().elem();
            total_loss += loss_scalar;
            num_batches += 1;

            if verbose && num_batches % 100 == 0 {
                println!(
                    "   Processed {} batches out of {} for epoch {} out of {}.",
                    num_batches, total_batches, epoch, genewalk_config.num_epochs
                );
            }

            let grads = loss.backward();
            let grads = GradientsParams::from_grads(grads, &model);
            model = optim.step(genewalk_config.learning_rate, model, grads);
        }

        if verbose {
            let epoch_duration = epoch_start.elapsed();
            let avg_loss = total_loss / num_batches as f64;
            println!(
                "  Epoch {} completed in {:.2}s - Total loss: {:.6}, Avg loss: {:.6}",
                epoch,
                epoch_duration.as_secs_f64(),
                total_loss,
                avg_loss
            );
        }
    }

    model.combined_embeddings_to_vec()
}
