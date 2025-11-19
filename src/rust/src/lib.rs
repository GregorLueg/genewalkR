pub mod embedding;
pub mod graph;

use extendr_api::prelude::*;
use node2vec_rs::train::TrainingConfig;
use std::ops::Deref;

////////////
// Params //
////////////

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod genewalkR;
}
