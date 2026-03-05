//! Utility functions such as cosine similarities, etc.

#![warn(missing_docs)]

use ann_search_rs::prelude::*;
use extendr_api::RMatrix;
use faer::{Mat, MatRef};
use num_traits::Float;
use rayon::prelude::*;

/////////////////////////
// Similarity measures //
/////////////////////////

/// Fast cosine similarity using auto-vectorization hints
///
/// Uses SIMD/auto-vectorisations under the hood and unsafe access of vectors
/// for speed.
///
/// ### Params
///
/// * `a` - Embedding vector a
/// * `b` - Embedding vector b
///
/// ### Returns
///
/// Cosine similarity between the two
#[inline(always)]
pub fn cosine_similarity<T: Float + SimdDistance>(a: &[T], b: &[T]) -> T {
    debug_assert_eq!(a.len(), b.len());

    let dot = T::dot_simd(a, b);
    let norm_a = T::calculate_l2_norm(a);
    let norm_b = T::calculate_l2_norm(b);

    dot / (norm_a * norm_b)
}

/// Cosine distances across gene and pathway embeddings
///
/// ### Params
///
/// * `genes` - Slice of gene embeddings
/// * `pathways` - Sloce of pathway embeddings
///
/// ### Returns
///
/// A faer matrix of similarities between the two of n_genes x n_pathways
pub fn compute_cross_cosine(genes: &[Vec<f32>], pathways: &[Vec<f32>]) -> Mat<f64> {
    let n_genes = genes.len();
    let n_pathways = pathways.len();

    // Compute similarity matrix with parallelization
    let similarities: Vec<f64> = (0..n_genes)
        .into_par_iter()
        .flat_map(|i| {
            (0..n_pathways)
                .map(|j| cosine_similarity(&genes[i], &pathways[j]) as f64)
                .collect::<Vec<f64>>()
        })
        .collect();

    Mat::from_fn(n_genes, n_pathways, |r, c| similarities[r * n_pathways + c])
}

///////////////////////
// P values and FDRs //
///////////////////////

/// Cosine distances across gene and pathway embeddings
///
/// ### Params
///
/// * `genes` - Slice of gene embeddings
/// * `pathways` - Sloce of pathway embeddings
///
/// ### Returns
///
/// A faer matrix of p-values between the two of n_genes x n_pathways based
/// on the null distribution
pub fn calculate_p_vals(
    observed_similarities: &MatRef<f64>,
    null_distribution: Vec<f64>,
) -> Mat<f64> {
    let n_genes = observed_similarities.nrows();
    let n_pathways = observed_similarities.ncols();
    let n_null = null_distribution.len() as f64;

    let mut sorted_null = null_distribution.clone();
    sorted_null.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let pvals: Vec<f64> = (0..n_genes)
        .into_par_iter()
        .flat_map(|i| {
            (0..n_pathways)
                .map(|j| {
                    let obs = observed_similarities[(i, j)];
                    let rank = sorted_null.partition_point(|&x| x < obs);
                    (n_null - rank as f64) / n_null
                })
                .collect::<Vec<f64>>()
        })
        .collect();

    Mat::from_fn(n_genes, n_pathways, |r, c| pvals[r * n_pathways + c])
}

///////////////
// R <> faer //
///////////////

/// Transform an R matrix to Vecs
///
/// ### Params
///
/// * `mat` - R matrix
///
/// ### Returns
///
/// Transformed into vectors
pub fn r_matrix_to_vec(mat: RMatrix<f64>) -> Vec<Vec<f32>> {
    let nrow = mat.nrows();
    let ncol = mat.ncols();

    (0..nrow)
        .map(|i| (0..ncol).map(|j| mat[[i, j]] as f32).collect())
        .collect()
}

////////////////
// Statistics //
////////////////

/// Geometric mean and CI matching GeneWalk's log_stats
///
/// ### Params
///
/// * `vals` - The p-values
///
/// ### Returns
///
/// Tuple of (geometric mean, lower CI, upper CI)
pub fn log_stats(vals: &[f64]) -> (f64, f64, f64) {
    let eps = 1e-16_f64;
    let n = vals.len() as f64;
    let log_vals: Vec<f64> = vals.iter().map(|&v| (v + eps).ln()).collect();
    let mean_log = log_vals.iter().sum::<f64>() / n;
    let g_mean = (mean_log.exp() - eps).clamp(0.0, 1.0);

    if vals.len() <= 1 {
        return (g_mean, eps, eps);
    }

    let variance = log_vals
        .iter()
        .map(|&x| (x - mean_log).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    let g_std = variance.sqrt().exp();
    let exponent = 1.96 / n.sqrt();

    (
        g_mean,
        (g_mean * g_std.powf(-exponent)).clamp(0.0, 1.0),
        (g_mean * g_std.powf(exponent)).clamp(0.0, 1.0),
    )
}
