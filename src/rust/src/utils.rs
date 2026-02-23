use extendr_api::RMatrix;
use faer::{Mat, MatRef};
use num_traits::Float;
use rayon::prelude::*;
use std::ops::AddAssign;

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
pub fn cosine_similarity<T: Float + AddAssign>(a: &[T], b: &[T]) -> T {
    debug_assert_eq!(a.len(), b.len());

    let mut dot = T::zero();
    let mut norm_a = T::zero();
    let mut norm_b = T::zero();

    let len = a.len();
    let chunks = len / 4;
    let remainder = len % 4;

    // usual 4 at once
    for i in 0..chunks {
        let idx = i * 4;
        unsafe {
            let a0 = *a.get_unchecked(idx);
            let a1 = *a.get_unchecked(idx + 1);
            let a2 = *a.get_unchecked(idx + 2);
            let a3 = *a.get_unchecked(idx + 3);

            let b0 = *b.get_unchecked(idx);
            let b1 = *b.get_unchecked(idx + 1);
            let b2 = *b.get_unchecked(idx + 2);
            let b3 = *b.get_unchecked(idx + 3);

            dot += a0 * b0 + a1 * b1 + a2 * b2 + a3 * b3;
            norm_a += a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
            norm_b += b0 * b0 + b1 * b1 + b2 * b2 + b3 * b3;
        }
    }

    // remainders
    for i in (len - remainder)..len {
        unsafe {
            let ai = *a.get_unchecked(i);
            let bi = *b.get_unchecked(i);
            dot += ai * bi;
            norm_a += ai * ai;
            norm_b += bi * bi;
        }
    }

    dot / (norm_a.sqrt() * norm_b.sqrt())
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

/// Calculate the FDR
///
/// ### Params
///
/// * `pvals` - P-values for which to calculate the FDR
///
/// ### Returns
///
/// The calculated FDRs
pub fn calc_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    if n == 0 {
        return vec![];
    }

    let n_f64 = n as f64;
    let mut indexed_pval: Vec<(usize, f64)> =
        pvals.iter().enumerate().map(|(i, &x)| (i, x)).collect();

    indexed_pval
        .sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut adj_pvals = vec![0.0; n];
    let mut current_min = ((n_f64 / n as f64) * indexed_pval[n - 1].1).min(1.0);
    adj_pvals[indexed_pval[n - 1].0] = current_min;

    for i in (0..n - 1).rev() {
        let adj_val = ((n_f64 / (i + 1) as f64) * indexed_pval[i].1)
            .min(current_min)
            .min(1.0);
        current_min = adj_val;
        adj_pvals[indexed_pval[i].0] = adj_val;
    }

    adj_pvals
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

///////////////
// Quantiles //
///////////////

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
