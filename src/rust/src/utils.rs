use extendr_api::RMatrix;
use faer::{Mat, MatRef};
use rayon::prelude::*;

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
pub fn cosine_similarity(a: &[f32], b: &[f32]) -> f32 {
    debug_assert_eq!(a.len(), b.len());

    let mut dot = 0.0_f32;
    let mut norm_a = 0.0_f32;
    let mut norm_b = 0.0_f32;

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

/// Helper function to transform faer to R
///
/// ### Params
///
/// * `mat` - The faer matrix
///
/// ### Returns
///
/// The R-matrix
pub fn faer_to_r_matrix(mat: &MatRef<f64>) -> RMatrix<f64> {
    let nrow = mat.nrows();
    let ncol = mat.ncols();
    RMatrix::new_matrix(nrow, ncol, |row, column| mat[(row, column)])
}

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
