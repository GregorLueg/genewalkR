use crate::graph::*;

use faer::linalg::matmul::matmul;
use faer::matrix_free::eigen::{
    partial_eigen_scratch, partial_self_adjoint_eigen, PartialEigenParams,
};
use faer::matrix_free::LinOp;
use faer::{
    col::generic::Col,
    dyn_stack::{MemBuffer, MemStack, StackReq},
};
use faer::{Accum, Mat, MatMut, MatRef, Par, Side};
use rand::{rngs::StdRng, seq::index, RngCore, SeedableRng};
use rayon::prelude::*;
use std::f64::consts::PI;
use std::num::NonZeroUsize;

/////////////
// Helpers //
/////////////

/// Get the default number of cores for `Par`
///
/// ### Returns
///
/// `Par` object for faer with default threads
fn default_par() -> Par {
    Par::Rayon(NonZeroUsize::new(rayon::current_num_threads()).unwrap_or(NonZeroUsize::MIN))
}

/// Extract rows of U at the given indices
///
/// ### Params
///
/// * `u` - Eigenvector matrix, n x k
/// * `bkgd` - Row indices to extract
///
/// ### Returns
///
/// Submatrix of shape (|bkgd|, k)
fn extract_u_bkgd(u: MatRef<f64>, bkgd: &[usize]) -> Mat<f64> {
    let k = u.ncols();
    let mut u_b = Mat::zeros(bkgd.len(), k);
    for (i, &node) in bkgd.iter().enumerate() {
        for j in 0..k {
            u_b[(i, j)] = u[(node, j)];
        }
    }
    u_b
}

/// Materialise the full dense kernel K = U * diag(f) * U^T
///
/// Useful for testing and for small graphs where you actually
/// want the dense matrix.
///
/// ### Params
///
/// * `decomp` - Spectral decomposition
/// * `f_lambda` - Transformed eigenvalues from `kernel_eigenvalues`
///
/// ### Returns
///
/// Dense n x n kernel matrix
pub fn materialise_kernel(decomp: &SpectralDecomp, f_lambda: &[f64]) -> Mat<f64> {
    let n = decomp.n;
    let k = decomp.eigenvalues.len();
    let u = decomp.eigenvectors.as_ref();

    // UF = U * diag(f)
    let mut uf = Mat::zeros(n, k);
    for j in 0..k {
        for i in 0..n {
            uf[(i, j)] = u[(i, j)] * f_lambda[j];
        }
    }

    // K = UF * U^T
    let mut kernel = Mat::zeros(n, n);
    matmul(
        kernel.as_mut(),
        Accum::Replace,
        uf.as_ref(),
        u.transpose(),
        1.0,
        default_par(),
    );
    kernel
}

///////////////
// Laplacian //
///////////////

/// Build the graph Laplacian from a symmetric adjacency matrix
///
/// For the unnormalised variant: L = D - A.
/// For the normalised variant: L_norm = I - D^{-1/2} A D^{-1/2}.
///
/// Self-loops are ignored. Isolated nodes get a zero diagonal in the
/// normalised case (matching igraph's convention).
///
/// ### Params
///
/// * `adj` - Symmetric adjacency matrix in CSR format. Column indices
///   within each row must be sorted ascending.
/// * `normalised` - Whether to return the normalised Laplacian
///
/// ### Returns
///
/// The Laplacian in CSR format, same shape as `adj`
pub fn laplacian(adj: &CompressedSparseData, normalised: bool) -> CompressedSparseData {
    assert!(adj.cs_type.is_csr(), "adjacency must be CSR");
    let n = adj.nrows();
    assert_eq!(adj.ncols(), n, "adjacency must be square");

    // weighted degrees (ignoring self-loops)
    let mut deg = vec![0.0f64; n];
    for i in 0..n {
        for idx in adj.indptr[i]..adj.indptr[i + 1] {
            if adj.indices[idx] != i {
                deg[i] += adj.data[idx];
            }
        }
    }

    let inv_sqrt_deg: Vec<f64> = if normalised {
        deg.iter()
            .map(|&d| if d > 0.0 { 1.0 / d.sqrt() } else { 0.0 })
            .collect()
    } else {
        vec![]
    };

    let mut data = Vec::with_capacity(adj.get_nnz() + n);
    let mut indices = Vec::with_capacity(adj.get_nnz() + n);
    let mut indptr = vec![0usize];

    for i in 0..n {
        let mut diag_inserted = false;

        for idx in adj.indptr[i]..adj.indptr[i + 1] {
            let j = adj.indices[idx];
            if j == i {
                continue;
            }

            // insert diagonal before first column index > i
            if !diag_inserted && j > i {
                indices.push(i);
                data.push(if normalised {
                    if deg[i] > 0.0 {
                        1.0
                    } else {
                        0.0
                    }
                } else {
                    deg[i]
                });
                diag_inserted = true;
            }

            indices.push(j);
            data.push(if normalised {
                -adj.data[idx] * inv_sqrt_deg[i] * inv_sqrt_deg[j]
            } else {
                -adj.data[idx]
            });
        }

        // diagonal at end if all neighbours had index < i
        if !diag_inserted {
            indices.push(i);
            data.push(if normalised {
                if deg[i] > 0.0 {
                    1.0
                } else {
                    0.0
                }
            } else {
                deg[i]
            });
        }

        indptr.push(data.len());
    }

    CompressedSparseData::new_csr(&data, &indices, &indptr, (n, n))
}

////////////////////////////
// Spectral decomposition //
////////////////////////////

/////////////
// Helpers //
/////////////

/// Wraps a sparse Laplacian as (shift * I - L) so that the smallest
/// eigenvalues of L become the largest eigenvalues of this operator.
///
/// Used with `partial_self_adjoint_eigen` for truncated spectral
/// decomposition without ever materialising the dense Laplacian.
#[derive(Debug)]
pub struct ShiftedLaplacian<'a> {
    /// Sparse Laplacian in CSR format
    pub lap: &'a CompressedSparseData,
    /// Spectral shift (>= lambda_max of L)
    pub shift: f64,
}

impl<'a> LinOp<f64> for ShiftedLaplacian<'a> {
    fn apply_scratch(&self, _rhs_ncols: usize, _par: Par) -> StackReq {
        StackReq::EMPTY
    }

    fn nrows(&self) -> usize {
        self.lap.nrows()
    }

    fn ncols(&self) -> usize {
        self.lap.ncols()
    }

    /// Computes out = (shift * I - L) * rhs
    fn apply(
        &self,
        mut out: MatMut<'_, f64>,
        rhs: MatRef<'_, f64>,
        _par: Par,
        _stack: &mut MemStack,
    ) {
        let n = self.lap.nrows();
        let ncols = rhs.ncols();

        for i in 0..n {
            for col in 0..ncols {
                out[(i, col)] = self.shift * rhs[(i, col)];
            }
            for idx in self.lap.indptr[i]..self.lap.indptr[i + 1] {
                let j = self.lap.indices[idx];
                let v = self.lap.data[idx];
                for col in 0..ncols {
                    out[(i, col)] -= v * rhs[(j, col)];
                }
            }
        }
    }

    /// For real f64, conjugate application is identical to apply
    fn conj_apply(
        &self,
        out: MatMut<'_, f64>,
        rhs: MatRef<'_, f64>,
        par: Par,
        stack: &mut MemStack,
    ) {
        self.apply(out, rhs, par, stack);
    }
}

//////////
// Main //
//////////

/// Strategy for the spectral decomposition
pub enum SpectralStrategy {
    /// Full O(n^3) decomposition via dense symmetric EVD.
    /// Exact. Needs O(n^2) memory for the dense Laplacian.
    Full,
    /// Truncated decomposition keeping the k smallest eigenvalues.
    /// Uses spectral flip + Lanczos (never builds the dense Laplacian).
    Truncated {
        /// Number of eigenvalues to keep
        k: usize,
        /// Convergence tolerance (defaults to 128 * f64::EPSILON)
        tolerance: Option<f64>,
    },
}

/// Spectral decomposition of a graph Laplacian
///
/// Stores eigenvalues and eigenvectors so that any kernel can be
/// applied without recomputing the decomposition and without ever
/// materialising the full dense kernel.
#[derive(Debug, Clone)]
pub struct SpectralDecomp {
    /// Eigenvalues of L in ascending order, length k <= n
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors as columns, n x k
    pub eigenvectors: Mat<f64>,
    /// Total nodes in the graph
    pub n: usize,
}

/// Eigendecompose the graph Laplacian
///
/// For `SpectralStrategy::Full`, converts the sparse Laplacian to dense and
/// calls faer's `self_adjoint_eigen`. For `SpectralStrategy::Truncated`,
/// uses the spectral flip trick with iterative Lanczos to find the k
/// smallest eigenvalues without materialising the dense matrix.
///
/// ### Params
///
/// * `lap` - Graph Laplacian in CSR format (must be symmetric)
/// * `normalised` - Whether the Laplacian is normalised (affects the
///   spectral flip shift value)
/// * `strategy` - Full or truncated decomposition
///
/// ### Returns
///
/// `SpectralDecomp` with eigenvalues in ascending order
pub fn spectral_decompose(
    lap: &CompressedSparseData,
    normalised: bool,
    strategy: SpectralStrategy,
) -> SpectralDecomp {
    let n = lap.nrows();

    match strategy {
        SpectralStrategy::Full => {
            let dense = lap.to_dense();
            let evd = dense
                .self_adjoint_eigen(Side::Lower)
                .expect("symmetric EVD failed");

            let eigenvalues: Vec<f64> = evd.S().column_vector().iter().copied().collect();

            let eigenvectors = evd.U().to_owned();

            SpectralDecomp {
                eigenvalues,
                eigenvectors,
                n,
            }
        }

        SpectralStrategy::Truncated { k, tolerance } => {
            assert!(k > 0 && k <= n, "k must be in [1, n]");
            let tol = tolerance.unwrap_or(128.0 * f64::EPSILON);

            // Shift: for normalised L, eigenvalues in [0, 2].
            // For unnormalised, use Gershgorin: lambda_max <= 2 * max_degree.
            let shift = if normalised {
                2.0
            } else {
                let mut max_diag = 0.0f64;
                for i in 0..n {
                    for idx in lap.indptr[i]..lap.indptr[i + 1] {
                        if lap.indices[idx] == i {
                            max_diag = max_diag.max(lap.data[idx]);
                        }
                    }
                }
                2.0 * max_diag
            };

            let shifted = ShiftedLaplacian { lap, shift };

            // initial vector: uniform, normalised
            let v0_val = 1.0 / (n as f64).sqrt();
            let v0 = Col::from_fn(n, |_| v0_val);

            let mut eigvecs = Mat::<f64>::zeros(n, k);
            let mut eigvals = vec![0.0f64; k];

            let params = PartialEigenParams {
                max_restarts: 1000,
                ..Default::default()
            };

            let par = default_par();

            let scratch = partial_eigen_scratch(&shifted, k, par, params);
            let mut mem = MemBuffer::new(scratch);

            let info = partial_self_adjoint_eigen(
                eigvecs.as_mut(),
                &mut eigvals,
                &shifted,
                v0.as_ref(),
                tol,
                par,
                MemStack::new(&mut mem),
                params,
            );

            let n_converged = info.n_converged_eigen;
            if n_converged < k {
                eprintln!(
                    "warning: only {}/{} eigenvalues converged after {} restarts",
                    n_converged, k, params.max_restarts
                );
            }

            // partial_self_adjoint_eigen returns eigenvalues of
            // (shift * I - L) sorted by decreasing magnitude.
            // Recover original: lambda_L = shift - lambda_shifted.
            // Since the shifted eigenvalues are in decreasing order,
            // the recovered values come out in ascending order -- which
            // is exactly what we want.
            for i in 0..k {
                eigvals[i] = shift - eigvals[i];
            }

            // clamp tiny negatives from numerical noise to zero
            for v in eigvals.iter_mut() {
                if *v < 0.0 && *v > -tol {
                    *v = 0.0;
                }
            }

            SpectralDecomp {
                eigenvalues: eigvals,
                eigenvectors: eigvecs,
                n,
            }
        }
    }
}

///////////////////////
// Kernel transforms //
///////////////////////

/// Kernel type with its parameters
pub enum KernelType {
    /// (sigma2 * L + add_diag * I)^{-1}
    RegularisedLaplacian { sigma2: f64, add_diag: f64 },
    /// Moore-Penrose pseudoinverse of L
    CommuteTime,
    /// exp(-sigma2/2 * L)
    Diffusion { sigma2: f64 },
    /// U * diag(cos(pi/4 * lambda)) * U^T -- normalised Laplacian only
    InverseCosine,
    /// (a*I - L)^p -- normalised Laplacian only, requires a >= 2
    PStep { a: f64, p: u32 },
}

/// Apply a kernel's scalar transform to the Laplacian eigenvalues
///
/// Each kernel type maps eigenvalue lambda to a scalar f(lambda).
/// The kernel matrix is then K = U * diag(f(lambda)) * U^T but we
/// never materialise it.
///
/// ### Params
///
/// * `eigenvalues` - Laplacian eigenvalues (ascending)
/// * `kernel` - Which kernel and its parameters
///
/// ### Returns
///
/// Vector of transformed eigenvalues, same length as input
pub fn kernel_eigenvalues(eigenvalues: &[f64], kernel: &KernelType) -> Vec<f64> {
    let eps = f64::EPSILON;

    eigenvalues
        .iter()
        .map(|&lam| match kernel {
            KernelType::RegularisedLaplacian { sigma2, add_diag } => {
                1.0 / (sigma2 * lam + add_diag)
            }
            KernelType::CommuteTime => {
                if lam.abs() > eps {
                    1.0 / lam
                } else {
                    0.0
                }
            }
            KernelType::Diffusion { sigma2 } => (-sigma2 / 2.0 * lam).exp(),
            KernelType::InverseCosine => (PI / 4.0 * lam).cos(),
            KernelType::PStep { a, p } => (a - lam).powi(*p as i32),
        })
        .collect()
}

///////////////////
// Scoring (raw) //
///////////////////

/// Compute raw diffusion scores: score = K[:, bkgd] * input
///
/// Uses the factored form U * diag(f) * U_bkgd^T * input so that
/// the dense kernel is never materialised.
///
/// ### Params
///
/// * `decomp` - Spectral decomposition of the Laplacian
/// * `f_lambda` - Transformed eigenvalues from `kernel_eigenvalues`
/// * `input` - Dense score matrix, shape (n_bkgd, n_inputs)
/// * `bkgd` - Indices mapping background rows to graph nodes
///
/// ### Returns
///
/// Dense matrix of shape (n, n_inputs) with scores for every node
pub fn diffuse_raw(
    decomp: &SpectralDecomp,
    f_lambda: &[f64],
    input: MatRef<f64>,
    bkgd: &[usize],
) -> Mat<f64> {
    let n = decomp.n;
    let k = decomp.eigenvalues.len();
    let n_inputs = input.ncols();
    let u = decomp.eigenvectors.as_ref();

    // U_b = U[bkgd, :], shape (n_bkgd, k)
    let u_b = extract_u_bkgd(u, bkgd);

    // T = U_b^T * input, shape (k, n_inputs)
    let mut t = Mat::zeros(k, n_inputs);
    matmul(
        t.as_mut(),
        Accum::Replace,
        u_b.transpose(),
        input,
        1.0,
        default_par(),
    );

    // T = diag(f) * T
    for j in 0..n_inputs {
        for i in 0..k {
            t[(i, j)] *= f_lambda[i];
        }
    }

    // score = U * T, shape (n, n_inputs)
    let mut score = Mat::zeros(n, n_inputs);
    matmul(
        score.as_mut(),
        Accum::Replace,
        u,
        t.as_ref(),
        1.0,
        default_par(),
    );

    score
}

/////////////////
// Scoring (Z) //
/////////////////

/// Compute z-score normalised diffusion scores
///
/// The z-score for each node is (raw - mean) / sqrt(var), where the
/// mean and variance are the analytical moments of the score under
/// a uniformly random permutation of the input labels. See
/// Picart-Armada et al., Bioinformatics 2020 for derivation.
///
/// All moments are computed from the spectral decomposition without
/// materialising the kernel.
///
/// ### Params
///
/// * `decomp` - Spectral decomposition of the Laplacian
/// * `f_lambda` - Transformed eigenvalues from `kernel_eigenvalues`
/// * `input` - Dense score matrix, shape (n_bkgd, n_inputs)
/// * `bkgd` - Indices mapping background rows to graph nodes
///
/// ### Returns
///
/// Dense matrix of shape (n, n_inputs) with z-scores for every node.
/// Nodes disconnected from all labelled nodes may get NaN.
pub fn diffuse_z(
    decomp: &SpectralDecomp,
    f_lambda: &[f64],
    input: MatRef<f64>,
    bkgd: &[usize],
) -> Mat<f64> {
    let n = decomp.n;
    let k = decomp.eigenvalues.len();
    let n_bkgd = bkgd.len();
    let n_inputs = input.ncols();
    let u = decomp.eigenvectors.as_ref();

    let u_b = extract_u_bkgd(u, bkgd);

    // --- row_sums of K[:, bkgd] ---
    // c[l] = sum_j U_b[j, l] (column sums of U_b)
    let mut c = vec![0.0f64; k];
    for j in 0..n_bkgd {
        for l in 0..k {
            c[l] += u_b[(j, l)];
        }
    }
    // fc[l] = f[l] * c[l]
    let fc: Vec<f64> = c.iter().zip(f_lambda).map(|(&ci, &fi)| fi * ci).collect();
    // row_sums[i] = sum_l U[i,l] * fc[l]
    let mut row_sums = vec![0.0f64; n];
    for i in 0..n {
        for l in 0..k {
            row_sums[i] += u[(i, l)] * fc[l];
        }
    }

    // --- row_sums2 of K[:, bkgd].^2 ---
    // M = diag(f) * (U_b^T * U_b) * diag(f), shape (k, k)
    let mut ub_t_ub = Mat::zeros(k, k);
    matmul(
        ub_t_ub.as_mut(),
        Accum::Replace,
        u_b.transpose(),
        u_b.as_ref(),
        1.0,
        default_par(),
    );
    // scale: M[l, m] = f[l] * ub_t_ub[l, m] * f[m]
    for l in 0..k {
        for m in 0..k {
            ub_t_ub[(l, m)] *= f_lambda[l] * f_lambda[m];
        }
    }
    let m_mat = ub_t_ub;

    // C = U * M, shape (n, k)
    let mut um = Mat::zeros(n, k);
    matmul(
        um.as_mut(),
        Accum::Replace,
        u,
        m_mat.as_ref(),
        1.0,
        default_par(),
    );
    // row_sums2[i] = sum_l C[i,l] * U[i,l]
    let mut row_sums2 = vec![0.0f64; n];
    for i in 0..n {
        for l in 0..k {
            row_sums2[i] += um[(i, l)] * u[(i, l)];
        }
    }

    // --- const_mean and const_var ---
    let nb = n_bkgd as f64;
    let mut const_mean = vec![0.0f64; n];
    let mut const_var = vec![0.0f64; n];
    for i in 0..n {
        const_mean[i] = row_sums[i] / nb;
        const_var[i] = (nb * row_sums2[i] - row_sums[i] * row_sums[i]) / ((nb - 1.0) * nb * nb);
    }

    // --- raw scores (reuse diffuse_raw logic) ---
    let raw = diffuse_raw(decomp, f_lambda, input, bkgd);

    // --- z-scores per input column ---
    let mut z = Mat::zeros(n, n_inputs);
    for col in 0..n_inputs {
        let mut s1 = 0.0f64;
        let mut s2 = 0.0f64;
        for j in 0..n_bkgd {
            let v = input[(j, col)];
            s1 += v;
            s2 += v * v;
        }

        for i in 0..n {
            let mean = const_mean[i] * s1;
            let var = const_var[i] * (nb * s2 - s1 * s1);
            z[(i, col)] = (raw[(i, col)] - mean) / var.sqrt();
        }
    }

    z
}

////////////////////////////////////
// Scoring (Monte Carlo heatrank) //
////////////////////////////////////

/// Compute empirical heatrank via permutation testing
///
/// For each input column, generates `n_perm` random permutations of
/// the input labels within the background, computes the permuted
/// diffusion score, and counts how many exceed the original score.
/// Returns the corrected heatrank (r + 1) / (n_perm + 1).
///
/// Parallelised over permutations using rayon.
///
/// ### Params
///
/// * `decomp` - Spectral decomposition of the Laplacian
/// * `f_lambda` - Transformed eigenvalues from `kernel_eigenvalues`
/// * `input` - Dense score matrix, shape (n_bkgd, n_inputs)
/// * `bkgd` - Indices mapping background rows to graph nodes
/// * `n_perm` - Number of permutations
/// * `seed` - RNG seed for reproducibility
///
/// ### Returns
///
/// Dense matrix of shape (n, n_inputs) with heatranks in (0, 1).
/// Smaller values indicate warmer (more significant) nodes.
pub fn diffuse_mc(
    decomp: &SpectralDecomp,
    f_lambda: &[f64],
    input: MatRef<f64>,
    bkgd: &[usize],
    n_perm: usize,
    seed: u64,
) -> Mat<f64> {
    let n = decomp.n;
    let k = decomp.eigenvalues.len();
    let n_bkgd = bkgd.len();
    let n_inputs = input.ncols();
    let u = decomp.eigenvectors.as_ref();

    // precompute U_bkgd (n_bkgd x k) and UF = U * diag(f) (n x k)
    let u_b = extract_u_bkgd(u, bkgd);
    let mut uf = Mat::zeros(n, k);
    for j in 0..k {
        for i in 0..n {
            uf[(i, j)] = u[(i, j)] * f_lambda[j];
        }
    }

    // original scores (n x n_inputs)
    let mut raw = Mat::zeros(n, n_inputs);
    for col in 0..n_inputs {
        // same t-vector construction as in the permutation loop
        let mut t = vec![0.0f64; k];
        for j in 0..n_bkgd {
            let v = input[(j, col)];
            if v != 0.0 {
                for l in 0..k {
                    t[l] += v * u_b[(j, l)];
                }
            }
        }
        // same score reconstruction as in the permutation loop
        for i in 0..n {
            let mut score = 0.0f64;
            for l in 0..k {
                score += uf[(i, l)] * t[l];
            }
            raw[(i, col)] = score;
        }
    }

    // per-permutation seeds for reproducibility across thread counts
    let mut master_rng = StdRng::seed_from_u64(seed);
    let perm_seeds: Vec<u64> = (0..n_perm).map(|_| master_rng.next_u64()).collect();

    let mut result = Mat::zeros(n, n_inputs);

    for col in 0..n_inputs {
        // extract nonzero positions and values for this input column
        let mut nz_positions = Vec::new();
        let mut nz_values = Vec::new();
        for j in 0..n_bkgd {
            let v = input[(j, col)];
            if v != 0.0 {
                nz_positions.push(j);
                nz_values.push(v);
            }
        }
        let nnz = nz_positions.len();
        if nnz == 0 {
            // no input signal: heatrank is uniform
            for i in 0..n {
                result[(i, col)] = 1.0;
            }
            continue;
        }

        let orig_col: Vec<f64> = (0..n).map(|i| raw[(i, col)]).collect();

        // accumulate exceedance counts in parallel
        let counts: Vec<u32> = perm_seeds
            .par_iter()
            .fold(
                || vec![0u32; n],
                |mut acc, &pseed| {
                    let mut rng = StdRng::seed_from_u64(pseed);

                    // sample nnz indices from 0..n_bkgd without replacement
                    let sampled = index::sample(&mut rng, n_bkgd, nnz);

                    // compute permuted score via factored form:
                    // t = sum_j nz_values[j] * U_b[perm[j], :]  (k-vector)
                    let mut t = vec![0.0f64; k];
                    for (j, perm_idx) in sampled.iter().enumerate() {
                        let val = nz_values[j];
                        for l in 0..k {
                            t[l] += val * u_b[(perm_idx, l)];
                        }
                    }

                    // perm_score = UF * t  (n-vector)
                    // count exceedances
                    for i in 0..n {
                        let mut perm_score = 0.0f64;
                        for l in 0..k {
                            perm_score += uf[(i, l)] * t[l];
                        }
                        if perm_score > orig_col[i] {
                            acc[i] += 1;
                        }
                    }

                    acc
                },
            )
            .reduce(
                || vec![0u32; n],
                |mut a, b| {
                    for i in 0..n {
                        a[i] += b[i];
                    }
                    a
                },
            );

        let denom = (n_perm + 1) as f64;
        for i in 0..n {
            result[(i, col)] = (counts[i] as f64 + 1.0) / denom;
        }
    }

    result
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    /// Path graph P3: 0 -- 1 -- 2
    fn make_p3_adjacency() -> CompressedSparseData {
        CompressedSparseData::new_csr(&[1.0, 1.0, 1.0, 1.0], &[1, 0, 2, 1], &[0, 1, 3, 4], (3, 3))
    }

    #[test]
    fn test_unnormalised_laplacian_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let dense = lap.to_dense();

        // Expected: [[1,-1,0],[-1,2,-1],[0,-1,1]]
        let expected = Mat::from_fn(3, 3, |i, j| match (i, j) {
            (0, 0) => 1.0,
            (0, 1) | (1, 0) => -1.0,
            (1, 1) => 2.0,
            (1, 2) | (2, 1) => -1.0,
            (2, 2) => 1.0,
            _ => 0.0,
        });

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (dense[(i, j)] - expected[(i, j)]).abs() < TOL,
                    "L[{},{}] = {}, expected {}",
                    i,
                    j,
                    dense[(i, j)],
                    expected[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_normalised_laplacian_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, true);
        let dense = lap.to_dense();

        let s = 1.0 / 2.0_f64.sqrt(); // 1/sqrt(2)

        // Expected: [[1, -1/sqrt(2), 0], [-1/sqrt(2), 1, -1/sqrt(2)], [0, -1/sqrt(2), 1]]
        let expected = Mat::from_fn(3, 3, |i, j| match (i, j) {
            (0, 0) | (1, 1) | (2, 2) => 1.0,
            (0, 1) | (1, 0) | (1, 2) | (2, 1) => -s,
            _ => 0.0,
        });

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (dense[(i, j)] - expected[(i, j)]).abs() < TOL,
                    "L_norm[{},{}] = {}, expected {}",
                    i,
                    j,
                    dense[(i, j)],
                    expected[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_full_spectral_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);

        // Eigenvalues of P3 unnormalised Laplacian: 0, 1, 3
        let expected_eigvals = [0.0, 1.0, 3.0];
        assert_eq!(decomp.eigenvalues.len(), 3);
        for (got, &want) in decomp.eigenvalues.iter().zip(&expected_eigvals) {
            assert!((got - want).abs() < 1e-8, "eigenvalue {} != {}", got, want);
        }

        // Verify U^T * U = I (orthonormality)
        let u = decomp.eigenvectors.as_ref();
        let mut utu = Mat::zeros(3, 3);
        matmul(
            utu.as_mut(),
            Accum::Replace,
            u.transpose(),
            u,
            1.0,
            default_par(),
        );
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (utu[(i, j)] - expected).abs() < 1e-8,
                    "U^T U [{},{}] = {}, expected {}",
                    i,
                    j,
                    utu[(i, j)],
                    expected
                );
            }
        }
    }

    #[test]
    fn test_regularised_laplacian_kernel_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);

        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );
        let kernel = materialise_kernel(&decomp, &f);

        // Known: (L + I)^{-1} = (1/8) * [[5,2,1],[2,4,2],[1,2,5]]
        let expected = Mat::from_fn(3, 3, |i, j| {
            let m = [[5.0, 2.0, 1.0], [2.0, 4.0, 2.0], [1.0, 2.0, 5.0]];
            m[i][j] / 8.0
        });

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (kernel[(i, j)] - expected[(i, j)]).abs() < 1e-8,
                    "K[{},{}] = {}, expected {}",
                    i,
                    j,
                    kernel[(i, j)],
                    expected[(i, j)]
                );
            }
        }
    }

    #[test]
    fn test_diffusion_kernel_transform() {
        // Just check the scalar transform
        let eigvals = vec![0.0, 1.0, 3.0];
        let f = kernel_eigenvalues(&eigvals, &KernelType::Diffusion { sigma2: 2.0 });

        assert!((f[0] - 1.0).abs() < TOL); // exp(0) = 1
        assert!((f[1] - (-1.0_f64).exp()).abs() < TOL); // exp(-1)
        assert!((f[2] - (-3.0_f64).exp()).abs() < TOL); // exp(-3)
    }

    #[test]
    fn test_pstep_kernel_transform() {
        let eigvals = vec![0.0, 1.0, 2.0];
        let f = kernel_eigenvalues(&eigvals, &KernelType::PStep { a: 2.0, p: 3 });

        assert!((f[0] - 8.0).abs() < TOL); // (2-0)^3
        assert!((f[1] - 1.0).abs() < TOL); // (2-1)^3
        assert!((f[2] - 0.0).abs() < TOL); // (2-2)^3
    }

    #[test]
    fn test_diffuse_raw_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        // input: node 0 is positive, all 3 nodes in background
        let input = Mat::from_fn(3, 1, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let bkgd = vec![0, 1, 2];

        let score = diffuse_raw(&decomp, &f, input.as_ref(), &bkgd);

        // Expected: column 0 of K = [5/8, 2/8, 1/8]
        assert!((score[(0, 0)] - 0.625).abs() < 1e-8);
        assert!((score[(1, 0)] - 0.25).abs() < 1e-8);
        assert!((score[(2, 0)] - 0.125).abs() < 1e-8);
    }

    #[test]
    fn test_diffuse_raw_subset_background() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        // background is only nodes 0 and 1; node 0 positive
        let input = Mat::from_fn(2, 1, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let bkgd = vec![0, 1];

        let score = diffuse_raw(&decomp, &f, input.as_ref(), &bkgd);

        // This is K[:, [0,1]] * [1, 0] = K[:, 0]
        // Same as full-background test above
        assert!((score[(0, 0)] - 0.625).abs() < 1e-8);
        assert!((score[(1, 0)] - 0.25).abs() < 1e-8);
        assert!((score[(2, 0)] - 0.125).abs() < 1e-8);
    }

    #[test]
    fn test_diffuse_z_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        let input = Mat::from_fn(3, 1, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let bkgd = vec![0, 1, 2];

        let z = diffuse_z(&decomp, &f, input.as_ref(), &bkgd);

        // node 0 (positive) should have highest z-score
        assert!(z[(0, 0)] > z[(1, 0)]);
        assert!(z[(0, 0)] > z[(2, 0)]);

        // node 0 should be positive (above mean), nodes 1 and 2 negative
        assert!(z[(0, 0)] > 0.0);
        assert!(z[(1, 0)] < 0.0);
        assert!(z[(2, 0)] < 0.0);

        // z[1] should be analytically -1/sqrt(2)
        assert!(
            (z[(1, 0)] - (-1.0 / 2.0_f64.sqrt())).abs() < 1e-8,
            "z[1] = {}, expected {}",
            z[(1, 0)],
            -1.0 / 2.0_f64.sqrt()
        );

        // all values should be finite
        for i in 0..3 {
            assert!(z[(i, 0)].is_finite(), "z[{}] is not finite", i);
        }
    }

    #[test]
    fn test_diffuse_mc_p3() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        let input = Mat::from_fn(3, 1, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let bkgd = vec![0, 1, 2];

        let hr = diffuse_mc(&decomp, &f, input.as_ref(), &bkgd, 10_000, 42);

        // heatranks in (0, 1)
        for i in 0..3 {
            assert!(hr[(i, 0)] > 0.0 && hr[(i, 0)] < 1.0);
        }

        // source node warmest
        assert!(
            hr[(0, 0)] < hr[(1, 0)],
            "source heatrank {} should be < node 1 heatrank {}",
            hr[(0, 0)],
            hr[(1, 0)]
        );
        assert!(
            hr[(0, 0)] < hr[(2, 0)],
            "source heatrank {} should be < node 2 heatrank {}",
            hr[(0, 0)],
            hr[(2, 0)]
        );

        // non-source nodes: only test they're worse than source.
        // strict ordering between them is not guaranteed by heatrank
        // (graph symmetry can produce identical exceedance counts)
        assert!(hr[(1, 0)] <= hr[(2, 0)]);
    }

    #[test]
    fn test_diffuse_mc_reproducibility() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        let input = Mat::from_fn(3, 1, |i, _| if i == 0 { 1.0 } else { 0.0 });
        let bkgd = vec![0, 1, 2];

        let hr1 = diffuse_mc(&decomp, &f, input.as_ref(), &bkgd, 1000, 42);
        let hr2 = diffuse_mc(&decomp, &f, input.as_ref(), &bkgd, 1000, 42);

        for i in 0..3 {
            assert_eq!(
                hr1[(i, 0)],
                hr2[(i, 0)],
                "non-deterministic result at node {}",
                i
            );
        }
    }

    #[test]
    fn test_raw_matches_materialised_kernel() {
        let adj = make_p3_adjacency();
        let lap = laplacian(&adj, false);
        let decomp = spectral_decompose(&lap, false, SpectralStrategy::Full);
        let f = kernel_eigenvalues(
            &decomp.eigenvalues,
            &KernelType::RegularisedLaplacian {
                sigma2: 1.0,
                add_diag: 1.0,
            },
        );

        let kernel = materialise_kernel(&decomp, &f);

        // try every unit vector as input
        let bkgd = vec![0, 1, 2];
        for src in 0..3 {
            let input = Mat::from_fn(3, 1, |i, _| if i == src { 1.0 } else { 0.0 });
            let factored = diffuse_raw(&decomp, &f, input.as_ref(), &bkgd);

            for i in 0..3 {
                assert!(
                    (factored[(i, 0)] - kernel[(i, src)]).abs() < 1e-10,
                    "mismatch at ({}, {}) for source {}",
                    i,
                    0,
                    src
                );
            }
        }
    }
}
