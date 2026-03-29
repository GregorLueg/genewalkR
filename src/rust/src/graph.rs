//! Helpers for the generations of the underlying node2vec graph or helpers
//! to store graph structures.

use faer::{Mat, MatMut, MatRef};
use node2vec_rs::graph::{compute_transition_prob, Node2VecGraph};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};

//////////////
// node2vec //
//////////////

/// Prepare and deduplicate edges
///
/// De-duplicates edges (keeps the largest weight) and generates the data ready
/// for the random walks.
///
/// ### Params
///
/// * `from` - Vector indicating the node indices from where the edge originates
/// * `to` - Vector indicating the node indices from where to the edge goes
/// * `weight` - Optional weight vector for the edge. If not supplied, will
///   default to `1.0`.
///
/// ### Returns
///
/// Prepared edges `Vec<(from, to, weight)>`
pub fn prepare_edges(
    from: Vec<u32>,
    to: Vec<u32>,
    weights: Option<Vec<f64>>,
    directed: bool,
) -> Vec<(u32, u32, f32)> {
    assert_eq!(from.len(), to.len());

    let mut edge_map: FxHashMap<(u32, u32), f32> =
        FxHashMap::with_capacity_and_hasher(from.len(), FxBuildHasher);

    for i in 0..from.len() {
        let mut u = from[i];
        let mut v = to[i];
        if !directed && u > v {
            std::mem::swap(&mut u, &mut v);
        }
        let edge = (u, v);

        let weight = weights.as_ref().map_or(1.0, |w| w[i] as f32);

        if let Some(existing) = edge_map.get_mut(&edge) {
            if weight > *existing {
                *existing = weight;
            }
        } else {
            edge_map.insert(edge, weight);
        }
    }

    edge_map
        .into_iter()
        .map(|((f, t), weight)| (f, t, weight))
        .collect()
}

/// Generate a random graph via the configuration model (stub-pairing)
///
/// Matches the Python `nx.configuration_model` approach: extracts the degree
/// sequence from the input graph, creates stubs, shuffles, pairs them, then
/// removes self-loops and deduplicates multi-edges.
///
/// Node identities are new (0..n_nodes) -- there is no correspondence between
/// the original node IDs and the random graph's node IDs. This is intentional:
/// we only need the degree distribution, not node identity.
///
/// ### Params
///
/// * `edges` - Original (deduplicated) edges to extract degree sequence from
/// * `directed` - Whether the graph is directed
/// * `seed` - Random seed for reproducibility
///
/// ### Returns
///
/// Tuple of (deduplicated random edges, vocab_size)
pub fn generate_configuration_model(
    edges: &[(u32, u32, f32)],
    directed: bool,
    seed: u64,
) -> (Vec<(u32, u32, f32)>, usize) {
    // Compute degree sequence from the original graph
    let mut degree_map: FxHashMap<u32, usize> = FxHashMap::default();
    for (u, v, _) in edges {
        *degree_map.entry(*u).or_insert(0) += 1;
        if !directed {
            *degree_map.entry(*v).or_insert(0) += 1;
        }
    }

    // Sorted descending to match Python's sorted(..., reverse=True)
    let mut degree_seq: Vec<usize> = degree_map.values().cloned().collect();
    degree_seq.sort_unstable_by(|a, b| b.cmp(a));

    let n_nodes = degree_seq.len();

    // Create stubs: node i appears degree_seq[i] times
    let mut stubs: Vec<u32> = Vec::new();
    for (node_id, &deg) in degree_seq.iter().enumerate() {
        for _ in 0..deg {
            stubs.push(node_id as u32);
        }
    }

    // Shuffle stubs
    let mut rng = StdRng::seed_from_u64(seed);
    stubs.shuffle(&mut rng);

    // Pair consecutive stubs, skip self-loops and multi-edges
    let mut edge_set: FxHashSet<(u32, u32)> = FxHashSet::default();
    let mut result_edges: Vec<(u32, u32, f32)> = Vec::new();

    for chunk in stubs.chunks(2) {
        if chunk.len() < 2 {
            break;
        }
        let u = chunk[0];
        let v = chunk[1];

        // Skip self-loops
        if u == v {
            continue;
        }

        // Canonicalise for undirected
        let key = if !directed && u > v { (v, u) } else { (u, v) };

        // Skip multi-edges
        if edge_set.insert(key) {
            result_edges.push((key.0, key.1, 1.0));
        }
    }

    (result_edges, n_nodes)
}

/// Build an adjacency map of unique neighbours (excluding self-loops)
///
/// Used for sampling the null distribution: for each node, we compute cosine
/// similarity with each unique neighbour. Each edge contributes twice (once
/// from each endpoint), matching the Python implementation.
///
/// ### Params
///
/// * `edges` - Edge list (assumed deduplicated)
///
/// ### Returns
///
/// Map from node ID to set of unique neighbour IDs
pub fn build_adjacency(edges: &[(u32, u32, f32)]) -> FxHashMap<u32, FxHashSet<u32>> {
    let mut adj: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
    for (u, v, _) in edges {
        if u != v {
            adj.entry(*u).or_default().insert(*v);
            adj.entry(*v).or_default().insert(*u);
        }
    }
    adj
}

/// Generate the graph from edges
///
/// ### Params
///
/// * `edges` - Slice of tuples with `(from, to, weight)`.
/// * `p` - p parameter in node2vec that controls probability to return to
///   origin node.
/// * `q` - q parameter in node2vec that controls probability to reach out
///   further in the graph.
/// * `directed` - Boolean indicating if the graph is directed or reverse edges
///   shall be added
///
/// ### Returns
///
/// The `Node2VecGraph` graph.
pub fn create_graph(edges: &[(u32, u32, f32)], p: f32, q: f32, directed: bool) -> Node2VecGraph {
    let mut adjacency = FxHashMap::default();

    for (from, to, weight) in edges {
        adjacency
            .entry(*from)
            .or_insert_with(Vec::new)
            .push((*to, *weight));

        if !directed {
            adjacency
                .entry(*to)
                .or_insert_with(Vec::new)
                .push((*from, *weight));
        }
    }

    let transition_probs = compute_transition_prob(&adjacency, p, q);

    Node2VecGraph {
        adjacency,
        transition_probs,
    }
}

/// Extract the largest connected component from an edge list
///
/// Relabels nodes to 0..n so the vocab is contiguous.
///
/// ### Params
///
/// * `edge` - The edge tuples `(from, to, weight)`
/// * `n_nodes` - Number of nodes in the original graph
///
/// ### Returns
///
/// The largest component of the graph
pub fn largest_connected_component(
    edges: &[(u32, u32, f32)],
    n_nodes: usize,
) -> (Vec<(u32, u32, f32)>, usize) {
    // Union-find
    let mut parent: Vec<usize> = (0..n_nodes).collect();
    let mut rank: Vec<usize> = vec![0; n_nodes];

    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut Vec<usize>, rank: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra] < rank[rb] {
            parent[ra] = rb;
        } else if rank[ra] > rank[rb] {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra] += 1;
        }
    }

    for (u, v, _) in edges {
        union(&mut parent, &mut rank, *u as usize, *v as usize);
    }

    // Count component sizes
    let mut comp_size: FxHashMap<usize, usize> = FxHashMap::default();
    for i in 0..n_nodes {
        *comp_size.entry(find(&mut parent, i)).or_insert(0) += 1;
    }

    let largest_root = comp_size
        .iter()
        .max_by_key(|(_, &size)| size)
        .map(|(&root, _)| root)
        .unwrap_or(0);

    // Collect nodes in the largest component, build relabelling map
    let mut old_to_new: FxHashMap<u32, u32> = FxHashMap::default();
    let mut new_id = 0u32;
    for i in 0..n_nodes {
        if find(&mut parent, i) == largest_root {
            old_to_new.insert(i as u32, new_id);
            new_id += 1;
        }
    }

    let new_edges: Vec<(u32, u32, f32)> = edges
        .iter()
        .filter_map(|(u, v, w)| match (old_to_new.get(u), old_to_new.get(v)) {
            (Some(&nu), Some(&nv)) => Some((nu, nv, *w)),
            _ => None,
        })
        .collect();

    (new_edges, new_id as usize)
}

//////////////////////////
// CompressedSparseData //
//////////////////////////

/// Type to describe the CompressedSparseFormat
#[derive(Debug, Clone)]
pub enum CompressedSparseFormat {
    /// CSC-formatted data
    Csc,
    /// CSR-formatted data
    Csr,
}

impl CompressedSparseFormat {
    /// Returns boolean if it's CSC
    pub fn is_csc(&self) -> bool {
        matches!(self, CompressedSparseFormat::Csc)
    }
    /// Returns boolean if it's CSR
    pub fn is_csr(&self) -> bool {
        matches!(self, CompressedSparseFormat::Csr)
    }
}

/// Structure to store compressed sparse data of either type
#[derive(Debug, Clone)]
pub struct CompressedSparseData {
    /// Data values
    pub data: Vec<f64>,
    /// Indices of the data (pending format, column or row indices)
    pub indices: Vec<usize>,
    /// Index pointers
    pub indptr: Vec<usize>,
    /// Storage type, either `Csr` or `CSC`
    pub cs_type: CompressedSparseFormat,
    /// Dimensionality of the sparse data
    pub shape: (usize, usize),
}

impl CompressedSparseData {
    /// Generate a nes CSC version of the matrix
    ///
    /// ### Params
    ///
    /// * `data` - The underlying data
    /// * `indices` - The index positions (in this case row indices)
    /// * `indptr` - The index pointer (in this case the column index pointers)
    /// * `data2` - An optional second layer
    #[allow(dead_code)]
    pub fn new_csc(
        data: &[f64],
        indices: &[usize],
        indptr: &[usize],
        shape: (usize, usize),
    ) -> Self {
        Self {
            data: data.to_vec(),
            indices: indices.to_vec(),
            indptr: indptr.to_vec(),
            cs_type: CompressedSparseFormat::Csc,
            shape,
        }
    }

    /// Generate a nes CSR version of the matrix
    ///
    /// ### Params
    ///
    /// * `data` - The underlying data
    /// * `indices` - The index positions (in this case row indices)
    /// * `indptr` - The index pointer (in this case the column index pointers)
    /// * `data2` - An optional second layer
    pub fn new_csr(
        data: &[f64],
        indices: &[usize],
        indptr: &[usize],
        shape: (usize, usize),
    ) -> Self {
        Self {
            data: data.to_vec(),
            indices: indices.to_vec(),
            indptr: indptr.to_vec(),
            cs_type: CompressedSparseFormat::Csr,
            shape,
        }
    }

    /// Transform from CSC to CSR or vice versa
    ///
    /// ### Returns
    ///
    /// The transformed/transposed version
    pub fn transform(&self) -> Self {
        match self.cs_type {
            CompressedSparseFormat::Csr => {
                // CSR transposed is naturally CSC
                // Just reinterpret the data with swapped dimensions
                let as_csc = CompressedSparseData {
                    data: self.data.clone(),
                    indices: self.indices.clone(),
                    indptr: self.indptr.clone(),
                    cs_type: CompressedSparseFormat::Csc,
                    shape: (self.shape.1, self.shape.0), // swap dimensions
                };
                // Convert back to CSR for consistency
                csc_to_csr(&as_csc)
            }
            CompressedSparseFormat::Csc => {
                // CSC transposed is naturally CSR
                // Just reinterpret the data with swapped dimensions
                let as_csr = CompressedSparseData {
                    data: self.data.clone(),
                    indices: self.indices.clone(),
                    indptr: self.indptr.clone(),
                    cs_type: CompressedSparseFormat::Csr,
                    shape: (self.shape.1, self.shape.0), // swap dimensions
                };
                // Convert back to CSC for consistency
                csr_to_csc(&as_csr)
            }
        }
    }

    /// Transpose the matrix, maintaining the same storage format
    ///
    /// * For CSR: creates CSR^T (which is naturally CSC, then converted back to
    ///   CSR)
    /// * For CSC: creates CSC^T (which is naturally CSR, then converted back to
    ///   CSC)
    ///
    /// ### Returns
    ///
    /// Transposed matrix in the same storage format
    pub fn transpose(&self) -> Self {
        transpose_sparse(self)
    }

    /// Returns the shape of the matrix
    ///
    /// ### Returns
    ///
    /// A tuple of `(nrow, ncol)`
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Returns the NNZ
    ///
    /// ### Returns
    ///
    /// The number of NNZ
    pub fn get_nnz(&self) -> usize {
        self.data.len()
    }

    /// Returns the number of rows
    ///
    /// ### Returns
    ///
    /// The number of rows
    pub fn nrows(&self) -> usize {
        self.shape.0
    }

    /// Returns the number of columns
    ///
    /// ### Returns
    ///
    /// The number of columns
    pub fn ncols(&self) -> usize {
        self.shape.1
    }

    /// Converts the compressed sparse matrix to a dense matrix
    ///
    /// ### Returns
    ///
    /// The dense matrix
    pub fn to_dense(&self) -> Mat<f64> {
        let (nrows, ncols) = self.shape;
        let mut dense: Mat<f64> = Mat::zeros(nrows, ncols);

        match self.cs_type {
            CompressedSparseFormat::Csr => {
                for i in 0..nrows {
                    let start = self.indptr[i];
                    let end = self.indptr[i + 1];
                    for idx in start..end {
                        let j = self.indices[idx];
                        let val = self.data[idx];
                        dense[(i, j)] = val;
                    }
                }
            }
            CompressedSparseFormat::Csc => {
                for j in 0..ncols {
                    let start = self.indptr[j];
                    let end = self.indptr[j + 1];
                    for idx in start..end {
                        let i = self.indices[idx];
                        let val = self.data[idx];
                        dense[(i, j)] = val;
                    }
                }
            }
        }

        dense
    }
}

/// Transpose a compressed sparse matrix (CSC→CSR or CSR→CSC).
///
/// This is the standard two-pass sparse transpose in O(nnz) time.
///
/// ### Params
///
/// * `sparse_data`: The input compressed sparse matrix to be transformed.
///
/// ### Returns
///
/// The transposed compressed sparse matrix.
pub fn transpose_sparse(sparse_data: &CompressedSparseData) -> CompressedSparseData {
    let nnz = sparse_data.get_nnz();
    let (nrow, ncol) = sparse_data.shape();

    // the "minor" dimension is what becomes the new indptr axis.
    let (new_major, new_type) = match sparse_data.cs_type {
        CompressedSparseFormat::Csc => (nrow, CompressedSparseFormat::Csr),
        CompressedSparseFormat::Csr => (ncol, CompressedSparseFormat::Csc),
    };

    // first pass: count entries per new-major index
    let mut new_indptr = vec![0usize; new_major + 1];
    for &idx in &sparse_data.indices {
        new_indptr[idx + 1] += 1;
    }
    for i in 0..new_major {
        new_indptr[i + 1] += new_indptr[i];
    }

    // second pass: scatter data
    let mut new_data: Vec<f64> = vec![0_f64; nnz];
    let mut new_indices: Vec<usize> = vec![0usize; nnz];
    unsafe {
        new_data.set_len(nnz);
        new_indices.set_len(nnz);
    }

    // Reuse new_indptr as the write cursor — we'll restore it afterwards.
    // Work on a mutable window so we don't need a separate `next` vec.
    // We iterate old major indices and scatter into new positions.
    let old_major_len = sparse_data.indptr.len() - 1;
    for major in 0..old_major_len {
        for idx in sparse_data.indptr[major]..sparse_data.indptr[major + 1] {
            let minor = sparse_data.indices[idx];
            let pos = new_indptr[minor];

            // SAFETY: pos < nnz guaranteed by the counting pass
            unsafe {
                *new_data.get_unchecked_mut(pos) = sparse_data.data[idx];
                *new_indices.get_unchecked_mut(pos) = major;
            }

            new_indptr[minor] += 1;
        }
    }

    // restore new_indptr: the scatter pass shifted every entry forward by its
    // count, so we shift the whole array right by one position.
    for i in (1..=new_major).rev() {
        new_indptr[i] = new_indptr[i - 1];
    }
    new_indptr[0] = 0;

    CompressedSparseData {
        data: new_data,
        indices: new_indices,
        indptr: new_indptr,
        cs_type: new_type,
        shape: (nrow, ncol),
    }
}

/// Transforms a CompressedSparseData that is CSC to CSR
///
/// ### Params
///
/// * `sparse_data` - The CompressedSparseData you want to transform
///
/// ### Returns
///
///
pub fn csc_to_csr(sparse_data: &CompressedSparseData) -> CompressedSparseData {
    // early return if already in the desired format
    if sparse_data.cs_type.is_csr() {
        return sparse_data.clone();
    }

    let (nrow, _) = sparse_data.shape();
    let nnz = sparse_data.get_nnz();
    let mut row_ptr = vec![0; nrow + 1];

    for &r in &sparse_data.indices {
        row_ptr[r + 1] += 1;
    }

    for i in 0..nrow {
        row_ptr[i + 1] += row_ptr[i];
    }

    let mut csr_data = vec![0_f64; nnz];
    let mut csr_col_ind = vec![0; nnz];
    let mut next = row_ptr[..nrow].to_vec();

    for col in 0..(sparse_data.indptr.len() - 1) {
        for idx in sparse_data.indptr[col]..sparse_data.indptr[col + 1] {
            let row = sparse_data.indices[idx];
            let pos = next[row];

            csr_data[pos] = sparse_data.data[idx];
            csr_col_ind[pos] = col;

            next[row] += 1;
        }
    }

    CompressedSparseData {
        data: csr_data,
        indices: csr_col_ind,
        indptr: row_ptr,
        cs_type: CompressedSparseFormat::Csr,
        shape: sparse_data.shape(),
    }
}

/// Transform CSR stored data into CSC stored data
///
/// This version does a full memory copy of the data.
///
/// ### Params
///
/// * `sparse_data` - The CompressedSparseData you want to transform. Needs
///   to be in CSR format.
///
/// ### Returns
///
/// The data in CSC format, i.e., `CompressedSparseData`
pub fn csr_to_csc(sparse_data: &CompressedSparseData) -> CompressedSparseData {
    // early return if already in the desired format
    if sparse_data.cs_type.is_csc() {
        return sparse_data.clone();
    }

    let nnz = sparse_data.get_nnz();
    let (_, ncol) = sparse_data.shape();
    let mut col_ptr = vec![0; ncol + 1];

    // count occurrences per column
    for &c in &sparse_data.indices {
        col_ptr[c + 1] += 1;
    }

    // cumulative sum to get column pointers
    for i in 0..ncol {
        col_ptr[i + 1] += col_ptr[i];
    }

    let mut csc_data = vec![0_f64; nnz];
    let mut csc_row_ind = vec![0; nnz];
    let mut next = col_ptr[..ncol].to_vec();

    // iterate through rows and place data in CSC format
    for row in 0..(sparse_data.indptr.len() - 1) {
        for idx in sparse_data.indptr[row]..sparse_data.indptr[row + 1] {
            let col = sparse_data.indices[idx];
            let pos = next[col];

            csc_data[pos] = sparse_data.data[idx];
            csc_row_ind[pos] = row;

            next[col] += 1;
        }
    }

    CompressedSparseData {
        data: csc_data,
        indices: csc_row_ind,
        indptr: col_ptr,
        cs_type: CompressedSparseFormat::Csc,
        shape: sparse_data.shape(),
    }
}

/// Sparse CSR times dense matrix: out = A * x
///
/// ### Params
///
/// * `csr` - Matrix in CSR format
/// * `x` - Dense input, shape (A.ncols, m)
/// * `out` - Dense output, shape (A.nrows, m), overwritten
pub fn spmm_csr(csr: &CompressedSparseData, x: MatRef<f64>, mut out: MatMut<f64>) {
    let nrows = csr.nrows();
    let m = x.ncols();

    for i in 0..nrows {
        for col in 0..m {
            out[(i, col)] = 0.0;
        }
        for idx in csr.indptr[i]..csr.indptr[i + 1] {
            let j = csr.indices[idx];
            let v = csr.data[idx];
            for col in 0..m {
                out[(i, col)] += v * x[(j, col)];
            }
        }
    }
}
