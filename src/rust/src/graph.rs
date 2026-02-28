use node2vec_rs::graph::{compute_transition_prob, Node2VecGraph};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};

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
