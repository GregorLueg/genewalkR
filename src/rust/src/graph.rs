use node2vec_rs::graph::{compute_transition_prob, Node2VecGraph};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rustc_hash::{FxBuildHasher, FxHashMap};

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
pub fn prepare_edges(
    from: Vec<u32>,
    to: Vec<u32>,
    weights: Option<Vec<f64>>,
) -> Vec<(u32, u32, f32)> {
    assert_eq!(from.len(), to.len());

    if let Some(ref w) = weights {
        assert_eq!(from.len(), w.len());
    }

    let mut edge_map: FxHashMap<(u32, u32), f32> =
        FxHashMap::with_capacity_and_hasher(from.len(), FxBuildHasher);

    for i in 0..from.len() {
        let edge = (from[i], to[i]);
        let weight = weights.as_ref().map_or(1.0, |w| w[i] as f32);

        if let Some(existing) = edge_map.get_mut(&edge) {
            if weight > *existing {
                eprintln!(
                    "Warning: Duplicate edge ({}, {}) found. Keeping higher weight {:.4} over {:.4}",
                    edge.0, edge.1, weight, *existing
                );
                *existing = weight;
            }
        } else {
            edge_map.insert(edge, weight);
        }
    }

    edge_map
        .into_iter()
        .map(|((from, to), weight)| (from, to, weight))
        .collect()
}

/// Generates a random network while preserving node degree distribution
///
/// ### Params
///
/// * `edges` - Original edges to extract degree sequence from
/// * `directed` - Whether the graph is directed
/// * `seed` - Random seed for reproducibility
///
/// ### Returns
///
/// Random edges with the same degree distribution as the input
pub fn generate_random_network(
    edges: &[(u32, u32, f32)],
    directed: bool,
    seed: u64,
) -> Vec<(u32, u32, f32)> {
    let mut rng = StdRng::seed_from_u64(seed);

    let mut degrees = FxHashMap::default();

    for (from, to, _) in edges {
        *degrees.entry(*from).or_insert(0) += 1;
        if !directed {
            *degrees.entry(*to).or_insert(0) += 1;
        }
    }

    let mut stubs: Vec<u32> = Vec::new();
    for (node, degree) in degrees {
        for _ in 0..degree {
            stubs.push(node);
        }
    }

    stubs.shuffle(&mut rng);

    // pair up stubs
    let mut random_edges = Vec::new();
    for i in (0..stubs.len()).step_by(2) {
        if i + 1 < stubs.len() {
            let from = stubs[i];
            let to = stubs[i + 1];
            // no self-loops
            if from != to {
                random_edges.push((from, to, 1.0));
            }
        }
    }

    random_edges
}

/// Generate the graph from edges
///
/// ### Params
///
/// * `edges` - Slice of tuples with `(from, to, weight)`.
/// * `p` - p parameter in node2vec that controls probability to return to
///   origin node.
/// * `q` - q parameter in node2vec that controls probability to reach out
///   futher in the graph.
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
