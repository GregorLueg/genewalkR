use node2vec_rs::graph::{compute_transition_prob, Node2VecGraph};
use rustc_hash::FxHashMap;

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
