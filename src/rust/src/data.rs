use extendr_api::*;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

/////////////
// Helpers //
/////////////

/// Enum to describe the different synthetic data types
#[derive(Debug, Clone, Default)]
pub enum Node2VecDataType {
    /// Barbell-like graph
    #[default]
    Barbell,
    /// Cavemen-like graph
    Cavemen,
    /// Stochastic graph with more in-community connections
    Stochastic,
}

/// Parse the wanted data
///
/// ### Params
///
/// * `s` - String to parse
///
/// ### Returns
///
/// Option of the `Node2VecDataType`
pub fn parse_node2vec_data(s: &str) -> Option<Node2VecDataType> {
    match s.to_lowercase().as_str() {
        "barbell" => Some(Node2VecDataType::Barbell),
        "cavemen" => Some(Node2VecDataType::Cavemen),
        "stochastic" => Some(Node2VecDataType::Stochastic),
        _ => None,
    }
}

/// Structure to store Node2VecData
///
/// ### Fields
///
/// * `from` - Name of the from nodes
/// * `to` - Name of the to nodes
/// * `edge_type` - The edge type
/// * `nodes` - The names of the nodes in the graph generally
/// * `cluster` - Which cluster does a node belong to
pub struct Node2VecData {
    from: Vec<String>,
    to: Vec<String>,
    edge_type: Vec<String>,
    nodes: Vec<String>,
    cluster: Vec<usize>,
}

impl Node2VecData {
    /// Create a new instance
    ///
    /// ### Params
    ///
    /// * `from` - The from node names
    /// * `to` - The to node names
    /// * `edge_type` - Type of the edge
    /// * `nodes` - Node names
    /// * `cluster` - To which cluster they belong
    pub fn new(
        from: Vec<String>,
        to: Vec<String>,
        edge_type: Vec<String>,
        nodes: Vec<String>,
        cluster: Vec<usize>,
    ) -> Self {
        Self {
            from,
            to,
            edge_type,
            nodes,
            cluster,
        }
    }

    /// Generate the R lists
    ///
    /// ### Returns
    ///
    /// `(list with edge info, list with node info)`
    pub fn generate_lists(self) -> (List, List) {
        (
            list!(from = self.from, to = self.to, edge_type = self.edge_type),
            list!(node = self.nodes, cluster = self.cluster),
        )
    }
}

/// Generate barbell like graph data
///
/// ### Params
///
/// * `n_nodes_per_cluster` - Number of nodes per cluster
///
/// ### Returns
///
/// Synthetic data as `Node2VecData` with barbell-like graph with two
/// communities
pub fn node2vec_barbell(n_nodes_per_cluster: usize) -> Node2VecData {
    let n_total = n_nodes_per_cluster * 2;

    let node_ids: Vec<String> = (1..=n_total).map(|i| format!("node_{:04}", i)).collect();

    let mut from = Vec::with_capacity(n_total);
    let mut to = Vec::with_capacity(n_total);

    // first clusters
    for i in 0..(n_nodes_per_cluster - 1) {
        for j in (i + 1)..n_nodes_per_cluster {
            from.push(node_ids[i].clone());
            to.push(node_ids[j].clone());
        }
    }

    // second cluster
    for i in 0..(n_nodes_per_cluster - 1) {
        for j in (i + 1)..n_nodes_per_cluster {
            from.push(node_ids[n_nodes_per_cluster + i].clone());
            to.push(node_ids[n_nodes_per_cluster + j].clone());
        }
    }

    // bridge edge
    from.push(node_ids[n_nodes_per_cluster - 1].clone());
    to.push(node_ids[n_nodes_per_cluster].clone());

    let edge_type = vec!["connects".to_string(); from.len()];
    let cluster: Vec<usize> = (0..n_total)
        .map(|i| if i < n_nodes_per_cluster { 1 } else { 2 })
        .collect();

    Node2VecData::new(from, to, edge_type, node_ids, cluster)
}

/// Generate caveman like graph data
///
/// ### Params
///
/// * `n_nodes_per_cluster` - Number of nodes per cluster.
/// * `n_clusters` - Number of clusters.
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// Synthetic data as `Node2VecData` with barbell-like graph with two
/// communities
pub fn node2vec_caveman(
    n_nodes_per_cluster: usize,
    n_clusters: usize,
    seed: usize,
) -> Node2VecData {
    let n_total = n_nodes_per_cluster * n_clusters;

    let mut rng = StdRng::seed_from_u64(seed as u64);

    let node_ids: Vec<String> = (1..=n_total).map(|i| format!("node_{:04}", i)).collect();

    let mut from = Vec::with_capacity(n_total);
    let mut to = Vec::with_capacity(n_total);

    // Complete graphs within clusters
    for cluster in 0..n_clusters {
        let offset = cluster * n_nodes_per_cluster;
        for i in 0..(n_nodes_per_cluster - 1) {
            for j in (i + 1)..n_nodes_per_cluster {
                from.push(node_ids[offset + i].clone());
                to.push(node_ids[offset + j].clone());
            }
        }
    }

    // Sparse connections between clusters
    for c1 in 0..(n_clusters - 1) {
        for c2 in (c1 + 1)..n_clusters {
            if rng.random::<f64>() < 0.3 {
                let node1 = c1 * n_nodes_per_cluster + rng.random_range(0..n_nodes_per_cluster);
                let node2 = c2 * n_nodes_per_cluster + rng.random_range(0..n_nodes_per_cluster);
                from.push(node_ids[node1].clone());
                to.push(node_ids[node2].clone());
            }
        }
    }

    let edge_type = vec!["connects".to_string(); from.len()];
    let cluster: Vec<usize> = (0..n_total)
        .map(|i| (i / n_nodes_per_cluster) + 1)
        .collect();

    Node2VecData::new(from, to, edge_type, node_ids, cluster)
}

/// Generate stochastic block-like graph data
///
/// ### Params
///
/// * `n_nodes_per_cluster` - Number of nodes per cluster.
/// * `n_clusters` - Number of clusters.
/// * `p_within` - Probability to create within cluster edge.
/// * `p_between` - Probability to create between cluster edge.
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// Synthetic data as `Node2VecData` with barbell-like graph with two
/// communities
pub fn node2vec_stochastic_block(
    n_nodes_per_cluster: usize,
    n_clusters: usize,
    p_within: f64,
    p_between: f64,
    seed: usize,
) -> Node2VecData {
    let n_total = n_nodes_per_cluster * n_clusters;

    let mut rng = StdRng::seed_from_u64(seed as u64);

    let node_ids: Vec<String> = (1..=n_total).map(|i| format!("node_{:04}", i)).collect();

    let mut from = Vec::with_capacity(n_total);
    let mut to = Vec::with_capacity(n_total);

    for i in 0..(n_total - 1) {
        for j in (i + 1)..n_total {
            let cluster_i = i / n_nodes_per_cluster;
            let cluster_j = j / n_nodes_per_cluster;

            let prob = if cluster_i == cluster_j {
                p_within
            } else {
                p_between
            };

            if rng.random::<f64>() < prob {
                from.push(node_ids[i].clone());
                to.push(node_ids[j].clone());
            }
        }
    }

    let edge_type = vec!["connects".to_string(); from.len()];
    let cluster: Vec<usize> = (0..n_total)
        .map(|i| (i / n_nodes_per_cluster) + 1)
        .collect();

    Node2VecData::new(from, to, edge_type, node_ids, cluster)
}
