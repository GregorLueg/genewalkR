use extendr_api::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rustc_hash::{FxHashMap, FxHashSet};

/////////////
// Helpers //
/////////////

/// Results for the
///
/// * `0` - Edge list with (from, to)
/// * `1` - Vector of the pathway nodes
type DagResults = (Vec<(String, String)>, Vec<PathwayNode>);

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

/// Pathway node structure
///
/// ### Params
///
/// * `id` - Name of the pathway
/// * `subtree` - Identifier of the sub tree
/// * `depth` - Depth of the node
/// * `ancestors` - Name of the ancestors
#[derive(Clone)]
pub struct PathwayNode {
    pub id: String,
    pub subtree: usize,
    pub depth: usize,
    pub ancestors: Vec<String>,
}

//////////////
// Node2Vec //
//////////////

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
/// Synthetic data as `Node2VecData` with caveman-like graph with sparse connections between communities
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
/// Synthetic data as `Node2VecData` with stochastic block-like graph with probabilistic community structure
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

//////////////////////////
// GeneWalk pathway DAG //
//////////////////////////

/// Build the pathway DAG
///
/// ### Params
///
/// * `n_pathways` - Number of pathways to generate
/// * `pathway_depth` - The depth of the pathways
/// * `pathway_branching` - Number of branches
/// * `n_communities` - Number of communities
pub fn build_pathway_dag(
    n_pathways: usize,
    pathway_depth: usize,
    pathway_branching: usize,
    n_communities: usize,
) -> DagResults {
    let mut pathway_edges = Vec::new();
    let mut pathways = Vec::new();
    let mut pathway_ancestors: FxHashMap<String, Vec<String>> = FxHashMap::default();

    let n_roots = std::cmp::max(2, n_communities);
    let mut current_id = 1;
    let mut subtree_id = 1;

    for _ in 0..n_roots {
        if current_id > n_pathways {
            break;
        }

        let root_name = format!("pathway_{:04}", current_id);
        pathways.push(PathwayNode {
            id: root_name.clone(),
            subtree: subtree_id,
            depth: 0,
            ancestors: Vec::new(),
        });
        pathway_ancestors.insert(root_name.clone(), Vec::new());

        let mut current_level = vec![root_name];
        current_id += 1;

        for depth in 1..pathway_depth {
            if current_id > n_pathways {
                break;
            }

            let mut next_level = Vec::new();

            for parent in &current_level {
                let n_children = fastrand::usize(
                    std::cmp::max(1, pathway_branching - 1)..=pathway_branching + 1,
                );

                for _ in 0..n_children {
                    if current_id > n_pathways {
                        break;
                    }

                    let child_name = format!("pathway_{:04}", current_id);
                    let parent_ancestors = pathway_ancestors.get(parent).unwrap();
                    let mut child_ancestors = parent_ancestors.clone();
                    child_ancestors.push(parent.clone());

                    pathways.push(PathwayNode {
                        id: child_name.clone(),
                        subtree: subtree_id,
                        depth,
                        ancestors: child_ancestors.clone(),
                    });

                    pathway_edges.push((parent.clone(), child_name.clone()));
                    pathway_ancestors.insert(child_name.clone(), child_ancestors.clone());

                    next_level.push(child_name.clone());
                    current_id += 1;

                    // occasional second parent for DAG structure
                    if fastrand::f64() < 0.2 && current_level.len() > 1 {
                        let eligible: Vec<_> =
                            current_level.iter().filter(|p| *p != parent).collect();

                        if !eligible.is_empty() {
                            let alt_parent = eligible[fastrand::usize(..eligible.len())];
                            pathway_edges.push((alt_parent.clone(), child_name.clone()));

                            let alt_ancestors = pathway_ancestors.get(alt_parent).unwrap();
                            let mut combined: FxHashSet<String> =
                                child_ancestors.iter().cloned().collect();
                            combined.extend(alt_ancestors.iter().cloned());
                            combined.insert(alt_parent.clone());

                            let updated: Vec<_> = combined.into_iter().collect();
                            pathway_ancestors.insert(child_name.clone(), updated.clone());

                            if let Some(node) = pathways.iter_mut().find(|n| n.id == child_name) {
                                node.ancestors = updated;
                            }
                        }
                    }
                }
            }

            if next_level.is_empty() {
                break;
            }
            current_level = next_level;
        }
        subtree_id += 1;
    }

    (pathway_edges, pathways)
}

/// Generate the gene to pathway associations
///
/// ### Params
///
/// * `gene_ids` - Names of the genes
/// * `gene_communities` - The indices of gene to community
/// * `pathway_ids` - The name of the pathways
/// * `pathways` - Slice of the pathway nodes
/// * `community_focal_pathways` -
pub fn build_gene_pathway_associations(
    gene_ids: &[String],
    gene_communities: &[i32],
    pathway_ids: &[String],
    pathways: &[PathwayNode],
    community_focal_pathways: &FxHashMap<i32, Vec<String>>,
    signal_strength: f64,
    connections_per_gene: usize,
) -> Vec<(String, String)> {
    let mut edges = Vec::new();
    let pathway_ancestor_map: FxHashMap<_, _> = pathways
        .iter()
        .map(|p| (p.id.clone(), p.ancestors.clone()))
        .collect();

    for (gene_idx, gene) in gene_ids.iter().enumerate() {
        let gene_comm = gene_communities[gene_idx];
        let focal_pathways = community_focal_pathways
            .get(&gene_comm)
            .map(|v| v.as_slice())
            .unwrap_or(pathway_ids);

        let mut connected_pathways = FxHashSet::default();

        for _ in 0..connections_per_gene {
            let pathway = if fastrand::f64() < signal_strength && !focal_pathways.is_empty() {
                focal_pathways[fastrand::usize(..focal_pathways.len())].clone()
            } else {
                pathway_ids[fastrand::usize(..pathway_ids.len())].clone()
            };
            connected_pathways.insert(pathway);
        }

        for pathway in connected_pathways {
            edges.push((gene.clone(), pathway.clone()));

            if let Some(ancestors) = pathway_ancestor_map.get(&pathway) {
                for ancestor in ancestors {
                    edges.push((gene.clone(), ancestor.clone()));
                }
            }
        }
    }

    edges.sort();
    edges.dedup();
    edges
}
