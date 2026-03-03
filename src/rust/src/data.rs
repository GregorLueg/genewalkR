use extendr_api::*;
use rand::seq::SliceRandom;
use rand::{rngs::StdRng, Rng, SeedableRng};
use rustc_hash::{FxHashMap, FxHashSet};

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
        "stochastic" | "stochastic_block" => Some(Node2VecDataType::Stochastic),
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
    p_between: f64,
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
            if rng.random::<f64>() < p_between {
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

/////////////////////////////
// GeneWalk synthetic data //
/////////////////////////////

/// Structure to store synthetic GeneWalk graph data for testing
///
/// ### Fields
///
/// * `ontology_edges` - Parent-child edges within the ontology graph
/// * `gene_ont_edges` - Edges linking genes to their ontology term annotations
/// * `ppi_edges` - Protein-protein interaction edges between genes
/// * `signal_genes` - Names of signal genes, annotated within a single ontology
///   subtree
/// * `noise_genes` - Names of noise genes, annotated across multiple ontology
///   subtrees
pub struct SyntheticGeneWalkData {
    pub ontology_edges: Vec<(String, String)>,
    pub gene_ont_edges: Vec<(String, String)>,
    pub ppi_edges: Vec<(String, String)>,
    pub signal_genes: Vec<String>,
    pub noise_genes: Vec<String>,
}

/// Build a synthetic GeneWalk dataset with controlled signal and noise
///
/// ### Params
///
/// * `n_signal_genes` - Number of signal genes, annotated within a single
///   ontology subtree.
/// * `n_noise_genes` - Number of noise genes, annotations spread across
///   multiple subtrees.
/// * `n_roots` - Number of root nodes in the ontology, each forming its own
///   subtree.
/// * `depth` - Depth of each ontology subtree.
/// * `branching` - Base branching factor per node in the ontology.
/// * `p_lateral` - Probability of a lateral edge between sibling nodes at the
///   same depth.
/// * `p_ppi` - Probability of a PPI edge between any two genes within the same
///   group.
/// * `min_annotations` - Minimum number of ontology term annotations per gene.
/// * `max_annotations` - Maximum number of ontology term annotations per gene.
/// * `min_noise_subtrees` - Minimum number of distinct ontology subtrees a
///   noise gene must span.
/// * `seed` - RNG seed for reproducibility
///
/// ### Returns
///
/// `SyntheticGeneWalkData` containing ontology edges, gene-ontology annotation
/// edges, PPI edges, and the signal and noise gene name lists
#[allow(clippy::too_many_arguments)]
pub fn build_synthetic_genewalk(
    n_signal_genes: usize,
    n_noise_genes: usize,
    n_roots: usize,
    depth: usize,
    branching: usize,
    p_lateral: f64,
    p_ppi: f64,
    min_annotations: usize,
    max_annotations: usize,
    min_noise_subtrees: usize,
    seed: u64,
) -> SyntheticGeneWalkData {
    let mut rng = StdRng::seed_from_u64(seed);

    // build ontology: levels[subtree][depth] = node names at that depth
    let mut levels: Vec<Vec<Vec<String>>> = Vec::new();
    let mut ontology_edges: Vec<(String, String)> = Vec::new();
    let mut term_id = 1usize;

    for _ in 0..n_roots {
        let mut subtree_levels: Vec<Vec<String>> = Vec::new();
        let root = format!("term_{:04}", term_id);
        term_id += 1;
        subtree_levels.push(vec![root]);

        for d in 1..=depth {
            let parents = subtree_levels[d - 1].clone();
            let mut current_level: Vec<String> = Vec::new();

            for parent in &parents {
                let n_children = branching - 1 + rng.random_range(0..=2);
                for _ in 0..n_children {
                    let child = format!("term_{:04}", term_id);
                    term_id += 1;
                    ontology_edges.push((parent.clone(), child.clone()));
                    current_level.push(child);
                }
            }

            for i in 0..current_level.len() {
                for j in (i + 1)..current_level.len() {
                    if rng.random::<f64>() < p_lateral {
                        ontology_edges.push((current_level[i].clone(), current_level[j].clone()));
                    }
                }
            }

            if !current_level.is_empty() {
                subtree_levels.push(current_level);
            }
        }

        levels.push(subtree_levels);
    }

    // collect all terms per subtree
    let all_terms_per_subtree: Vec<Vec<String>> = levels
        .iter()
        .map(|subtree_levels| subtree_levels.iter().flatten().cloned().collect())
        .collect();

    let all_terms: Vec<String> = all_terms_per_subtree.iter().flatten().cloned().collect();

    // map each term back to its subtree index for noise gene validation
    let mut term_to_subtree: FxHashMap<String, usize> = FxHashMap::default();
    for (subtree_idx, terms) in all_terms_per_subtree.iter().enumerate() {
        for term in terms {
            term_to_subtree.insert(term.clone(), subtree_idx);
        }
    }

    let mut gene_ont_edges: Vec<(String, String)> = Vec::new();
    let mut ppi_edges: Vec<(String, String)> = Vec::new();

    // signal genes -> concentrated strictly in Subtree 0
    let signal_subtree = &all_terms_per_subtree[0];
    let mut signal_genes: Vec<String> = Vec::new();

    for i in 0..n_signal_genes {
        let gene = format!("gene_signal_{:04}", i + 1);
        let n_ann = min_annotations + rng.random_range(0..=(max_annotations - min_annotations));
        let mut seen = FxHashSet::default();

        for _ in 0..n_ann {
            let term = &signal_subtree[rng.random_range(0..signal_subtree.len())];
            if seen.insert(term.clone()) {
                gene_ont_edges.push((gene.clone(), term.clone()));
            }
        }

        signal_genes.push(gene);
    }

    // dense PPIs for the signal community
    for i in 0..signal_genes.len() {
        for j in (i + 1)..signal_genes.len() {
            if rng.random::<f64>() < p_ppi {
                ppi_edges.push((signal_genes[i].clone(), signal_genes[j].clone()));
            }
        }
    }

    // we distribute a pool of anchor genes equal to `n_signal_genes`
    // across the remaining subtrees so they don't blow up the graph size.
    if n_roots > 1 {
        let anchors_per_st = (n_signal_genes / (n_roots - 1)).max(2);

        #[allow(clippy::needless_range_loop)]
        for st in 1..n_roots {
            let st_terms = &all_terms_per_subtree[st];
            let mut local_anchors = Vec::new();

            for i in 0..anchors_per_st {
                let gene = format!("gene_anchor_{}_{:04}", st, i + 1);
                let n_ann =
                    min_annotations + rng.random_range(0..=(max_annotations - min_annotations));
                let mut seen = FxHashSet::default();

                for _ in 0..n_ann {
                    let term = &st_terms[rng.random_range(0..st_terms.len())];
                    if seen.insert(term.clone()) {
                        gene_ont_edges.push((gene.clone(), term.clone()));
                    }
                }
                local_anchors.push(gene);
            }

            // PPIs for local anchors
            for i in 0..local_anchors.len() {
                for j in (i + 1)..local_anchors.len() {
                    if rng.random::<f64>() < p_ppi {
                        ppi_edges.push((local_anchors[i].clone(), local_anchors[j].clone()));
                    }
                }
            }
        }
    }

    // noise genes -> spread randomly across multiple subtrees
    let n_subtrees = all_terms_per_subtree.len();
    let effective_min_subtrees = min_noise_subtrees.min(n_subtrees);
    let mut noise_genes: Vec<String> = Vec::new();

    for i in 0..n_noise_genes {
        let gene = format!("gene_noise_{:04}", i + 1);
        let n_ann = min_annotations + rng.random_range(0..=(max_annotations - min_annotations));
        let mut seen = FxHashSet::default();

        let mut forced_subtrees: Vec<usize> = (0..n_subtrees).collect();
        forced_subtrees.shuffle(&mut rng);
        forced_subtrees.truncate(effective_min_subtrees);

        for &st in &forced_subtrees {
            let terms = &all_terms_per_subtree[st];
            let term = &terms[rng.random_range(0..terms.len())];
            if seen.insert(term.clone()) {
                gene_ont_edges.push((gene.clone(), term.clone()));
            }
        }

        // fill remaining annotations randomly across all terms
        let remaining = n_ann.saturating_sub(seen.len());
        for _ in 0..remaining {
            let term = &all_terms[rng.random_range(0..all_terms.len())];
            if seen.insert(term.clone()) {
                gene_ont_edges.push((gene.clone(), term.clone()));
            }
        }

        noise_genes.push(gene);
    }

    // random PPIs across noise genes (no community structure)
    for i in 0..noise_genes.len() {
        for j in (i + 1)..noise_genes.len() {
            if rng.random::<f64>() < p_ppi {
                ppi_edges.push((noise_genes[i].clone(), noise_genes[j].clone()));
            }
        }
    }

    SyntheticGeneWalkData {
        ontology_edges,
        gene_ont_edges,
        ppi_edges,
        signal_genes,
        noise_genes,
    }
}

////////////////////////
// Differential graph //
////////////////////////

/// A pair of graphs for differential node embedding analysis
///
/// Holds two graphs that are largely topologically identical but differ in
/// defined regions, together with ground truth labels for which shared nodes
/// are expected to show context drift between graphs.
///
/// ### Fields
///
/// * `from_g1` - Edge origins in graph 1
/// * `to_g1` - Edge destinations in graph 1
/// * `from_g2` - Edge origins in graph 2
/// * `to_g2` - Edge destinations in graph 2
/// * `nodes_g1` - All node names in graph 1
/// * `cluster_g1` - Cluster membership for each node in graph 1
/// * `nodes_g2` - All node names in graph 2
/// * `cluster_g2` - Cluster membership for each node in graph 2
/// * `shared_nodes` - Nodes present in both graphs, defines the index for
///   `is_differential`
/// * `is_differential` - Ground truth flag indexed to `shared_nodes`; true
///   where topology changes
/// * `g1_only` - Nodes exclusive to graph 1
/// * `g2_only` - Nodes exclusive to graph 2
pub struct DifferentialGraphData {
    from_g1: Vec<String>,
    to_g1: Vec<String>,
    from_g2: Vec<String>,
    to_g2: Vec<String>,
    nodes_g1: Vec<String>,
    cluster_g1: Vec<usize>,
    nodes_g2: Vec<String>,
    cluster_g2: Vec<usize>,
    shared_nodes: Vec<String>,
    is_differential: Vec<bool>,
    g1_only: Vec<String>,
    g2_only: Vec<String>,
}

/// R differential graph data
///
/// ### Fields
///
/// * `0` - from, to edges for G1
/// * `1` - nodes of G1 and their community membership
/// * `2` - from, to edges for G2
/// * `3` - nodes of G2 and their community membership
/// * `4` - shared data with the shared nodes, which ones are different and
///   which ones are unique per graph
pub type RDifferentialGraphData = (List, List, List, List, List);

impl DifferentialGraphData {
    /// Create a new instance
    ///
    /// ### Params
    ///
    /// * `from_g1` - Edge origins in graph 1
    /// * `to_g1` - Edge destinations in graph 1
    /// * `from_g2` - Edge origins in graph 2
    /// * `to_g2` - Edge destinations in graph 2
    /// * `nodes_g1` - All node names in graph 1
    /// * `cluster_g1` - Cluster membership for graph 1 nodes
    /// * `nodes_g2` - All node names in graph 2
    /// * `cluster_g2` - Cluster membership for graph 2 nodes
    /// * `shared_nodes` - Intersection of both node sets
    /// * `is_differential` - Ground truth indexed to `shared_nodes`
    /// * `g1_only` - Nodes exclusive to graph 1
    /// * `g2_only` - Nodes exclusive to graph 2
    ///
    /// ### Returns
    ///
    /// New `DifferentialGraphData` instance
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        from_g1: Vec<String>,
        to_g1: Vec<String>,
        from_g2: Vec<String>,
        to_g2: Vec<String>,
        nodes_g1: Vec<String>,
        cluster_g1: Vec<usize>,
        nodes_g2: Vec<String>,
        cluster_g2: Vec<usize>,
        shared_nodes: Vec<String>,
        is_differential: Vec<bool>,
        g1_only: Vec<String>,
        g2_only: Vec<String>,
    ) -> Self {
        Self {
            from_g1,
            to_g1,
            from_g2,
            to_g2,
            nodes_g1,
            cluster_g1,
            nodes_g2,
            cluster_g2,
            shared_nodes,
            is_differential,
            g1_only,
            g2_only,
        }
    }

    /// Generate R lists from the differential graph data
    ///
    /// ### Returns
    ///
    /// Tuple of `(edges_g1, nodes_g1, edges_g2, nodes_g2, ground_truth)` where
    /// each element is a named R list. `ground_truth` contains `shared_nodes`,
    /// `is_differential`, `g1_only`, and `g2_only`.
    pub fn generate_lists(self) -> RDifferentialGraphData {
        (
            list!(from = self.from_g1, to = self.to_g1),
            list!(node = self.nodes_g1, cluster = self.cluster_g1),
            list!(from = self.from_g2, to = self.to_g2),
            list!(node = self.nodes_g2, cluster = self.cluster_g2),
            list!(
                shared_nodes = self.shared_nodes,
                is_differential = self.is_differential,
                g1_only = self.g1_only,
                g2_only = self.g2_only,
            ),
        )
    }
}

/// Add all pairwise edges for a clique to the edge lists
///
/// ### Params
///
/// * `from` - Edge origin list to append to
/// * `to` - Edge destination list to append to
/// * `nodes` - Nodes forming the clique
fn add_clique(from: &mut Vec<String>, to: &mut Vec<String>, nodes: &[String]) {
    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            from.push(nodes[i].clone());
            to.push(nodes[j].clone());
        }
    }
}

/// Generate a synthetic differential graph pair for testing context-aware node
/// embeddings
///
/// Produces two graphs with three communities. Community 1 is a negative
/// control and is topologically identical across both graphs. Community 2
/// contains a hub node that is demoted to a peripheral node in graph 2. Two
/// bridge nodes span community 2 and community 3 in graph 1, but are fully
/// embedded within community 3 in graph 2. A set of exclusive nodes appears in
/// only one graph.
///
/// Cluster assignments: 0 = exclusive, 1 = stable community, 2 = community 2,
/// 3 = community 3.
///
/// ### Params
///
/// * `n_stable` - Number of nodes in the stable negative-control community
/// * `n_comm2` - Number of regular nodes in community 2, excluding the hub
/// * `n_comm3` - Number of nodes in community 3, excluding bridge nodes
/// * `n_exclusive` - Number of G1-only and G2-only peripheral nodes each
///
/// ### Returns
///
/// `DifferentialGraphData`
pub fn differential_graph_synthetic(
    n_stable: usize,
    n_comm2: usize,
    n_comm3: usize,
    n_exclusive: usize,
) -> DifferentialGraphData {
    let stable: Vec<String> = (0..n_stable).map(|i| format!("stable_{:03}", i)).collect();
    let comm2: Vec<String> = (0..n_comm2).map(|i| format!("comm2_{:03}", i)).collect();
    let hub = "hub_gene".to_string();
    let comm3: Vec<String> = (0..n_comm3).map(|i| format!("comm3_{:03}", i)).collect();
    let bridge_1 = "bridge_001".to_string();
    let bridge_2 = "bridge_002".to_string();
    let g1_excl: Vec<String> = (0..n_exclusive)
        .map(|i| format!("g1_only_{:03}", i))
        .collect();
    let g2_excl: Vec<String> = (0..n_exclusive)
        .map(|i| format!("g2_only_{:03}", i))
        .collect();

    let mut from_g1: Vec<String> = Vec::new();
    let mut to_g1: Vec<String> = Vec::new();
    let mut from_g2: Vec<String> = Vec::new();
    let mut to_g2: Vec<String> = Vec::new();

    // community 1: stable clique, identical in both graphs
    add_clique(&mut from_g1, &mut to_g1, &stable);
    add_clique(&mut from_g2, &mut to_g2, &stable);

    // community 2: regular clique nodes
    add_clique(&mut from_g1, &mut to_g1, &comm2);
    add_clique(&mut from_g2, &mut to_g2, &comm2);

    // hub: fully connected within community 2 in G1, peripheral (one edge) in G2
    for node in &comm2 {
        from_g1.push(hub.clone());
        to_g1.push(node.clone());
    }
    if let Some(first) = comm2.first() {
        from_g2.push(hub.clone());
        to_g2.push(first.clone());
    }

    // community 3: stable clique, identical in both graphs
    add_clique(&mut from_g1, &mut to_g1, &comm3);
    add_clique(&mut from_g2, &mut to_g2, &comm3);

    // bridge nodes in G1: span community 2 (including hub) and community 3
    let n_bridge_comm2 = (n_comm2 / 2).max(1).min(n_comm2);
    let n_bridge_comm3 = (n_comm3 / 2).max(1).min(n_comm3);
    for bridge in &[&bridge_1, &bridge_2] {
        for i in 0..n_bridge_comm2 {
            from_g1.push((*bridge).clone());
            to_g1.push(comm2[i].clone());
        }
        from_g1.push((*bridge).clone());
        to_g1.push(hub.clone());
        for i in 0..n_bridge_comm3 {
            from_g1.push((*bridge).clone());
            to_g1.push(comm3[i].clone());
        }
    }
    from_g1.push(bridge_1.clone());
    to_g1.push(bridge_2.clone());

    // bridge nodes in G2: fully embedded within community 3
    for bridge in &[&bridge_1, &bridge_2] {
        for node in &comm3 {
            from_g2.push((*bridge).clone());
            to_g2.push(node.clone());
        }
    }
    from_g2.push(bridge_1.clone());
    to_g2.push(bridge_2.clone());

    // exclusive peripheral nodes
    for node in &g1_excl {
        from_g1.push(node.clone());
        to_g1.push(stable[0].clone());
    }
    for node in &g2_excl {
        from_g2.push(node.clone());
        to_g2.push(comm2.first().unwrap_or(&hub).clone());
    }

    // shared nodes and ground truth, ordered: stable → comm2 → hub → comm3 → bridges
    let mut shared_nodes: Vec<String> = Vec::new();
    shared_nodes.extend(stable.iter().cloned());
    shared_nodes.extend(comm2.iter().cloned());
    shared_nodes.push(hub.clone());
    shared_nodes.extend(comm3.iter().cloned());
    shared_nodes.push(bridge_1.clone());
    shared_nodes.push(bridge_2.clone());

    let mut is_differential: Vec<bool> = Vec::new();
    is_differential.extend(vec![false; n_stable]);
    is_differential.extend(vec![false; n_comm2]);
    is_differential.push(true); // hub
    is_differential.extend(vec![false; n_comm3]);
    is_differential.push(true); // bridge_1
    is_differential.push(true); // bridge_2

    // node lists with cluster assignments; exclusive nodes get cluster 0
    let mut nodes_g1 = shared_nodes.clone();
    let mut cluster_g1: Vec<usize> = Vec::new();
    cluster_g1.extend(vec![1; n_stable]);
    cluster_g1.extend(vec![2; n_comm2 + 1]); // +1 for hub
    cluster_g1.extend(vec![3; n_comm3 + 2]); // +2 for bridges
    nodes_g1.extend(g1_excl.iter().cloned());
    cluster_g1.extend(vec![0; n_exclusive]);

    let mut nodes_g2 = shared_nodes.clone();
    let mut cluster_g2: Vec<usize> = Vec::new();
    cluster_g2.extend(vec![1; n_stable]);
    cluster_g2.extend(vec![2; n_comm2 + 1]);
    cluster_g2.extend(vec![3; n_comm3 + 2]);
    nodes_g2.extend(g2_excl.iter().cloned());
    cluster_g2.extend(vec![0; n_exclusive]);

    DifferentialGraphData::new(
        from_g1,
        to_g1,
        from_g2,
        to_g2,
        nodes_g1,
        cluster_g1,
        nodes_g2,
        cluster_g2,
        shared_nodes,
        is_differential,
        g1_excl,
        g2_excl,
    )
}
