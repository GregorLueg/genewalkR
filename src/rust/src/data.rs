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

    // ---------------------------------------------------------
    // 1. SIGNAL GENES: Concentrated strictly in Subtree 0
    // ---------------------------------------------------------
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

    // Dense PPIs for the signal community
    for i in 0..signal_genes.len() {
        for j in (i + 1)..signal_genes.len() {
            if rng.random::<f64>() < p_ppi {
                ppi_edges.push((signal_genes[i].clone(), signal_genes[j].clone()));
            }
        }
    }

    // ---------------------------------------------------------
    // 2. BACKGROUND ANCHOR GENES: Prevent isolate/degree-1 terms
    // ---------------------------------------------------------
    // We distribute a pool of anchor genes equal to `n_signal_genes`
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

    // ---------------------------------------------------------
    // 3. NOISE GENES: Spread randomly across multiple subtrees
    // ---------------------------------------------------------
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

        // Fill remaining annotations randomly across all terms
        let remaining = n_ann.saturating_sub(seen.len());
        for _ in 0..remaining {
            let term = &all_terms[rng.random_range(0..all_terms.len())];
            if seen.insert(term.clone()) {
                gene_ont_edges.push((gene.clone(), term.clone()));
            }
        }

        noise_genes.push(gene);
    }

    // Random PPIs across noise genes (no community structure)
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
