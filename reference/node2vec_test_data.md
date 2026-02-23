# Generate Simple Node2Vec Test Data

Creates graphs with known structure to verify node2vec is working
correctly. The following options exists:

- `"barbell"` - Two dense clusters connected by single edge.

- `"cavemen"` - Multiple isolated cliques with sparse inter-connections.

- `"stochastic_block"` - Communities with controlled edge probabilities.

## Usage

``` r
node2vec_test_data(
  test_type = c("barbell", "cavemen", "stochastic_block"),
  n_nodes_per_cluster = 10L,
  n_clusters = 3L,
  p_within = 0.8,
  p_between = 0.05,
  seed = 42L
)
```

## Arguments

- test_type:

  String. One of `c("barbell", "cavemen", "stochastic_block")`

- n_nodes_per_cluster:

  Integer. Nodes per cluster/community.

- n_clusters:

  Integer. Number of clusters.

- p_within:

  Numeric. Probability of edges within cluster (0-1).

- p_between:

  Numeric. Probability of edges between clusters (0-1).

- seed:

  Integer. Random seed.

## Value

A list containing:

- edges:

  data.table with from, to, edge_type

- node_labels:

  data.table with node, cluster
