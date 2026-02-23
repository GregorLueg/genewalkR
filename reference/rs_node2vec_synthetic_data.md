# Generate synthetic data for node2vec

Generate synthetic data for node2vec

## Usage

``` r
rs_node2vec_synthetic_data(
  test_data,
  n_nodes_per_cluster,
  n_clusters,
  p_within,
  p_between,
  seed
)
```

## Arguments

- test_data:

  String. One of c("barbell", "caveman", "stochastic_block").

- n_nodes_per_cluster:

  Integer. Nodes per cluster.

- n_clusters:

  Integer. Number of clusters.

- p_within:

  Numeric. Within-cluster edge probability (0-1).

- p_between:

  Numeric. Between-cluster edge probability (0-1).

- seed:

  Integer. Random seed.

## Value

A list with edges and nodes.
