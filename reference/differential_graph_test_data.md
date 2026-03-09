# Generate Synthetic Differential Graph Data

Creates a pair of graphs with known topological differences to verify
context-aware differential embedding analysis. The two graphs share a
common backbone but differ in defined regions:

- Community 1 is a stable negative control, identical across both
  graphs.

- Community 2 contains a hub node demoted to a peripheral node in graph
  2.

- Two bridge nodes span communities 2 and 3 in graph 1 but are fully
  embedded within community 3 in graph 2.

- A set of exclusive nodes appears in only one graph.

## Usage

``` r
differential_graph_test_data(
  n_stable = 50L,
  n_comm2 = 50L,
  n_comm3 = 50L,
  n_exclusive = 5L
)
```

## Arguments

- n_stable:

  Integer. Number of nodes in the stable negative-control community.

- n_comm2:

  Integer. Number of regular nodes in community 2, excluding the hub.

- n_comm3:

  Integer. Number of nodes in community 3, excluding bridge nodes.

- n_exclusive:

  Integer. Number of exclusive nodes per graph.

## Value

A named list containing:

- g1_edges - data.table with columns `from`, `to` for graph 1.

- g1_nodes - data.table with columns `node`, `cluster` for graph 1.

- g2_edges - data.table with columns `from`, `to` for graph 2.

- g2_nodes - data.table with columns `node`, `cluster` for graph 2.

- data_info - data.table with columns `node`, `is_differential`,
  `node_status` (`"shared"`, `"g1_only"`, `"g2_only"`).
