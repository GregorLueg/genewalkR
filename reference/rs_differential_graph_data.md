# Generate a synthetic differential graph pair for testing context-aware node embeddings

Produces two graphs sharing the same backbone but differing in defined
regions. Community 1 is a stable negative control, identical across both
graphs. Community 2 contains a hub node demoted to a peripheral node in
graph 2. Two bridge nodes span communities 2 and 3 in graph 1 but are
fully embedded within community 3 in graph 2. Exclusive nodes appear in
only one graph.

## Usage

``` r
rs_differential_graph_data(n_stable, n_comm2, n_comm3, n_exclusive)
```

## Arguments

- n_stable:

  Integer. Number of nodes in the stable negative-control community.

- n_comm2:

  Integer. Number of regular nodes in community 2, excluding the hub.

- n_comm3:

  Integer. Number of nodes in community 3, excluding bridge nodes.

- n_exclusive:

  Integer. Number of exclusive nodes in each graph.

## Value

A named list with: `g1_edges` (from, to), `g1_nodes` (node, cluster),
`g2_edges` (from, to), `g2_nodes` (node, cluster), `data_info`
(shared_nodes, is_differential, g1_only, g2_only).
