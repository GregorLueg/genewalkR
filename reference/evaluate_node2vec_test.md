# Evaluate node2vec test results

Compares embedding similarities against ground truth clusters.

## Usage

``` r
evaluate_node2vec_test(embeddings, node_labels)
```

## Arguments

- embeddings:

  Matrix of node embeddings (nodes x dimensions).

- node_labels:

  data.table with columns: node, cluster.

## Value

A list with:

- within_cluster_sim - Mean similarity within clusters.

- between_cluster_sim - Mean similarity between clusters.

- within_dist - The within cluster similarities as vector.

- between_dist - The between cluster similarities as vector.
