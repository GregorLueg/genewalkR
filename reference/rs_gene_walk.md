# Generate GeneWalk node embeddings

Uses a SIMD-accelerated CPU implementation of word2vec with negative
sampling to learn node representations from biased random walks
(node2vec).

## Usage

``` r
rs_gene_walk(
  from,
  to,
  weights,
  gene_walk_params,
  embd_dim,
  directed,
  seed,
  verbose
)
```

## Arguments

- from:

  Integer vector. Node indices for edge origins.

- to:

  Integer vector. Node indices for edge destinations.

- weights:

  Optional numeric vector. Edge weights, defaults to 1.0.

- gene_walk_params:

  Named list. Training parameters (p, q, walks_per_node, walk_length,
  num_workers, n_epochs, num_negatives, window_size, lr, dim).

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed.

- verbose:

  Boolean. Controls verbosity.

## Value

A numeric matrix of n_nodes x embedding dimensions.
