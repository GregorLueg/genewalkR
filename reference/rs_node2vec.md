# Generate a node2vec embeddings

Trains node2vec on the original provided graph.

## Usage

``` r
rs_node2vec(
  from,
  to,
  weights,
  node2vec_params,
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

- node2vec_params:

  Named list. Training parameters (p, q, walks_per_node, walk_length,
  num_workers, n_epochs, num_negatives, window_size, lr, dim).

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed (incremented per rep).

- verbose:

  Boolean. Controls verbosity.

## Value

An embedding of dimension n_nodes x embedding dim.
