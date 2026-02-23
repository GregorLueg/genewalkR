# Generate permuted embeddings for null distribution

Generates degree-preserving random networks and trains node2vec on each,
returning the raw embedding matrices for downstream statistical testing.

## Usage

``` r
rs_gene_walk_perm(
  from,
  to,
  weights,
  gene_walk_params,
  n_perm,
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

  Named list. Training parameters.

- n_perm:

  Integer. Number of permutations.

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed.

- verbose:

  Boolean. Controls verbosity.

## Value

A list of n_perm embedding matrices (each n_nodes x embd_dim).
