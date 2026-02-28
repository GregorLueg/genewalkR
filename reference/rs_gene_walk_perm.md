# Generate null distributions from degree-matched random networks

For each permutation, generates a random graph via the configuration
model (matching the original degree distribution), trains node2vec on
it, then collects cosine similarities between each node and its unique
neighbours in the random graph. This matches the original Python
GeneWalk procedure.

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

  Integer. Number of null permutations.

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed.

- verbose:

  Boolean. Controls verbosity.

## Value

A list of n_perm numeric vectors of null cosine similarities.
