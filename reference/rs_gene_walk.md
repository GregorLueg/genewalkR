# Generate GeneWalk node embeddings (multiple reps)

Trains node2vec on the original GeneWalk network `n_graph` times with
different random seeds, returning a list of embedding matrices. This
provides the variance across embedding reps needed for the final
log_stats aggregation.

## Usage

``` r
rs_gene_walk(
  from,
  to,
  weights,
  gene_walk_params,
  n_graph,
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

- n_graph:

  Integer. Number of independent embedding reps.

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed (incremented per rep).

- verbose:

  Boolean. Controls verbosity.

## Value

A list of n_graph numeric matrices, each n_nodes x embedding dim.
