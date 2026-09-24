# Generate metapath2vec embeddings

Runs metapath-constrained random walks over a heterogeneous graph and
trains the skip-gram model on them (metapath2vec, or metapath2vec++ with
per-type negative sampling).

## Usage

``` r
rs_metapath2vec(
  node_ids,
  node_types,
  from,
  to,
  weights,
  metapath,
  metapath_plus,
  metapath2vec_params,
  embd_dim,
  directed,
  seed,
  verbose
)
```

## Arguments

- node_ids:

  Character vector. Node identifiers.

- node_types:

  Character vector. Node type per node, same length as `node_ids`.

- from:

  Integer vector. 1-based indices into `node_ids` for edge origins.

- to:

  Integer vector. 1-based indices into `node_ids` for edge destinations.

- weights:

  Optional numeric vector. Edge weights.

- metapath:

  String. Hyphen-separated metapath closing on its starting type, e.g.
  `"gene-pathway-gene"`.

- metapath_plus:

  Boolean. Per-type negative sampling (metapath2vec++).

- metapath2vec_params:

  Named list. Training parameters (walks_per_node, walk_length,
  num_workers, n_epochs, n_negatives, window_size, lr, sample).

- embd_dim:

  Integer. Embedding dimension.

- directed:

  Boolean. Treat graph as directed.

- seed:

  Integer. Random seed.

- verbose:

  Boolean. Controls verbosity.

## Value

A list with:

- embedding - Matrix of n_nodes x embd_dim, rows in `node_names` order.

- node_names - Node identifiers in row order.

- node_types - Node type per row.

- visited - Logical per row. `FALSE` if no surviving walk touched the
  node, i.e. its row is the random initialisation.

- walk_stats - Named list with the walk generation statistics.
