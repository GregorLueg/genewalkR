# Wrapper function for the metapath2vec parameters

Parameters for
[`metapath2vec()`](https://gregorlueg.github.io/genewalkR/reference/metapath2vec.md).
There is no `p` or `q`: metapath walks are first-order and follow the
schema, not a biased return/explore rule.

## Usage

``` r
params_metapath2vec(
  walks_per_node = 40L,
  walk_length = 40L,
  n_epochs = 20L,
  n_negatives = 5L,
  window_size = 2L,
  lr = 0.01,
  sample = 0.001,
  num_workers = NULL
)
```

## Arguments

- walks_per_node:

  Integer. Number of random walks per node of the metapath's starting
  type. Defaults to `40L`.

- walk_length:

  Integer. Length of each random walk. Rounded up so the walk closes on
  a full number of schema cycles. Defaults to `40L`.

- n_epochs:

  Integer. Number of training epochs. Defaults to `20L`.

- n_negatives:

  Integer. Number of negative samples. Defaults to `5L`.

- window_size:

  Integer. Context window size. Defaults to `2L`.

- lr:

  Numeric. Learning rate. Defaults to `0.01`.

- sample:

  Numeric. Subsampling threshold for frequent nodes. Defaults to
  `0.001`.

- num_workers:

  Integer or `NULL`. Number of worker threads. If kept to `NULL`, it
  resolves to `available cores - 2 (min 1)`. Defaults to `NULL`.

## Value

A named list with the following elements:

- walks_per_node - Integer. Number of random walks per node of the
  metapath's starting type. Defaults to `40L`.

- walk_length - Integer. Length of each random walk. Rounded up so the
  walk closes on a full number of schema cycles. Defaults to `40L`.

- num_workers - Integer or `NULL`. Number of worker threads. If kept to
  `NULL`, it resolves to `available cores - 2 (min 1)`. Defaults to
  `NULL`.

- n_epochs - Integer. Number of training epochs. Defaults to `20L`.

- n_negatives - Integer. Number of negative samples. Defaults to `5L`.

- window_size - Integer. Context window size. Defaults to `2L`.

- lr - Numeric. Learning rate. Defaults to `0.01`.

- sample - Numeric. Subsampling threshold for frequent nodes. Defaults
  to `0.001`.
