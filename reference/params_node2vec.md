# Wrapper function for the node2vec parameters

These are general node2vec parameters, designed for speed. They will
default to `available cores - 2 (min 1)`.

## Usage

``` r
params_node2vec(
  p = 1,
  q = 1,
  walks_per_node = 40L,
  walk_length = 40L,
  batch_size = 256L,
  n_epochs = 20L,
  n_negatives = 5L,
  window_size = 2L,
  lr = 0.01,
  num_workers = NULL
)
```

## Arguments

- p:

  Numeric. Return parameter for biased random walks. Defaults to `1.0`.

- q:

  Numeric. In-out parameter for biased random walks. Defaults to `1.0`.

- walks_per_node:

  Integer. Number of random walks per node. Defaults to `40L`.

- walk_length:

  Integer. Length of each random walk. Defaults to `40L`.

- batch_size:

  Integer. Batch size for training. Defaults to `256L`.

- n_epochs:

  Integer. Number of training epochs. Defaults to `20L`.

- n_negatives:

  Integer. Number of negative samples. Defaults to `5L`.

- window_size:

  Integer. Context window size. Defaults to `2L`.

- lr:

  Numeric. Learning rate. Defaults to `1e-2`.

- num_workers:

  Optional integer. If kept to `NULL`, it will default to
  `available cores - 2 (min 1)`.

## Value

A list with the node2vec parameters.
