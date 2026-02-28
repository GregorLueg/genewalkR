# Wrapper function for the GeneWalk-specific node2vec parameters

Wrapper function for the GeneWalk-specific node2vec parameters

## Usage

``` r
params_genewalk(
  p = 1,
  q = 1,
  walks_per_node = 100L,
  walk_length = 10L,
  batch_size = 256L,
  n_epochs = 5L,
  n_negatives = 5L,
  window_size = 1L,
  lr = 0.01,
  num_workers = 1L
)
```

## Arguments

- p:

  Numeric. Return parameter for biased random walks. Defaults to `1.0`.

- q:

  Numeric. In-out parameter for biased random walks. Defaults to `1.0`.

- walks_per_node:

  Integer. Number of random walks per node. Defaults to `100L`.

- walk_length:

  Integer. Length of each random walk. Defaults to `10L`.

- batch_size:

  Integer. Batch size for training. Defaults to `256L`.

- n_epochs:

  Integer. Number of training epochs. Defaults to `5L`.

- n_negatives:

  Integer. Number of negative samples. Defaults to `5L`.

- window_size:

  Integer. Context window size. Defaults to `1L`. For GeneWalk you want
  embedding from direct neighbours.

- lr:

  Numeric. Learning rate. Defaults to `1e-2`.

- num_workers:

  Integer. Number of worker threads. Defaults to `1L` (for determinism.
  If set to ≥1, the fitting is much faster at cost of deterministic
  results.

## Value

A list with the (gene-walk) specific node2vec parameters.
