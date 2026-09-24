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

  Integer. Context window size. For GeneWalk you want embeddings from
  direct neighbours. Defaults to `1L`.

- lr:

  Numeric. Learning rate. Defaults to `0.01`.

- num_workers:

  Integer. Number of worker threads. Anything above `1L` makes the
  fitting much faster, at the cost of reproducibility. Defaults to `1L`.

## Value

A named list with the following elements:

- p - Numeric. Return parameter for biased random walks. Defaults to
  `1.0`.

- q - Numeric. In-out parameter for biased random walks. Defaults to
  `1.0`.

- walks_per_node - Integer. Number of random walks per node. Defaults to
  `100L`.

- walk_length - Integer. Length of each random walk. Defaults to `10L`.

- num_workers - Integer. Number of worker threads. Anything above `1L`
  makes the fitting much faster, at the cost of reproducibility.
  Defaults to `1L`.

- batch_size - Integer. Batch size for training. Defaults to `256L`.

- n_epochs - Integer. Number of training epochs. Defaults to `5L`.

- n_negatives - Integer. Number of negative samples. Defaults to `5L`.

- window_size - Integer. Context window size. For GeneWalk you want
  embeddings from direct neighbours. Defaults to `1L`.

- lr - Numeric. Learning rate. Defaults to `0.01`.
