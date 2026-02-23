# Node2Vec

Runs the Node2Vec algorithm

## Usage

``` r
node2vec(
  graph_dt,
  embd_dim = 8L,
  node2vec_params = params_node2vec(),
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- graph_dt:

  data.table. The edge data.frame. Needs to have the columns `"from"`
  and `"to"`. Can also optionally a `"weight"` column.

- embd_dim:

  Integer. Size of the embedding dimensions to create. Defaults to `8L`
  in line with the authors recommendations.

- node2vec_params:

  Named list. Contains the node2vec parameters. The list has the
  following elements:

  - p - Numeric. Return parameter for biased random walks. Defaults to
    `1.0`.

  - q - Numeric. In-out parameter for biased random walks. Defaults to
    `1.0.`.

  - walks_per_node - Integer. Number of random walks per node. Defaults
    to ` 40L.`.

  - walk_length - Integer. Length of each random walk. Defaults to
    `40L`.

  - num_workers - Number of worker threads during batching. Defaults to
    `4L`. To note: during graph generation, Rust will use via the Rayon
    backend all available threads.

  - batch_size - Integer. Batch size for training. Defaults to `256L`.

  - n_epochs - Integer. Number of training epochs. Defaults to `15L`.

  - n_negatives - Integer. Number of negative samples. Defaults to `5L`.

  - window_size - Integer. Context window size. Defaults to `2L`.

  - lr - Numeric. Learning rate. Defaults to `1e-2`.

- directed:

  Boolean. Indicates if this is a directed or undirected network.
  Defaults to `FALSE`.

- seed:

  Integer. Seed for reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function
