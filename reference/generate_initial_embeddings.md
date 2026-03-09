# Generate initial node embeddings for both graphs

Runs node2vec independently on both graphs using identical
hyperparameters, storing the resulting embedding matrices in the object.
Graph 1 serves as the reference for subsequent Procrustes alignment. For
details on the node2vec algorithm, see Grover & Leskovec. For the
approach generally speaking, see Amir, et al.

## Usage

``` r
generate_initial_embeddings(
  object,
  embd_dim = 64L,
  node2vec_params = params_node2vec(),
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  An `EmbedDrift` object, see
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md).

- embd_dim:

  Integer. Embedding dimension size. Defaults to `64L`.

- node2vec_params:

  Named list. Node2vec training parameters, see
  [`params_node2vec()`](https://gregorlueg.github.io/genewalkR/reference/params_node2vec.md).

- directed:

  Boolean. Whether to treat both graphs as directed. Defaults to
  `FALSE`.

- seed:

  Integer. Random seed for reproducibility.

- .verbose:

  Boolean. Controls verbosity. Defaults to `TRUE`.

## Value

The `EmbedDrift` object with `embd_1`, `embd_2`, and `params` populated.

## References

Grover & Leskovec, arXiv, 2016  
Amir, et al., Cancer Cell, 2025
