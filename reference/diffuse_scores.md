# Run diffusion scoring

Diffuses input scores over the precomputed kernel. Requires
[`build_kernel()`](https://gregorlueg.github.io/genewalkR/reference/build_kernel.md)
to have been called first.

## Usage

``` r
diffuse_scores(
  object,
  input,
  background = NULL,
  method = "raw",
  n_perm = 1000L,
  seed = 42L
)
```

## Arguments

- object:

  A `DiffusionScores` object with a built kernel.

- input:

  Named numeric vector or named matrix. Names must correspond to node
  names in the graph. For a vector, treated as a single input column.
  For a matrix, rows are background nodes and columns are separate input
  sets.

- background:

  Character vector. Node names forming the background set. Defaults to
  all nodes in `input`.

- method:

  Character. One of `"raw"`, `"z"`, or `"mc"`. Defaults to `"raw"`.

- n_perm:

  Integer. Number of permutations for `method = "mc"`. Defaults to
  `1000L`.

- seed:

  Integer. RNG seed for `method = "mc"`. Defaults to `42L`.

## Value

The `DiffusionScores` object with `scores` populated as a `data.table`
with columns `node` and one score column per input.
