# Build the diffusion kernel

Constructs the graph Laplacian and runs the spectral decomposition, then
applies the chosen kernel transform. The result is stored as an external
pointer for subsequent calls to `diffuse()`.

## Usage

``` r
build_kernel(
  object,
  kernel = c("regularised_laplacian", "commute_time", "diffusion", "inverse_cosine",
    "pstep"),
  kernel_params = params_kernel(),
  normalised = TRUE,
  strategy = "truncated",
  k = 100L,
  .verbose = TRUE
)
```

## Arguments

- object:

  A `DiffusionScores` object, see
  [`DiffusionScores()`](https://gregorlueg.github.io/genewalkR/reference/DiffusionScores.md).

- kernel:

  Character. One of `"regularised_laplacian"`, `"commute_time"`,
  `"diffusion"`, `"inverse_cosine"`, `"pstep"`. Defaults to
  `"regularised_laplacian"`.

- kernel_params:

  Named list. Kernel-specific parameters, see
  [`params_kernel()`](https://gregorlueg.github.io/genewalkR/reference/params_kernel.md).

- normalised:

  Logical. Use the normalised Laplacian. Defaults to `TRUE`.

- strategy:

  Character. `"full"` or `"truncated"`. Defaults to `"truncated"`.

- k:

  Integer. Number of eigenvalues for the truncated strategy. Ignored if
  `strategy = "full"`. Defaults to `100L`.

- .verbose:

  Logical. Controls verbosity. Defaults to `TRUE`.

## Value

The `DiffusionScores` object with `kernel_pointer` and `params`
populated.
