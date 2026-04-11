# Build the diffusion kernel from sparse adjacency components

Takes the CSC slots directly from R's dgCMatrix (igraph output).

## Usage

``` r
rs_diffusion_kernel(
  i,
  p,
  x,
  n,
  kernel,
  kernel_params,
  normalised,
  strategy,
  k,
  node_names,
  verbose
)
```

## Arguments

- i:

  Integers. Row indices (0-based, from dgCMatrix@i)

- p:

  Integers. Column pointers (from dgCMatrix@p)

- x:

  Numeric. Values (from dgCMatrix@x)

- n:

  Integer. Matrix dimension

- kernel:

  String. Kernel type string. One of
  `c("regularised_laplacian", "commute_time", "inverse_cosine", "pstep")`.

- kernel_params:

  Named list. Contains the kernel-specific parameters.

- normalised:

  Boolean. Use normalised Laplacian

- strategy:

  String. One of `"full"` or `"truncated"`. Shall the full
  Eigendecomposition of a truncated Eigendecomposition be applied on the
  data. The latter is useful for large graphs to reduce memory pressure.

- k:

  Integer. Number of eigenvalues for truncated

- node_names:

  String. Character vector of node names.

- verbose:

  Boolean. Verbosity of the function.

## Value

External pointer to `DiffusionKernel`
