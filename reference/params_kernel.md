# Wrapper function for the diffusion kernel parameters

These are parameters for the graph diffusion kernels. Not all parameters
are used by all kernels, please see Details.

## Usage

``` r
params_kernel(sigma2 = 1, add_diag = 1, a = 3, p = 5L)
```

## Arguments

- sigma2:

  Numeric. Bandwidth parameter controlling the extent of diffusion. Used
  by `regularised_laplacian` and `diffusion` kernels. Defaults to `1.0`.

- add_diag:

  Numeric. Regularisation term added to the diagonal. Used by
  `regularised_laplacian`. Must be positive for every connected
  component to be invertible. Defaults to `1.0`.

- a:

  Numeric. Regularisation parameter for the `pstep` kernel. Must be
  `>= 2` (eigenvalues of the normalised Laplacian are in `[0, 2]`.
  Defaults to `3.0`.

- p:

  Integer. Number of steps for the `pstep` kernel. Must be a positive
  integer. Defaults to `5L`.

## Value

A named list.

## Details

Parameter usage by kernel type:

- `regularised_laplacian`: `sigma2`, `add_diag`

- `commute_time`: none

- `diffusion`: `sigma2`

- `inverse_cosine`: none

- `pstep`: `a`, `p`
