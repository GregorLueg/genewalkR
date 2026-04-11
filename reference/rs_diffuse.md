# Run diffusion scoring on a precomputed kernel

Run diffusion scoring on a precomputed kernel

## Usage

``` r
rs_diffuse(
  kernel,
  scores,
  n_bkgd,
  n_inputs,
  bkgd_indices,
  method,
  n_perm,
  seed
)
```

## Arguments

- kernel:

  External pointer to `DiffusionKernel`

- scores:

  Numeric. The scores.

- n_bkgd:

  Integer. Number of background nodes.

- n_inputs:

  Integer. Number of input columns

- bkgd_indices:

  Integer. 0-based node indices for the background

- method:

  String. One of `c("raw", "z", "mc")`.

- n_perm:

  Integer. Number of permutations (MC only)

- seed:

  Integer RNG seed (MC only)

## Value

Scores per given node
