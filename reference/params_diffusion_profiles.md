# Wrapper function for the diffusion profile parameters

Parameters for
[`generate_profiles()`](https://gregorlueg.github.io/genewalkR/reference/generate_profiles.md),
i.e. the constrained personalised PageRank behind the diffusion profiles
of Ruiz et al.

## Usage

``` r
params_diffusion_profiles(alpha = 0.85, max_iter = 100L, tol = 1e-06)
```

## Arguments

- alpha:

  Numeric. Probability of continuing the walk instead of restarting at
  the seed. Defaults to `0.85`.

- max_iter:

  Integer. Maximum number of power iterations. Defaults to `100L`.

- tol:

  Numeric. Convergence threshold on the L1 change between two
  iterations. Defaults to `1e-06`.

## Value

A named list with the following elements:

- alpha - Numeric. Probability of continuing the walk instead of
  restarting at the seed. Defaults to `0.85`.

- max_iter - Integer. Maximum number of power iterations. Defaults to
  `100L`.

- tol - Numeric. Convergence threshold on the L1 change between two
  iterations. Defaults to `1e-06`.
