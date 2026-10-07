# Generate diffusion profiles

Runs one constrained personalised PageRank per seed (Ruiz et al.). The
graph is built once in Rust and the seeds run in parallel. Each call
replaces the stored profiles, so pass drugs and diseases together if you
want to compare them afterwards.

## Usage

``` r
generate_profiles(
  object,
  seeds,
  diffusion_profile_params = params_diffusion_profiles(),
  .verbose = TRUE
)
```

## Arguments

- object:

  A `DiffusionProfiles` object, see
  [`DiffusionProfiles()`](https://gregorlueg.github.io/genewalkR/reference/DiffusionProfiles.md).

- seeds:

  Character vector. Node ids to start the walks from; one profile per
  seed.

- diffusion_profile_params:

  Named list. The PageRank parameters, see
  [`params_diffusion_profiles()`](https://gregorlueg.github.io/genewalkR/reference/params_diffusion_profiles.md).

- .verbose:

  Boolean. Controls verbosity. Defaults to `TRUE`.

## Value

The `DiffusionProfiles` object with `profiles` populated as a numeric
matrix of nodes x seeds. Each column sums to 1.
