# Generate diffusion profiles

**\[experimental\]** Constrained personalised PageRank over a
heterogeneous graph (Ruiz et al., 2021). The graph is built once and one
profile per seed is computed in parallel. The seed acts as a source; all
other nodes of a sink type absorb mass.

## Usage

``` r
rs_diffusion_profiles(
  node_types,
  from,
  to,
  weights,
  type_weight_names,
  type_weight_values,
  sink_types,
  seeds,
  directed,
  diffusion_profile_params
)
```

## Arguments

- node_types:

  Character vector. Node type per node.

- from:

  Integer vector. 1-based node indices for edge origins.

- to:

  Integer vector. 1-based node indices for edge destinations.

- weights:

  Optional numeric vector. Edge weights, defaults to 1.

- type_weight_names:

  Optional character vector. Node types for `type_weight_values`. `NULL`
  gives the plain random walk.

- type_weight_values:

  Optional numeric vector. Weight per node type.

- sink_types:

  Character vector. Node types that act as sinks.

- seeds:

  List of integer vectors. 1-based node indices; one profile per
  element. A set restarts uniformly over its nodes.

- directed:

  Boolean. Treat the graph as directed.

- diffusion_profile_params:

  Named list with `alpha`, `max_iter` and `tol`.

## Value

Numeric matrix of n_nodes x n_seeds. Columns sum to 1.
