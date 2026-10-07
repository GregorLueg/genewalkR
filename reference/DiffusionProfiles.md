# DiffusionProfiles

S7 class for the diffusion profiles of Ruiz et al.: one constrained
personalised PageRank per seed over a heterogeneous graph. The seed is a
source; every other node of a sink type absorbs the walker. With type
weights, the walker at node `i` first picks a neighbouring node type
proportional to its weight, then a neighbour of that type.

## Usage

``` r
DiffusionProfiles(
  graph_dt,
  node_dt,
  sink_types,
  type_weights = NULL,
  directed = FALSE
)
```

## Arguments

- graph_dt:

  data.table. The edge table. Needs the columns `"from"` and `"to"`,
  optionally `"weight"` (non-negative). Every endpoint needs to be in
  `node_dt$id`.

- node_dt:

  data.table. The node table with the columns `"id"` and `"type"`.

- sink_types:

  Character vector. Node types that absorb the walker, e.g.
  `c("drug", "disease")`. Can be empty.

- type_weights:

  Optional named numeric vector. Positive weight per node type; needs to
  cover every type in `node_dt`. `NULL` gives the plain (edge-weighted)
  random walk.

- directed:

  Boolean. Treat the graph as directed. Defaults to `FALSE`.

## Value

An initialised `DiffusionProfiles` object.

## Properties

- graph_dt:

  data.table. The edge table with `from`, `to` and optionally `weight`.

- node_dt:

  data.table. The node table with `id` and `type`.

- sink_types:

  Character vector. Node types that act as sinks.

- type_weights:

  Named numeric vector of weights per node type, or `NULL` for the plain
  random walk.

- profiles:

  Numeric matrix of nodes x seeds. `NULL` until
  [`generate_profiles()`](https://gregorlueg.github.io/genewalkR/reference/generate_profiles.md)
  is called.

- params:

  Named list of the parameters used.

## References

Ruiz, Zitnik and Leskovec, Identification of disease treatment
mechanisms through the multiscale interactome, Nat Commun 2021.
