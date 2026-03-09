# Calculate embedding drift between two graphs

Aligns the embedding of graph 2 onto the reference embedding of graph 1
via orthogonal Procrustes, then computes per-node cosine similarities as
a measure of context drift. Nodes shared between both graphs receive a
cosine similarity in `[-1, 1]`. Nodes exclusive to graph 1 are assigned
`-1.1` and nodes exclusive to graph 2 are assigned `+1.1` as sentinel
values indicating absence rather than a computed similarity.

## Usage

``` r
calculate_drift(object)
```

## Arguments

- object:

  An `EmbedDrift` object, see
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md).

## Value

The `EmbedDrift` object with `aligned_embd` and `stats` populated.
`stats` is a `data.table` with columns `node`, `cosine_similarity`, and
`node_status` (`"shared"`, `"g1_only"`, `"g2_only"`).
