# Get the aligned embedding

Get the aligned embedding

## Usage

``` r
get_aligned_embd(object)
```

## Arguments

- object:

  An `EmbedDrift` object, see
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md).

## Value

A numeric matrix of the Procrustes-aligned embedding of graph 1 onto
graph 2, restricted to shared nodes.
