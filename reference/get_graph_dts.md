# Get graph edge tables

Get graph edge tables

## Usage

``` r
get_graph_dts(object, which = NULL)
```

## Arguments

- object:

  An `EmbedDrift` object, see
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md).

- which:

  Integer. Which graph table to return: `1`, `2`, or `NULL` for both as
  a named list. Defaults to `NULL`.

## Value

A data.frame, or a named list with elements `graph_dt_1` and
`graph_dt_2`.
