# Get node embeddings

Get node embeddings

## Usage

``` r
get_embeddings(object, which = NULL)
```

## Arguments

- object:

  An `EmbedDrift` object, see
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md).

- which:

  Integer. Which embedding to return: `1`, `2`, or `NULL` for both as a
  named list. Defaults to `NULL`.

## Value

A matrix, or a named list with elements `embd_1` and `embd_2`.
