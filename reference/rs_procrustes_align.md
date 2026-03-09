# Orthogonal Procrustes alignment and cosine similarity scoring

Aligns `embd1` onto `embd2` via orthogonal Procrustes (SVD) and returns
cosine similarities per node. Both matrices must have identical
dimensions and matching row order, i.e. row i in `embd1` corresponds to
row i in `embd2`.

## Usage

``` r
rs_procrustes_align(embd1, embd2)
```

## Arguments

- embd1:

  Numeric matrix. k x d embedding matrix for graph 1.

- embd2:

  Numeric matrix. k x d embedding matrix for graph 2.

## Value

Named list with `aligned` (k x d matrix, `embd1` rotated onto `embd2`)
and `cosine_similarities` (numeric vector of length k).
