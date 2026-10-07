# Distances between the columns of two matrices

**\[experimental\]** SIMD-accelerated distances between every column of
`mat_a` and every column of `mat_b`.

## Usage

``` r
rs_profile_distances(mat_a, mat_b, metric)
```

## Arguments

- mat_a:

  Numeric matrix. Columns are the samples.

- mat_b:

  Numeric matrix with the same number of rows as `mat_a`.

- metric:

  String. One of `c("correlation", "canberra", "l1", "l2", "cosine")`.

## Value

Numeric matrix of ncol(mat_a) x ncol(mat_b).
