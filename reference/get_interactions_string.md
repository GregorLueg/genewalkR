# Get the STRING network

Get the STRING network

## Usage

``` r
get_interactions_string(threshold = NULL, head_only = FALSE)
```

## Arguments

- threshold:

  Optional numeric. If you wish to filter on the STRING score.

- head_only:

  Logical. If TRUE, returns only the first 5 rows.

## Value

data.table with the String network data.

## Details

The data is extracted from the OpenTargets release `25.12`
