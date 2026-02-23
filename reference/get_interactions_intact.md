# Get the Intact network

Get the Intact network

## Usage

``` r
get_interactions_intact(
  threshold = NULL,
  physical_interactions = FALSE,
  head_only = FALSE
)
```

## Arguments

- threshold:

  Optional numeric. If you wish to filter on the STRING score.

- physical_interactions:

  Boolean. Shall only physical interactions be returned.

- head_only:

  Logical. If TRUE, returns only the first 5 rows.

## Value

data.table with the Intact network data.

## Details

The data is extracted from the OpenTargets release `25.12`
