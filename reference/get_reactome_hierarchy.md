# Get the Reactome pathway hierarchy

Get the Reactome pathway hierarchy

## Usage

``` r
get_reactome_hierarchy(
  relationship = c("child_of", "descendant_of", "parent_of", "ancestor_of"),
  head_only = FALSE
)
```

## Arguments

- relationship:

  String vector. Defines which relationships to return. Multiple allowed
  and choices from
  `c("child_of", "descendant_of", "parent_of", "ancestor_of")`

- head_only:

  Logical. If TRUE, returns only the first 5 rows.

## Value

data.table with the Reactome to Reactome edges with edge type info.

## Details

The data is extracted from the OpenTargets release `25.12`
