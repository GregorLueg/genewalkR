# Get the Combined network

This returns a combined network akin to the workflow from
Barrio-Hernandez et al. based on high quality edges from STRING,
physical interactions from Intact and the full sets from SIGNOR,
Reactome and Pathway Commons interactions.

## Usage

``` r
get_interactions_combined(head_only = FALSE)
```

## Arguments

- head_only:

  Logical. If TRUE, returns only the first 5 rows.

## Value

data.table with the combined, collapsed network data.

## Details

The data is extracted from the OpenTargets release `25.12`

## References

Barrio-Hernandez et al., Nat Gen, 2023
