# Get the Gene Ontology

Get the Gene Ontology

## Usage

``` r
get_gene_ontology_hierarchy(head_only = FALSE)
```

## Arguments

- head_only:

  Logical. If TRUE, returns only the first 5 rows.

## Value

data.table with the GO to GO edges with edge type info.

## Details

The data is extracted from the obo file downloaded from [GO OBO
Foundry](https://obofoundry.org/ontology/go.html). Last download:
`January '26`.
