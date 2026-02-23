# Generate a new DataBuilder helper class

Generate a new DataBuilder helper class

## Usage

``` r
new_data_builder(ppis, gene_to_pathways, pathway_hierarchy)
```

## Arguments

- ppis:

  data.table. The full interaction data to subset subsequently.

- gene_to_pathways:

  data.table. The gene to pathway associations.

- pathway_hierarchy:

  data.table. The pathway ontology

## Value

Returns the initialised `DataBuilder` object.
