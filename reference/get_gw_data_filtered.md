# Get the filtered GW data for a list of genes

Get the filtered GW data for a list of genes

## Usage

``` r
get_gw_data_filtered(x, gene_ids)

# S3 method for class 'DataBuilder'
get_gw_data_filtered(x, gene_ids)
```

## Arguments

- x:

  An object from which to extract the genewalk data

- gene_ids:

  String. The gene ids to which to filter the data.

## Value

A list with the following elements

- gwn - The final GeneWalk network.

- genes_to_pathways - The genes to pathway information.

- represented_genes - All included genes.

- represented_pathways - All included pathways.
