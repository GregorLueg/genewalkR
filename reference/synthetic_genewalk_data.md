# Generate synthetic GeneWalk data

Generates a synthetic network that mimics real GeneWalk input data for
testing and benchmarking. The network consists of three components: a
hierarchical ontology (tree structure with lateral edges),
gene-to-ontology annotations, and protein-protein interactions (PPI).

Signal genes are annotated exclusively within a single ontology subtree
and form PPI edges among themselves, creating coherent graph
neighbourhoods. Noise genes receive the same number of annotations but
are forced to span multiple subtrees, and their PPI edges are random
across the noise group. Both groups are degree-matched so that
significance is driven by annotation coherence rather than connectivity
alone.

## Usage

``` r
synthetic_genewalk_data(data_params = params_genewalk_data(), seed = 42L)
```

## Arguments

- data_params:

  Named list. Return from
  [`params_genewalk_data()`](https://gregorlueg.github.io/genewalkR/reference/params_genewalk_data.md).

- seed:

  Integer. Random seed. Defaults to `42L`.

## Value

A list with the following items:

- full_data:

  A data.table with columns `from` and `to` containing all edges
  (ontology hierarchy, gene-to-pathway, and PPI).

- gene_to_pathways:

  A data.table with columns `from` (gene) and `to` (pathway) for
  gene-to-ontology annotation edges only.

- gene_ids:

  Character vector of all gene identifiers (signal and noise).

- pathway_ids:

  Character vector of all ontology term identifiers.
