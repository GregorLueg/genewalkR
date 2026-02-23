# Generate synthetic GeneWalk data with controlled signal structure

Generate synthetic GeneWalk data with controlled signal structure

## Usage

``` r
rs_build_synthetic_genewalk(
  n_signal_genes,
  n_noise_genes,
  n_roots,
  depth,
  branching,
  p_lateral,
  p_ppi,
  min_annotations,
  max_annotations,
  min_noise_subtrees,
  seed
)
```

## Arguments

- n_signal_genes:

  Integer. Genes annotated to a single ontology subtree.

- n_noise_genes:

  Integer. Genes with annotations scattered across subtrees.

- n_roots:

  Integer. Number of ontology root terms.

- depth:

  Integer. Depth of each ontology subtree.

- branching:

  Integer. Average branching factor per node.

- p_lateral:

  Numeric. Probability of lateral edges within each ontology level.

- p_ppi:

  Numeric. PPI connection probability within gene groups.

- min_annotations:

  Integer. Minimum annotations per gene.

- max_annotations:

  Integer. Maximum annotations per gene.

- min_noise_subtrees:

  Integer. Minimum number of different subtrees each noise gene must
  span.

- seed:

  Integer. Random seed.

## Value

A list with: ontology_edges, gene_ont_edges, ppi_edges, signal_genes,
noise_genes.
