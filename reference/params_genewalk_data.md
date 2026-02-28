# GeneWalk synthetic data params

GeneWalk synthetic data params

## Usage

``` r
params_genewalk_data(
  n_signal_genes = 150L,
  n_noise_genes = 150L,
  n_roots = 10L,
  depth = 3L,
  branching = 3L,
  p_lateral = 0.1,
  p_ppi = 0.05,
  min_annotations = 10L,
  max_annotations = 25L,
  min_noise_subtrees = 15L
)
```

## Arguments

- n_signal_genes:

  Integer. Number of signal genes. Defaults to `300L`.

- n_noise_genes:

  Integer. Number of noise genes. Defaults to `300L`.

- n_roots:

  Integer. Number of ontology root nodes. Defaults to `10L`.

- depth:

  Integer. Depth of each ontology subtree. Defaults to `3L`.

- branching:

  Integer. Base branching factor per ontology node. Defaults to `3L`.

- p_lateral:

  Numeric. Probability of a lateral edge between sibling nodes. Defaults
  to `0.1`.

- p_ppi:

  Numeric. Probability of a PPI edge between any two genes in the same
  group. Defaults to `0.05`.

- min_annotations:

  Integer. Minimum ontology term annotations per gene. Defaults to
  `10L`.

- max_annotations:

  Integer. Maximum ontology term annotations per gene. Defaults to
  `25L`.

- min_noise_subtrees:

  Integer. Minimum distinct ontology subtrees a noise gene must span.
  Defaults to `10L`.

## Value

A list with GeneWalk parameters.
