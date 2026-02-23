# Generate pathway structure and gene-pathway associations

Generate pathway structure and gene-pathway associations

## Usage

``` r
rs_generate_pathway_structure(
  n_pathways,
  pathway_depth,
  pathway_branching,
  n_communities,
  gene_ids,
  gene_communities,
  n_focal_pathways,
  signal_strength,
  connections_per_gene,
  seed
)
```

## Arguments

- n_pathways:

  Integer. Number of pathways to generate.

- pathway_depth:

  Integer. Maximum depth of pathway hierarchy.

- pathway_branching:

  Integer. Average branching factor for pathway tree.

- n_communities:

  Integer. Number of gene communities.

- gene_ids:

  Character vector. Gene identifiers.

- gene_communities:

  Integer vector. Community assignment for each gene.

- n_focal_pathways:

  Integer. Number of focal pathways per community.

- signal_strength:

  Numeric. Probability of connecting to focal pathways (0-1).

- connections_per_gene:

  Integer. Number of pathway connections per gene.

- seed:

  Integer. Random seed.

## Value

A list with pathway edges, metadata, gene-pathway associations ano
community to subtree parts.
