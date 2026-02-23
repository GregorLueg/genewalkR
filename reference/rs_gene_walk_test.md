# Calculate the test statistics

Calculates test statistics for gene-pathway pairs. The null distribution
is derived from gene-pathway cosine similarities in permuted embeddings
(full cross-cosine, not just connected pairs).

## Usage

``` r
rs_gene_walk_test(
  gene_embds,
  pathway_embds,
  permuted_embds,
  gene_indices,
  pathway_indices,
  connected_pathways,
  verbose
)
```

## Arguments

- gene_embds:

  Matrix of n_genes x embedding dimensions.

- pathway_embds:

  Matrix of n_pathways x embedding dimensions.

- permuted_embds:

  List of permuted embedding matrices (n_nodes x dim).

- gene_indices:

  Integer vector. 1-based row indices for genes in permuted embeddings.

- pathway_indices:

  Integer vector. 1-based row indices for pathways in permuted
  embeddings.

- connected_pathways:

  List. Gene to pathway connections (1-indexed).

- verbose:

  Controls verbosity.

## Value

A list with per-pair statistics (see original docs).
