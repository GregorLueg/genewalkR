# Calculate GeneWalk test statistics

Pools all null distributions into one, then for each observed embedding
rep: computes p-values, per-gene FDR, and global FDR. Finally aggregates
across reps with geometric mean and 95\\

Matches the paper: "we repeat the above-described network representation
learning and significance testing procedures of the GWN nreps_graph
times and provide the mean and 95\\

## Usage

``` r
rs_gene_walk_test(
  gene_embds_list,
  pathway_embds_list,
  null_similarities,
  connected_pathways,
  verbose
)
```

## Arguments

- gene_embds_list:

  List of n_graph gene embedding matrices (n_genes x dim).

- pathway_embds_list:

  List of n_graph pathway embedding matrices (n_pathways x dim).

- null_similarities:

  List of n_perm numeric vectors (null cosine similarities to be
  pooled).

- connected_pathways:

  List. Per-gene integer vectors of connected pathway indices
  (1-indexed).

- verbose:

  Controls verbosity.

## Value

A list with per-pair statistics.
