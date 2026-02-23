# Calculate the pathway and gene statistics

Calculates statistical significance for gene-pathway associations using
the GeneWalk algorithm. For each gene-pathway pair, this function
computes p-values and FDR corrections (both globally across all pairs
and per-gene) across multiple network permutations. The final output
includes mean values and 95% confidence intervals for all statistics.

The function processes all permutations in Rust for optimal performance,
calculating cosine similarities between gene and pathway embeddings,
determining p-values based on null distributions from randomised
networks, and applying Benjamini-Hochberg FDR corrections.

## Usage

``` r
calculate_genewalk_stats(object, .verbose = TRUE)
```

## Arguments

- object:

  The `GeneWalk` object containing actual and permuted embeddings. See
  [`GeneWalk()`](https://gregorlueg.github.io/genewalkR/reference/GeneWalk.md).

- .verbose:

  Logical. If `TRUE`, prints progress messages during calculation.
  Default is `TRUE`.

## Value

The input object with an added `stats` property containing a data.table
with the following columns for each gene-pathway pair:

- gene - Gene identifier

- pathway - Pathway identifier

- avg_pval - Mean p-value across permutations

- pval_ci_lower, pval_ci_upper - 95% CI for p-values

- avg_global_fdr - Mean global FDR across permutations

- global_fdr_ci_lower,global_fdr_ci_upper - 95% CI for global FDR

- avg_gene_fdr - Mean gene-wise FDR across permutations

- gene_fdr_ci_lower,gene_fdr_ci_upper - 95% CI for gene-wise FDR

## References

Ietswaart, et al., Genome Biol, 2021
