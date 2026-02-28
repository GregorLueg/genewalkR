# Plot similarity distributions

Plots the cosine similarity distributions of the permuted embeddings
alongside the actual values as density histograms.

## Usage

``` r
plot_similarities(object)
```

## Arguments

- object:

  The `GeneWalk` class. Permutations and stats need to be available in
  the object.

## Value

A ggplot2 object with the similarity distributions.
