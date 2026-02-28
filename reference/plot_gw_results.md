# Plot the gene walk results

Plots the GeneWalk results. Creates an interactive plotly plot to
explore the results interactively.

## Usage

``` r
plot_gw_results(object, fdr_treshold = 0.25)
```

## Arguments

- object:

  The `GeneWalk` class. The stats need to be available in the object.

- fdr_treshold:

  Numeric. The FDR approach you wish to apply. GeneWalk is quite
  stringent, so, setting this to lower values can make sense.

## Value

ggplot2 object with the plot for the gene
