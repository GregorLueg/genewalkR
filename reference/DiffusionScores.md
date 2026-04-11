# DiffusionScores

S7 class wrapping the Rust spectral diffusion backend. Holds the input
graph, the precomputed kernel (as an external pointer), and the
resulting per-node diffusion scores.

## Usage

``` r
DiffusionScores(graph)
```

## Arguments

- graph:

  An igraph object.

## Value

An initialised `DiffusionScores` object.

## Properties

- graph:

  An igraph object. The input graph.

- kernel_pointer:

  External pointer to the Rust `DiffusionKernel`. `NULL` until
  [`build_kernel()`](https://gregorlueg.github.io/genewalkR/reference/build_kernel.md)
  is called.

- scores:

  A data.table of diffusion results. Empty until `diffuse()` is called.

- params:

  Named list of parameters used during kernel construction and
  diffusion.
