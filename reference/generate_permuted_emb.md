# Generate the random, permuted embedding

This function will generate generate permuted subgraphs that show the
same node degree and run the node2vec algorithm on top of these
subgraphs. The cosine distance between the real edges will be calculated
and stored in the object. These will be used subsequently for
statistical significance testing. For more details, please refer to
Ietswaart, et al. The function will return a warning and the class as is
if you have not run
[`generate_initial_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_initial_emb.md).

## Usage

``` r
generate_permuted_emb(object, n_perm = 3L, seed = 42L, .verbose = TRUE)
```

## Arguments

- object:

  The `GeneWalk` class, please see
  [`GeneWalk()`](https://gregorlueg.github.io/genewalkR/reference/GeneWalk.md).

- n_perm:

  Inter. Number of permutations to generate. Defaults to `3L`, the
  recommendation of the authors.

- seed:

  Integer. Seed for reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function

## Value

The object with the added random embeddings for statistical testing.

## Details

The function will be using the same parameters for the node2vec training
as were used for
[`generate_initial_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_initial_emb.md).

## References

Ietswaart, et al., Genome Biol, 2021
