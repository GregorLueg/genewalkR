# Calculate distances between diffusion profiles

Distances between two sets of stored profiles, e.g. drugs against
diseases. Lower means more similar. Ruiz et al. rank drug-disease pairs
by the correlation distance on the multiscale interactome and by
Canberra when comparing against gene expression signatures.

## Usage

``` r
calculate_profile_distances(
  object,
  from,
  to,
  metric = c("correlation", "canberra", "l1", "l2", "cosine")
)
```

## Arguments

- object:

  A `DiffusionProfiles` object with generated profiles.

- from:

  Character vector. Seeds for the rows of the result.

- to:

  Character vector. Seeds for the columns of the result.

- metric:

  Character. One of `"correlation"`, `"canberra"`, `"l1"`, `"l2"` or
  `"cosine"`. Defaults to `"correlation"`.

## Value

A numeric matrix of `length(from) x length(to)` distances.

## Details

Every profile sums to 1 and the seed keeps roughly `1 - alpha` of the
mass on itself. What each metric weights:

- `"correlation"`: Pearson on the full vector. Driven by the large
  entries, i.e. the immediate neighbourhood of the seed and hubs.

- `"cosine"`: Without centring. Since all profiles sum to 1 and are
  mostly near zero, it ranks almost identically to `"correlation"`.

- `"l1"`: Twice the total variation distance between the two
  distributions. Less dominated by the top entries than the above.

- `"l2"`: Dominated by the largest entries; favours big, well-connected
  nodes.

- `"canberra"`: Relative differences per node, so the many near-zero
  entries count as much as the large ones. Useful when the profile is
  compared against a signature on a different scale.

## References

Ruiz, Zitnik and Leskovec, Identification of disease treatment
mechanisms through the multiscale interactome, Nat Commun 2021.
