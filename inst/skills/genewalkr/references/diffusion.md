# Diffusion

Two separate tools. `DiffusionScores` spreads input scores over a homogeneous
network through a graph kernel (guilt by association, network propagation).
`DiffusionProfiles` runs one constrained personalised PageRank per seed over a
heterogeneous, typed graph (Ruiz et al.) and compares the resulting profiles.

## DiffusionScores

```r
# graph_from_data_frame() names the vertices
g <- igraph::graph_from_data_frame(edges, directed = FALSE)

obj <- DiffusionScores(g)
obj <- build_kernel(
  obj,
  kernel = "regularised_laplacian",
  kernel_params = params_kernel(sigma2 = 1, add_diag = 1),
  normalised = TRUE,
  strategy = "full"
)
obj <- diffuse_scores(obj, input = seeds, method = "z")
get_scores(obj)
```

### The kernel

`build_kernel()` builds the Laplacian, eigendecomposes it and transforms the
eigenvalues. Do it once per graph and kernel, then call `diffuse_scores()` as
often as you like.

- Kernels: `"regularised_laplacian"` (default), `"diffusion"`, `"pstep"`,
  `"commute_time"`, `"inverse_cosine"`. Which `params_kernel()` elements each
  one reads is in `references/conventions.md`.
- `strategy = "truncated"` (default) keeps the top `k = 100L` eigenpairs;
  `"full"` does the complete decomposition, which is cubic in the node count.
  The vignettes use `"full"` on graphs of a few hundred nodes.
- The kernel lives behind an external pointer. It doesn't survive `saveRDS()`
  or a restart; rebuild it.
- The graph needs vertex names (`igraph::V(g)$name`). Inputs, background and
  outputs are all keyed by them.
- The adjacency is unweighted: an igraph `weight` edge attribute is ignored.

### Scoring

`input` is a named numeric vector, or a matrix with node names as rownames and
one column per input set. Names have to be a subset of `background`, which
defaults to every node. Nodes in the background but missing from `input` score
0.

| `method` | What you get |
|---|---|
| `"raw"` | kernel times input |
| `"z"` | raw score standardised with the analytical mean and variance under permuted inputs |
| `"mc"` | heatrank from `n_perm` permutations (`seed` applies): `(r + 1) / (n_perm + 1)` in `(0, 1)`, lower = warmer |

Raw scores have no null. `"z"` and `"mc"` score each node against permuted
inputs, which is what you want when ranking candidates. Mind the direction:
higher z is warmer, lower heatrank is warmer.

`get_scores()` returns a `data.table` with `node` and one column per input,
named after the matrix colnames or `input_1`, `input_2`... for a vector. Each
call replaces the stored scores.

## DiffusionProfiles

Ruiz, Zitnik and Leskovec, Nat Commun 2021. For each seed (a drug, a disease,
a gene signature node), a personalised PageRank whose walker is absorbed by
every other node of a sink type, so it can't shortcut through another drug or
disease. With `type_weights`, the walker first picks a neighbouring node type
in proportion to its weight, then a neighbour of that type.

```r
edges <- dedup_edges(edges, node_dt = nodes)

dp <- DiffusionProfiles(
  graph_dt = edges,                  # from, to, optional weight
  node_dt = nodes,                   # id, type
  sink_types = c("drug", "disease"),
  type_weights = c(drug = 1, protein = 1, `function` = 5, disease = 1)
)
dp <- generate_profiles(dp, seeds = c(drugs, diseases))
calculate_profile_distances(dp, from = drugs, to = diseases)
```

- `sink_types` can be `character()`: then it's a plain personalised PageRank.
- `type_weights = NULL` gives the plain edge-weighted walk. If given, it needs
  a positive weight for every type in `node_dt`.
- **Deduplicate first.** Parallel edges add up in the transition matrix, and
  `A-B` plus `B-A` counts twice in an undirected graph. The constructor warns;
  `dedup_edges()` fixes it (`weight_agg = "max"` by default).
- **Each `generate_profiles()` call replaces the stored profiles.** Pass every
  seed you want to compare in one call.
- Profiles are a nodes x seeds matrix; every column sums to 1, and roughly
  `1 - alpha` of the mass stays on the seed.

### Distances

`calculate_profile_distances(dp, from, to, metric)` returns a
`length(from) x length(to)` matrix, lower = closer.

| `metric` | Weights |
|---|---|
| `"correlation"` (default) | large entries: the seed's neighbourhood and hubs. Ruiz et al. use it for drug-disease ranking |
| `"cosine"` | ranks almost identically to correlation on these profiles |
| `"l1"` | twice the total variation distance; less top-heavy |
| `"l2"` | the largest entries; favours big, well-connected nodes |
| `"canberra"` | relative differences per node, so near-zero entries count. Ruiz et al. use it against expression signatures |

The diffusion profile vignette runs MYC targets against all of Reactome and
compares metrics against an over-representation reference. Worth reading
before choosing one.
