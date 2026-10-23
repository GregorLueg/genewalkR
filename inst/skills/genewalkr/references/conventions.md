# genewalkR conventions

## Parameter bundles

Six of them. Each returns a plain named list checked by a matching assertion, so
a misspelt or mistyped element fails at the call site.

| Bundle | Feeds | Defaults worth knowing |
|---|---|---|
| `params_genewalk()` | `generate_initial_emb()` | `walks_per_node = 100L`, `walk_length = 10L`, `window_size = 1L`, `n_epochs = 5L`, `num_workers = 1L` |
| `params_node2vec()` | `node2vec()`, `generate_initial_embeddings()` | `walks_per_node = 40L`, `walk_length = 40L`, `window_size = 2L`, `n_epochs = 20L`, `num_workers = NULL` (cores minus two) |
| `params_metapath2vec()` | `metapath2vec()` | as node2vec, no `p`/`q`, plus `sample = 1e-3` |
| `params_kernel()` | `build_kernel()` | `sigma2 = 1`, `add_diag = 1`, `a = 3`, `p = 5L` |
| `params_diffusion_profiles()` | `generate_profiles()` | `alpha = 0.85`, `max_iter = 100L`, `tol = 1e-6` |
| `params_genewalk_data()` | `synthetic_genewalk_data()` | 150 signal and 150 noise genes |

`params_genewalk()` and `params_node2vec()` produce the same list shape and pass
the same assertion, so either is accepted wherever the other is. The GeneWalk
defaults (short walks, `window_size = 1L`) are deliberate: GeneWalk wants
direct neighbour context. Don't swap in `params_node2vec()` for GeneWalk
without a reason.

Not every kernel uses every `params_kernel()` element: `regularised_laplacian`
uses `sigma2` and `add_diag`, `diffusion` uses `sigma2`, `pstep` uses `a`
(must be `>= 2`) and `p`, `commute_time` and `inverse_cosine` use none.

## The S7 pipelines

| Class | Chain | Results |
|---|---|---|
| `GeneWalk` | `GeneWalk()` -> `generate_initial_emb()` -> `generate_permuted_emb()` -> `calculate_genewalk_stats()` | `get_stats()`, `get_embedding()`, `plot_similarities()`, `plot_gw_results()` |
| `EmbedDrift` | `EmbedDrift()` -> `generate_initial_embeddings()` -> `calculate_drift()` | `get_stats()`, `get_embeddings()`, `get_aligned_embd()` |
| `DiffusionScores` | `DiffusionScores()` -> `build_kernel()` -> `diffuse_scores()` | `get_scores()` |
| `DiffusionProfiles` | `DiffusionProfiles()` -> `generate_profiles()` | `get_profiles()`, `calculate_profile_distances()` |

Mind the near-identical names: `generate_initial_emb()` is GeneWalk,
`generate_initial_embeddings()` is EmbedDrift. `get_embedding()` is GeneWalk,
`get_embeddings()` is EmbedDrift.

Each verb returns a new object; reassign. Running a stage out of order:

- `generate_permuted_emb()` without an initial embedding warns and returns the
  object unchanged. Same for `calculate_genewalk_stats()` without permutations
  and `calculate_drift()` without embeddings.
- `diffuse_scores()` without a kernel errors. So does
  `calculate_profile_distances()` without profiles.
- Getters on an unpopulated slot warn and return an empty `data.table`, an empty
  list or `NULL`.

A warning here means the stage didn't run. Don't carry on as if it did.

`node2vec()` and `metapath2vec()` are plain functions: edge table in, matrix
out, no class.

`GeneWalkGenerator` is R6, not S7. Its methods mutate in place and chain with
`$`: `gen$add_pathways(c("go", "reactome"))$add_ppi("combined")$build()`.

## Edge and node tables

Everything takes a `data.table` (not a plain `data.frame`, except
`EmbedDrift()`, which accepts either).

- Edges: columns `from` and `to`, optionally a non-negative `weight`. Use
  character IDs.
- `GeneWalk` additionally needs `type` on every edge: `"interaction"` (gene to
  gene), `"part_of"` (gene to pathway) or `"hierarchy"` (pathway to pathway),
  and all three have to be present. The constructor errors otherwise.
- `metapath2vec()` and `DiffusionProfiles()` need a node table with `id`
  (unique) and `type`. Every edge endpoint has to be in `id`.
- Graphs are undirected unless `directed = TRUE`.

Duplicate edges are not dropped for you in `DiffusionProfiles()`; it warns and
each copy adds to the transition weight. Run `dedup_edges()` first.

## Embedding output

Embedding matrices have node IDs as rownames and `emb_dim_1`... as colnames
(`dim_1`... for `EmbedDrift`). Index by rowname, never by position: the row
order is not the order of your input.

## Verbosity and seeds

- `.verbose` (leading dot) is a single boolean.
- `seed` is an integer, default `42L`. With `num_workers > 1L` the seed does not
  make training reproducible.

## The `rs_*` functions

Raw Rust bindings, exported for other packages. They
expect pre-validated inputs and the index base differs between them (most take
1-based indices, `rs_diffuse()` takes 0-based background indices). Use the R
wrappers unless you have a concrete reason not to.
