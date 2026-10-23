---
name: genewalkr
description: How to USE the genewalkR R package for graph methods in computational biology. Covers installation, GeneWalk (Ietswaart et al.) over GO/Reactome plus PPI networks, node2vec and metapath2vec embeddings, embedding drift between two networks (Procrustes alignment), heat diffusion kernels (raw, z-score, Monte Carlo) and diffusion profiles (constrained personalised PageRank, Ruiz et al.), and the bundled DuckDB of gene, pathway and interaction data. Use this whenever a task involves genewalkR or any of its classes and functions (GeneWalk, GeneWalkGenerator, EmbedDrift, DiffusionScores, DiffusionProfiles, node2vec, metapath2vec, build_kernel, diffuse_scores, generate_profiles, get_interactions_*, params_* bundles), or when a user in R asks for GeneWalk, gene-function relevance in a gene list, node2vec or metapath2vec embeddings, network propagation or guilt-by-association scoring, drug-disease diffusion profiles, or comparing gene context between two networks. Apply even when the package is not named but the user is working in an R session where genewalkR is loaded.
---

# Using genewalkR

genewalkR is an R package over a Rust crate (extendr). Graph methods for
computational biology: GeneWalk, node2vec, metapath2vec, embedding drift, heat
diffusion and diffusion profiles. The numerics run in Rust; the R layer
validates, marshals and stores results.

The API is small (~80 exports), but don't guess signatures. Check
`references/api-index.md` for what exists, then `?fn` for how to call it.

## The three things to know before writing any code

**1. Tunables come in `params_*()` bundles.** One validated list per method:

```r
node2vec(edges, embd_dim = 32L, node2vec_params = params_node2vec(
  walks_per_node = 20L, n_epochs = 10L
))
```

Each bundle has a matching checkmate assertion, so a malformed list fails at the
call site. There are six of them, see `references/conventions.md`.

**2. Analysis objects are S7 and each step returns the object.** Constructor,
then verbs, then getters. Reassign every step: `obj <- build_kernel(obj)`.
Nothing mutates in place. A step whose prerequisite is missing warns and hands
back the object unchanged (or errors), it does not compute the prerequisite for
you.

**3. Print the object.** `GeneWalk`, `EmbedDrift`, `DiffusionScores`,
`DiffusionProfiles` and `GeneWalkGenerator` all print which stages have run.

## Where to look

| Task | Read |
|---|---|
| Installing, Rust toolchain, the data download | `references/install.md` |
| `params_*()`, S7 chains, getters, edge table formats, threading | `references/conventions.md` |
| GeneWalk: building networks from the DB, running it, reading the stats | `references/genewalk.md` |
| node2vec, metapath2vec, embedding drift between two graphs | `references/embeddings.md` |
| Diffusion kernels and scores, diffusion profiles, `dedup_edges()` | `references/diffusion.md` |
| The DuckDB: genes, GO, Reactome, STRING, IntAct, SIGNOR, PathwayCommons | `references/data.md` |
| Does function X exist? What's it called? | `references/api-index.md` |

## Smoke test

Synthetic data, no downloads. Use it to confirm an install works before
debugging anything else.

```r
library(genewalkR)

syn <- synthetic_genewalk_data(
  data_params = params_genewalk_data(n_signal_genes = 30L, n_noise_genes = 30L)
)

gw <- GeneWalk(
  graph_dt = syn$full_data,
  gene_to_pathway_dt = syn$gene_to_pathways,
  gene_ids = syn$gene_ids,
  pathway_ids = syn$pathway_ids
)

gw <- generate_initial_emb(
  gw,
  n_graph = 2L,
  genewalk_params = params_genewalk(walks_per_node = 20L),
  .verbose = FALSE
)
gw <- generate_permuted_emb(gw, n_perm = 2L, .verbose = FALSE)
gw <- calculate_genewalk_stats(gw, .verbose = FALSE)

gw # all three stages should say "yes"
head(get_stats(gw))
```

## Traps that bite everywhere

- **Integers are integers.** Every count goes through `checkmate::qassert(x,
  "I1")`, so `embd_dim = 8` fails where `8L` works.
- **Gene identifiers are Ensembl.** The DB and `myc_genes$ensembl_gene` use
  `ENSG...` IDs. Feed symbols to `$create_for_genes()` and every gene gets
  dropped with a warning, then `GeneWalk()` errors on empty `gene_ids`.
  Translate with `get_gene_info()` first.
- **Edge weights go in a column called `weight`.** Every edge-table function
  reads that name. Anything else (`weights`, `score`) is silently ignored and
  the graph is treated as unweighted. `DiffusionScores` is the exception: it
  takes an igraph and builds an unweighted adjacency, ignoring edge weights.
- **Multi-threaded training is not reproducible.** `params_genewalk()` defaults
  to `num_workers = 1L` for that reason. `params_node2vec()` and
  `params_metapath2vec()` default to all cores minus two. Set
  `num_workers = 1L` when the exact numbers matter.
- **The diffusion kernel is an external pointer.** It does not survive
  `saveRDS()` or a session restart. Rebuild with `build_kernel()`.
- **This is not the Python GeneWalk.** Same method, separate implementation.
  Don't port the Python CLI arguments or output column names across.
