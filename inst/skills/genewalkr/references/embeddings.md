# Embeddings: node2vec, metapath2vec, embedding drift

## node2vec

Grover and Leskovec 2016. Biased random walks plus a skip-gram model. Plain
function, no class:

```r
embd <- node2vec(
  graph_dt = edges,              # data.table with from, to
  embd_dim = 32L,
  node2vec_params = params_node2vec(walks_per_node = 20L, n_epochs = 10L),
  directed = FALSE,
  seed = 42L
)
```

Returns a nodes x `embd_dim` matrix with node IDs as rownames. `p` (return) and
`q` (in-out) in the bundle bias the walk; `p = q = 1` is DeepWalk.

It's unweighted at the moment, see the weight trap in `SKILL.md`.

Sanity check on known structure: `node2vec_test_data()` generates barbell,
caveman or stochastic block graphs with cluster labels, and
`evaluate_node2vec_test(embd, node_labels)` returns mean within- and
between-cluster cosine similarity. Within should beat between by a clear margin.

## metapath2vec

Dong, Chawla and Swami, KDD 2017. node2vec for heterogeneous graphs: walks
follow a schema over node types.

```r
embd <- metapath2vec(
  graph_dt = edges,
  node_dt = nodes,               # data.table with id, type
  metapath = c("gene", "pathway", "gene"),
  embd_dim = 32L,
  metapath2vec_params = params_metapath2vec(walks_per_node = 20L),
  metapath_plus = FALSE          # TRUE: per-type negative sampling (++)
)
attr(embd, "node_type")
attr(embd, "walk_stats")
attr(embd, "unvisited_nodes")
```

Rules the function enforces:

- The metapath closes on its starting type and has at least three entries.
- Every metapath type exists in `node_dt$type`, and types must not contain `-`.
- Every edge endpoint is in `node_dt$id`.

Walks start only on nodes of the first type. Nodes no walk reaches are dropped
from the matrix (`filter_unvisited = TRUE`) and listed in `unvisited_nodes`.
Rows come back ordered by node type, not input order.

Read the warnings. More than 10% of walks dropped, or walks reaching under half
the requested length, means the schema doesn't fit the graph:
`c("gene", "pathway", "pathway", "gene")` on a Reactome graph dead-ends: genes
are mostly annotated to leaf pathways, so the parent pathway the walk lands on
has no genes to hop back to. Check `walk_stats` and pick a schema the edges
support. `gene-pathway-gene` is the canonical choice and ignores gene-gene
edges entirely.

A single-type metapath warns too: that's DeepWalk on a subgraph, use
`node2vec()`.

## Embedding drift

Embed two networks over overlapping nodes (say, co-expression in two
conditions), align them with orthogonal Procrustes, and flag nodes whose
neighbourhood moved. The gene context drift idea from Cancer Cell 2025, see
`?generate_initial_embeddings` for the reference.

```r
obj <- EmbedDrift(graph_dt_1 = g1_edges, graph_dt_2 = g2_edges)
obj <- generate_initial_embeddings(
  obj,
  embd_dim = 8L,
  node2vec_params = params_node2vec(window_size = 2L)
)
obj <- calculate_drift(obj)
drift <- get_stats(obj)
```

`get_stats()` gives `node`, `cosine_similarity` and `node_status`:

- `"shared"`: cosine between the aligned embeddings, in `[-1, 1]`. Low means the
  node's context changed.
- `"g1_only"` and `"g2_only"`: sentinel values `-1.1` and `+1.1`, not
  similarities. Filter them out before ranking or plotting distributions.

Graph 1 is embedded with `seed`, graph 2 with `seed + 1L`, as independent
fits. There is no p-value: rank the shared nodes by `cosine_similarity`. With
`num_workers > 1L` the numbers move between runs.

`differential_graph_test_data()` builds a graph pair with known differential
nodes (`data_info$is_differential`) for checking a parameter choice before
running on real data.

Getters: `get_embeddings(obj, which = NULL)` (both, or `1L`/`2L`),
`get_graph_dts()`, `get_aligned_embd()`.
