# Metapath2Vec

Runs metapath2vec over a heterogeneous graph. Random walks are
constrained to follow a metapath schema over node types, e.g.
`c("gene", "pathway", "gene")`, and a skip-gram model is trained on
them. With `metapath_plus = TRUE`, negative samples are drawn from the
context node's own type (metapath2vec++).

Walks start only on nodes of the metapath's first type. Nodes no walk
ever reaches keep their random initialisation, so by default they are
removed and listed in the `unvisited_nodes` attribute. The function
warns when more than 10% of walks are dropped or walks reach less than
half the requested length on average: the schema does not fit the graph.
Details are in the `walk_stats` attribute. A metapath over a single type
warns as well; that is DeepWalk on a subgraph, so use
[`node2vec()`](https://gregorlueg.github.io/genewalkR/reference/node2vec.md).

## Usage

``` r
metapath2vec(
  graph_dt,
  node_dt,
  metapath,
  embd_dim = 8L,
  metapath2vec_params = params_metapath2vec(),
  metapath_plus = FALSE,
  filter_unvisited = TRUE,
  directed = FALSE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- graph_dt:

  data.table. The edge table. Needs to have the columns `"from"` and
  `"to"`, and can optionally have a `"weight"` column.

- node_dt:

  data.table. The node table with the columns `"id"` and `"type"`. Every
  edge endpoint needs to be in `id`.

- metapath:

  Character vector. The metapath over node types, closing on its
  starting type, e.g. `c("gene", "pathway", "gene")`.

- embd_dim:

  Integer. Size of the embedding dimensions to create. Defaults to `8L`.

- metapath2vec_params:

  Named list. The training parameters, see
  [`params_metapath2vec()`](https://gregorlueg.github.io/genewalkR/reference/params_metapath2vec.md).

- metapath_plus:

  Boolean. Use per-type negative sampling (metapath2vec++). Defaults to
  `FALSE`.

- filter_unvisited:

  Boolean. Remove the rows of nodes that no walk visited. Defaults to
  `TRUE`.

- directed:

  Boolean. Indicates if this is a directed or undirected network.
  Defaults to `FALSE`.

- seed:

  Integer. Seed for reproducibility.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A numeric matrix of `n_nodes x embd_dim` with the node ids as rownames.
Rows are ordered by node type. Carries three attributes:

- node_type - Named factor with the type of each row.

- unvisited_nodes - Character vector with the ids of the nodes no walk
  visited. Removed from the matrix if `filter_unvisited = TRUE`.

- walk_stats - List with `start_nodes`, `attempted`, `truncated`,
  `dropped`, `mean_length` and `walk_length`.

## References

Dong, Chawla and Swami, metapath2vec: Scalable Representation Learning
for Heterogeneous Networks, KDD 2017.
