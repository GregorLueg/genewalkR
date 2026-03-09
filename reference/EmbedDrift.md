# EmbedDrift

S7 class for context-aware differential node embedding analysis between
two graphs. Holds both input graphs, their embeddings, the
Procrustes-aligned embedding, and the resulting cosine similarity
statistics.

## Usage

``` r
EmbedDrift(graph_dt_1, graph_dt_2)
```

## Arguments

- graph_dt_1:

  A data.frame with columns `from`, `to`, and optionally `weight`
  representing the first graph.

- graph_dt_2:

  A data.frame with columns `from`, `to`, and optionally `weight`
  representing the second graph.

## Value

An initialised `EmbedDrift` object.

## Properties

- graph_dt_1:

  A data.frame containing the edges of graph 1. Must have columns
  `from`, `to`, and optionally `weight`.

- graph_dt_2:

  A data.frame containing the edges of graph 2. Must have columns
  `from`, `to`, and optionally `weight`.

- embd_1:

  Numeric matrix. Node embedding for graph 1, rows are nodes.

- embd_2:

  Numeric matrix. Node embedding for graph 2, rows are nodes.

- aligned_embd:

  Numeric matrix. Procrustes-aligned embedding of graph 1 onto graph 2,
  restricted to shared nodes.

- stats:

  A data.frame holding per-node cosine similarity scores and node
  metadata after running the differential analysis.

- params:

  Named list of parameters used during embedding and alignment.
