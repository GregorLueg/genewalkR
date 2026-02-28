# Node2vec

## Running node2vec

### Loading in libraries

``` r
library(genewalkR)
library(data.table)
library(magrittr)
library(ggplot2)
```

### Synthetic data

genewalkR provides quite a few synthetic graph data sets for testing and
exploration purposes, namely:

- **Barbell graph**: a graph of two distinct clusters
- **Cavemen graph**: a graph with distinct communities with a low
  probability of connection between them.
- **Stochastic blocks**: a more densely connected graph with
  communities; however, more cross-edges compared to the previous
  graphs.

#### Barbell graph

Let’s load in the barbell graph.

``` r
# this creates a random data set with 20 nodes per barbell cluster
barbell_data <- node2vec_test_data(
  test_type = "barbell",
  n_nodes_per_cluster = 20L
)
```

##### Exploring the data

Let’s transform this into an igraph first and just look at the data like
this

``` r
g <- igraph::graph_from_data_frame(d = barbell_data$edges, directed = FALSE)

plot(g, vertex.label = NA)
```

![](node2vec_files/figure-html/generate%20igraph%20from%20barbell%20data-1.png)

##### Running node2vec

You can appreciate a nice barbell graph with two distinct groups. Now
let’s run node2vec on this data. The node2vec implementation is based on
[node2vec-rs](https://github.com/GregorLueg/node2vec-rs), a Rust-based
implementation of node2vec. In its default setting, we use a HIGHLY
optimised version that runs on the CPU with SIMD-acceleration and a
tolerated race condition (undefined behaviour in low level terms) that
makes this very fast. The issue with this is that the underlying
embeddings will not be deterministic; however, the structure of the
underlying graph will be identified. If you really, really, really need
deterministic results here, set the `n_workers = 1L`.

``` r
# run node2vec
barbell_node2vec_results <- node2vec(
  graph_dt = barbell_data$edges,
  .verbose = TRUE
)

# the parameters are controlled via a supplied parameter list, see code here
str(params_node2vec())
#> List of 10
#>  $ p             : num 1
#>  $ q             : num 1
#>  $ walks_per_node: int 40
#>  $ walk_length   : int 40
#>  $ num_workers   : int 2
#>  $ batch_size    : int 256
#>  $ n_epochs      : int 20
#>  $ n_negatives   : int 5
#>  $ window_size   : int 2
#>  $ lr            : num 0.01

# this list can be supplied and tells node2vec how to run
```

If you want more details on how node2vec actually works, please, check
the original paper from [Grover and
Leskovec](https://arxiv.org/abs/1607.00653). Let’s project that data
onto two dimensions via PCA and visualise the results

``` r
# run pca on the embeddings
pca_results <- prcomp(barbell_node2vec_results, scale. = TRUE)

pca_dt <- as.data.table(pca_results$x[, 1:2], keep.rownames = "node") %>%
  merge(., barbell_data$node_labels, by = "node")

pca_dt[, cluster := factor(cluster)]

ggplot(data = pca_dt, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = cluster)) 
```

![](node2vec_files/figure-html/node2vec%20pca-1.png)

As you can appreciate, node2vec embeddings nicely identify the two
groups in the underlying graph.

### Other graphs

Let’s have a look at some other graphs with well known structures

#### Cavemen graph

Visualisation of the graph structure in igraph first:

``` r
# generate the cavemen
caveman_data <- node2vec_test_data(
  test_type = "cavemen", n_clusters = 5L, p_between = 0.5
)

g <- igraph::graph_from_data_frame(d = caveman_data$edges, directed = FALSE)

plot(g, vertex.label = NA)
```

![](node2vec_files/figure-html/cavemen%20igraph-1.png)

Let’s check out node2vec on this graph

``` r
# run node2vec
cavemen_node2vec_results <- node2vec(
  graph_dt = caveman_data$edges,
  .verbose = TRUE
)

# run pca on the embeddings
pca_results <- prcomp(cavemen_node2vec_results, scale. = TRUE)

pca_dt <- as.data.table(pca_results$x[, 1:2], keep.rownames = "node") %>%
  merge(., caveman_data$node_labels, by = "node")

pca_dt[, cluster := factor(cluster)]

ggplot(data = pca_dt, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = cluster)) 
```

![](node2vec_files/figure-html/cavemen%20node2vec%20visualisation-1.png)

#### Stochastic blocks

In this graph, we have less defined communities with more between edges
across the communities.

``` r
# generate the stochastic block graph
stochastic_block_data <- node2vec_test_data(
  test_type = "stochastic_block", n_clusters = 3L, p_between = 0.05
)

g <- igraph::graph_from_data_frame(
  d = stochastic_block_data$edges, 
  directed = FALSE
)

plot(g, vertex.label = NA)
```

![](node2vec_files/figure-html/stochastic%20block%20igraph-1.png)

Let’s check out node2vec on this graph

``` r
# run node2vec
stochastic_blocks_node2vec_results <- node2vec(
  graph_dt = stochastic_block_data$edges,
  .verbose = TRUE
)

# run pca on the embeddings
pca_results <- prcomp(stochastic_blocks_node2vec_results, scale. = TRUE)

pca_dt <- as.data.table(pca_results$x[, 1:2], keep.rownames = "node") %>%
  merge(., stochastic_block_data$node_labels, by = "node")

pca_dt[, cluster := factor(cluster)]

ggplot(data = pca_dt, mapping = aes(x = PC1, y = PC2)) +
  geom_point(mapping = aes(color = cluster)) 
```

![](node2vec_files/figure-html/stochastic%20block%20node2vec%20visualisation-1.png)

In conclusion, node2vec generates an embedding of the underlying network
structure that you can then use to visualise the data or do additional
analysis (supervised tasks or like in this package statistical testing
for gene \<\> pathway associations in a bag of genes). Pending on your
needs, you can go high on the available workers and make this fast, or
if you really, really, really need similar embeddings every single time,
you can also choose to disable the parallelism during training.
