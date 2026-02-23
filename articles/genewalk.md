# Running GeneWalk

## Exploring Gene Walk

This vignette will first show how Gene Walk behaves on synthetic data,
what assumptions are baked in and then we will move onto using it with
real data. If you want to understand the method in more detail, please
check out [Ietswaart et
al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02264-8).

### Loading in libraries

``` r
library(genewalkR)
library(ggplot2)
library(data.table)
library(magrittr)
```

### Synthetic data

#### Exploring the data and initialising the object

Let’s start with some synthetic data to understand how this works… The
synthetic GeneWalk network is designed testing and benchmarking, with
signal genes forming coherent, degree-matched graph neighbourhoods and
noise genes spanning multiple ontology subtrees at random. The signal
genes will have some VERY clear signal. Noise genes also still due to
the nature of the synthetic data; however, to a lesser extent.

``` r
gene_walk_syn_data <- synthetic_genewalk_data()

str(gene_walk_syn_data)
#> List of 4
#>  $ full_data       :Classes 'data.table' and 'data.frame':   14224 obs. of  3 variables:
#>   ..$ from: chr [1:14224] "term_0001" "term_0001" "term_0002" "term_0002" ...
#>   ..$ to  : chr [1:14224] "term_0002" "term_0003" "term_0004" "term_0005" ...
#>   ..$ type: chr [1:14224] "hierarchy" "hierarchy" "hierarchy" "hierarchy" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ gene_to_pathways:Classes 'data.table' and 'data.frame':   9028 obs. of  3 variables:
#>   ..$ from: chr [1:9028] "gene_signal_0001" "gene_signal_0001" "gene_signal_0001" "gene_signal_0001" ...
#>   ..$ to  : chr [1:9028] "term_0025" "term_0019" "term_0024" "term_0017" ...
#>   ..$ type: chr [1:9028] "part_of" "part_of" "part_of" "part_of" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ gene_ids        : chr [1:600] "gene_signal_0001" "gene_signal_0002" "gene_signal_0003" "gene_signal_0004" ...
#>  $ pathway_ids     : chr [1:355] "term_0001" "term_0002" "term_0003" "term_0004" ...
```

The data contains everything we need to initialise a new GeneWalk:

- The full_data with the term ontology, gene to gene edges and gene to
  terms, mimicking relevant inputs.
- The genes to pathways for testing later.
- The gene identifiers including in the run. In this case, we have a set
  of “signal” genes that serve as a positive control and “noise genes”
  that are randomly distributed.
- The term/pathway identifiers.

Let’s initialise the class.

``` r
genewalk_obj <- GeneWalk(
  graph_dt = gene_walk_syn_data$full_data,
  gene_to_pathway_dt = gene_walk_syn_data$gene_to_pathways,
  gene_ids = gene_walk_syn_data$gene_ids,
  pathway_ids = gene_walk_syn_data$pathway_ids
)

genewalk_obj
#> GeneWalk
#>   Represented genes gene_signal_0001 | gene_signal_0002 | gene_signal_0003 ; Total of 600 genes. 
#>   Number of edges: 14224 
#>   Embedding generated: no 
#>   Permutations generated: no 
#>   Statistics calculated: no
```

#### Running gene walk on the synthetic data

This function here generates the node2vec embedding based on the network

``` r
genewalk_obj <- generate_initial_emb(
  genewalk_obj,
  node2vec_params = params_node2vec(),
  .verbose = TRUE
)
```

We need to generate a background distribution for testing purposes. This
will generate three random permuted networks as in the paper. These will
be used for statistical testing.

``` r
genewalk_obj <- generate_permuted_emb(genewalk_obj, .verbose = TRUE)
```

We now need to compare the similarity of the Cosine similarities between
the gene embeddings to the pathway embeddings and check how often they
are larger than the ones from the background embedding, giving us the
p-values.

``` r
genewalk_obj <- calculate_genewalk_stats(
  genewalk_obj,
  .verbose = TRUE
)
```

Now we can extract the statistics:

``` r
statistics <- get_stats(genewalk_obj)

head(statistics)
#>                gene   pathway similarity     avg_pval pval_ci_lower
#>              <char>    <char>      <num>        <num>         <num>
#> 1:  gene_noise_0299 term_0186  0.9897635 3.606504e-13  3.845573e-20
#> 2: gene_signal_0060 term_0018  0.9879104 3.606504e-13  3.845573e-20
#> 3: gene_signal_0108 term_0013  0.9866632 3.606504e-13  3.845573e-20
#> 4: gene_signal_0154 term_0007  0.9860725 3.606504e-13  3.845573e-20
#> 5: gene_signal_0062 term_0023  0.9859357 3.606504e-13  3.845573e-20
#> 6:  gene_noise_0166 term_0335  0.9855766 3.606504e-13  3.845573e-20
#>    pval_ci_upper avg_global_fdr global_fdr_ci_lower global_fdr_ci_upper
#>            <num>          <num>               <num>               <num>
#> 1:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#> 2:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#> 3:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#> 4:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#> 5:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#> 6:  3.382297e-06   3.755742e-12        4.057602e-21         0.003476339
#>    avg_gene_fdr gene_fdr_ci_lower gene_fdr_ci_upper
#>           <num>             <num>             <num>
#> 1: 9.089340e-13      1.583855e-20      5.216142e-05
#> 2: 8.895869e-13      1.616901e-20      4.894331e-05
#> 3: 8.481449e-13      1.692655e-20      4.249831e-05
#> 4: 9.453335e-13      1.525276e-20      5.858976e-05
#> 5: 8.022017e-13      1.785590e-20      3.604005e-05
#> 6: 8.895869e-13      1.616901e-20      4.894331e-05
```

The different metrics in the data Let’s explore the signal in the
synthetic data a bit… What are we observing … ?

``` r
# add the labels for signal and noise
statistics[, signal := grepl("signal", gene)][, 
  signal := factor(signal, levels = c("TRUE", "FALSE"))
]

# plot the gene <> term similarities
ggplot(
  data = statistics,
  mapping = aes(x = similarity)
) + 
  geom_histogram(bins = 45, fill = "lightgrey") +
  facet_wrap(~ signal) +
  xlab("Cosine similarity") +
  ylab("Count") +
  theme_bw()
```

![](genewalk_files/figure-html/plot%20the%20similarities%20-%20synthetic-1.png)

As we can appreciate, the Cosine similarities between the genes and
pathways in the signal data set are all much higher compared to the
FALSE ones. Let’s look at the p-values

``` r
ggplot(
  data = statistics,
  mapping = aes(x = avg_pval)
) + 
  geom_histogram(bins = 45, fill = "lightgrey") +
  facet_wrap(~ signal) +
  xlab("Cosine similarity") +
  ylab("Count") +
  theme_bw()
```

![](genewalk_files/figure-html/plot%20the%20p-values%20-%20synthetic-1.png)

Similar pattern. The `TRUE` signal ones are nearly all significant. Less
so for the noisy signal. The signal in the noise data is there, but this
is expected due to the data. Let’s look at genes passing the FDR
thresholds

``` r
table(
  grepl("signal", statistics$gene),
  statistics$avg_global_fdr < 0.05
)
#>        
#>         FALSE TRUE
#>   FALSE  2817 2312
#>   TRUE    420 3479
```

Basically, all of our signal genes show significance associations
between gene and pathways; less so, for the noise genes.

### Real data

Let’s explore real data now. The package provides a builder factory to
generate the objects. Within the package, there is a DuckDB that
contains

**Network resources**

- The STRING network extracted from OpenTargets.
- The SIGNOR network extracted from OpenTargets.
- The Reactome gene to gene network extracted from OpenTargets.
- The Intact network extracted from OpenTargets.
- A combined network from the sources above, based on the approach from
  [Barrio-Hernandez, et
  al.](https://www.nature.com/articles/s41588-023-01327-9).

**Pathway terms**

- The Gene Ontology data extracted from the OBO files from the OBO
  foundry and OpenTargets.
- The Reactome pathway ontology and their gene to pathway associations
  from OpenTargets.

#### Using the builder factory

Let’s use the Gene Ontology and combined network for this example.

``` r
gw_factory <- GeneWalkGenerator$new()

gw_factory$add_pathways() # will add GO to the builder
gw_factory$add_ppi(source = "combined") # will add the combined one
gw_factory$build() # will load the data into the factory
#> Downloading database...
#> Download complete
#> Built network with 1343553 edges and 61047 nodes
```

The idea of the factory is to easily iterate through various bags of
genes of interest. Now let’s use the factory to look specifically at the
MYC target genes (provided in the package) and generate a
GeneWalkNetwork for them.

``` r
data(myc_genes)

myc_gwn <- gw_factory$create_for_genes(genes = myc_genes$ensembl_gene)

myc_gwn
#> GeneWalk
#>   Represented genes ENSG00000004779 | ENSG00000013275 | ENSG00000041357 ; Total of 200 genes. 
#>   Number of edges: 86056 
#>   Embedding generated: no 
#>   Permutations generated: no 
#>   Statistics calculated: no
```

#### Running gene walk on actual data

We can now use the same steps as above

``` r
data(myc_genes)

myc_gwn <- gw_factory$create_for_genes(genes = myc_genes$ensembl_gene)

# this will take a while
myc_gwn <- generate_initial_emb(
  myc_gwn,
  node2vec_params = params_node2vec(
    n_epochs = 5L,
    walks_per_node = 20L,
    walk_length = 20L
  ),
  .verbose = TRUE
)
# if you are on a many core machines, consider upping the num_workers here
# the default is 4L, but you can up this on your machine. also, for 
# speed purposes, the walks per node and walk length has been reduced
# to 20L each.

# this even longer... go get coffee
myc_gwn <- generate_permuted_emb(
  myc_gwn,
  .verbose = TRUE
)

# let's calculate the statistics
myc_gwn <- calculate_genewalk_stats(
  myc_gwn,
  .verbose = TRUE
)

# let's extract the results
myc_gwn_res <- get_stats(myc_gwn)
```

And look at some of the genes

``` r
# translate go ids to names and do the same for the gene symbols
gene_symbol_translation <- setNames(
  myc_genes$gene_symbol,
  myc_genes$ensembl_gene
)

# get the go data
go_info <- get_gene_ontology_info()

go_id_translation <- setNames(
  go_info$go_name,
  go_info$go_id
)

myc_gwn_res_translated <- copy(myc_gwn_res)[, `:=`(
  gene = gene_symbol_translation[gene],
  pathway = go_id_translation[pathway]
)]

head(myc_gwn_res_translated, 15L)
#>       gene                                                           pathway
#>     <char>                                                            <char>
#>  1:   TCP1                                             protein stabilization
#>  2:  PSMA1       proteasomal ubiquitin-independent protein catabolic process
#>  3:   TCP1                                          unfolded protein binding
#>  4:  PSMB2       proteasomal ubiquitin-independent protein catabolic process
#>  5: HNRNPU                                                    lncRNA binding
#>  6:  H2AZ1                                                       euchromatin
#>  7:  PSMB2                                        spermatoproteasome complex
#>  8:  PSMA2       proteasomal ubiquitin-independent protein catabolic process
#>  9:  PSMB3       proteasomal ubiquitin-independent protein catabolic process
#> 10:  PSMA4       proteasomal ubiquitin-independent protein catabolic process
#> 11:   MCM5                                       3'-5' DNA helicase activity
#> 12:   RPS5                                                cytosolic ribosome
#> 13:  PSMA6       proteasomal ubiquitin-independent protein catabolic process
#> 14:  PSMB3 proteasome-mediated ubiquitin-dependent protein catabolic process
#> 15: SNRPD1                                              spliceosomal complex
#>     similarity     avg_pval pval_ci_lower pval_ci_upper avg_global_fdr
#>          <num>        <num>         <num>         <num>          <num>
#>  1:  0.9951420 0.000000e+00  0.000000e+00  0.000000e+00   0.000000e+00
#>  2:  0.9941487 1.477964e-10  1.321279e-16  1.653229e-04   1.104924e-08
#>  3:  0.9929487 2.016823e-07  1.282315e-07  3.172055e-07   1.186081e-04
#>  4:  0.9923260 2.016823e-07  1.282315e-07  3.172055e-07   1.186081e-04
#>  5:  0.9919932 2.541038e-07  2.541038e-07  2.541038e-07   1.186081e-04
#>  6:  0.9918907 2.541038e-07  2.541038e-07  2.541038e-07   1.186081e-04
#>  7:  0.9913864 2.541038e-07  2.541038e-07  2.541038e-07   1.186081e-04
#>  8:  0.9913607 2.541038e-07  2.541038e-07  2.541038e-07   1.186081e-04
#>  9:  0.9911833 2.908762e-07  2.231834e-07  3.791006e-07   1.186081e-04
#> 10:  0.9910046 3.947798e-07  2.348029e-07  6.637526e-07   1.505605e-04
#> 11:  0.9908384 3.947798e-07  2.348029e-07  6.637526e-07   1.505605e-04
#> 12:  0.9903601 5.564262e-07  3.417126e-07  9.060541e-07   1.523113e-04
#> 13:  0.9903355 5.564262e-07  3.417126e-07  9.060541e-07   1.523113e-04
#> 14:  0.9900156 5.817524e-07  3.288170e-07  1.029253e-06   1.523113e-04
#> 15:  0.9899194 5.817524e-07  3.288170e-07  1.029253e-06   1.523113e-04
#>     global_fdr_ci_lower global_fdr_ci_upper avg_gene_fdr gene_fdr_ci_lower
#>                   <num>               <num>        <num>             <num>
#>  1:        0.000000e+00        0.0000000000 0.000000e+00      0.000000e+00
#>  2:        1.443486e-16        0.8457693680 1.757330e-09      1.389487e-16
#>  3:        8.349323e-05        0.0001684912 3.025235e-06      1.923473e-06
#>  4:        8.349323e-05        0.0001684912 4.827972e-06      4.827972e-06
#>  5:        8.349323e-05        0.0001684912 2.566448e-05      2.566448e-05
#>  6:        8.349323e-05        0.0001684912 5.082076e-06      5.082076e-06
#>  7:        8.349323e-05        0.0001684912 4.827972e-06      4.827972e-06
#>  8:        8.349323e-05        0.0001684912 9.927350e-06      8.226341e-06
#>  9:        8.349323e-05        0.0001684912 6.518872e-06      4.395730e-06
#> 10:        6.917179e-05        0.0003277124 8.657224e-06      3.782931e-06
#> 11:        6.917179e-05        0.0003277124 6.690153e-06      3.781395e-06
#> 12:        6.846424e-05        0.0003388446 1.391065e-05      8.542815e-06
#> 13:        6.846424e-05        0.0003388446 2.837773e-05      1.742734e-05
#> 14:        6.846424e-05        0.0003388446 8.213264e-06      3.588934e-06
#> 15:        6.846424e-05        0.0003388446 2.734236e-05      1.545440e-05
#>     gene_fdr_ci_upper
#>                 <num>
#>  1:      0.000000e+00
#>  2:      2.222554e-02
#>  3:      4.758083e-06
#>  4:      4.827972e-06
#>  5:      2.566448e-05
#>  6:      5.082076e-06
#>  7:      4.827972e-06
#>  8:      1.198009e-05
#>  9:      9.667493e-06
#> 10:      1.981203e-05
#> 11:      1.183641e-05
#> 12:      2.265135e-05
#> 13:      4.620876e-05
#> 14:      1.879602e-05
#> 15:      4.837490e-05
```

#### Exploring the results

A simple way to visualise the initial results is to use the
plot_results() function. This will tell you the connectivity of the
genes within your bag of genes (degree + 1 on the x-axis), the number of
`gene <> pathway` connections that are significant for this specific
gene and the number of `gene <> pathway` connections.

``` r
plot_gw_results(myc_gwn)
```

Other options are to plot for individual genes the significantly
associated terms.

### Using your own data

What if you do not want to use the provided data … ? In this case, you
can use this simple wrapper class to help you. Let’s show case this in
terms of the internally stored Reactome data.

``` r
# these are helpers to get the reactome data from the DuckDB in the
# package
reactome_genes <- get_gene_to_reactome()
reactome_ppi <- get_interactions_reactome()
reactome_hierarchy <- get_reactome_hierarchy(relationship = "parent_of")

data_builder <- new_data_builder(
  ppis = reactome_ppi,
  gene_to_pathways = reactome_genes,
  pathway_hierarchy = reactome_hierarchy
)
```

If you want to now get the full graph dt including everything for
verification:

``` r
get_gw_data(data_builder) %>% head()
#>               from              to        type
#>             <char>          <char>      <char>
#> 1: ENSG00000116030 ENSG00000069667 interaction
#> 2: ENSG00000115904 ENSG00000113070 interaction
#> 3: ENSG00000115904 ENSG00000168546 interaction
#> 4: ENSG00000115904 ENSG00000091181 interaction
#> 5: ENSG00000116030 ENSG00000084676 interaction
#> 6: ENSG00000115904 ENSG00000169855 interaction
```

If you want to subset to your bag of genes of interest, you can run the
following:

``` r
genes_of_interest <- reactome_genes[to == "R-HSA-1989781", from]

gwr_data <- get_gw_data_filtered(data_builder, genes_of_interest)

str(gwr_data)
#> List of 4
#>  $ gwn                 :Classes 'data.table' and 'data.frame':   4007 obs. of  3 variables:
#>   ..$ from: chr [1:4007] "ENSG00000132170" "ENSG00000196498" "ENSG00000196498" "ENSG00000204231" ...
#>   ..$ to  : chr [1:4007] "ENSG00000132522" "ENSG00000025434" "ENSG00000177565" "ENSG00000141027" ...
#>   ..$ type: chr [1:4007] "interaction" "interaction" "interaction" "interaction" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ genes_to_pathways   :Classes 'data.table' and 'data.frame':   898 obs. of  3 variables:
#>   ..$ from: chr [1:898] "ENSG00000001167" "ENSG00000001167" "ENSG00000001167" "ENSG00000001167" ...
#>   ..$ to  : chr [1:898] "R-HSA-9614657" "R-HSA-1989781" "R-HSA-380994" "R-HSA-381183" ...
#>   ..$ type: chr [1:898] "part_of" "part_of" "part_of" "part_of" ...
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ represented_genes   : chr [1:116] "ENSG00000132170" "ENSG00000196498" "ENSG00000204231" "ENSG00000163586" ...
#>  $ represented_pathways: chr [1:2825] "R-HSA-9614657" "R-HSA-1989781" "R-HSA-380994" "R-HSA-381183" ...
```

The data can be easily supplied now:

``` r
custom_gwn <- with(gwr_data, 
  GeneWalk(
    graph_dt = gwn, 
    gene_to_pathway_dt = genes_to_pathways, 
    gene_ids = represented_genes,
    pathway_ids = represented_pathways
  )
)
```

From here on, you can run the whole approach on your custom network.
