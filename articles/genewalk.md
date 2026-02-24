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
#> 1: gene_signal_0060 term_0018  0.9918368 0.000000e+00  0.000000e+00
#> 2:  gene_noise_0299 term_0186  0.9876428 3.606504e-13  3.845573e-20
#> 3: gene_signal_0001 term_0007  0.9853014 3.606504e-13  3.845573e-20
#> 4:  gene_noise_0228 term_0223  0.9850876 3.606504e-13  3.845573e-20
#> 5:  gene_noise_0166 term_0335  0.9838016 3.606504e-13  3.845573e-20
#> 6: gene_signal_0198 term_0007  0.9837879 3.606504e-13  3.845573e-20
#>    pval_ci_upper avg_global_fdr global_fdr_ci_lower global_fdr_ci_upper
#>            <num>          <num>               <num>               <num>
#> 1:  0.000000e+00   0.000000e+00        0.000000e+00         0.000000000
#> 2:  3.382297e-06   3.412309e-12        4.448859e-21         0.002617267
#> 3:  3.382297e-06   3.412309e-12        4.448859e-21         0.002617267
#> 4:  3.382297e-06   3.412309e-12        4.448859e-21         0.002617267
#> 5:  3.382297e-06   3.412309e-12        4.448859e-21         0.002617267
#> 6:  3.382297e-06   3.412309e-12        4.448859e-21         0.002617267
#>    avg_gene_fdr gene_fdr_ci_lower gene_fdr_ci_upper
#>           <num>             <num>             <num>
#> 1: 0.000000e+00      0.000000e+00      0.000000e+00
#> 2: 9.089340e-13      1.583855e-20      5.216142e-05
#> 3: 9.453335e-13      1.525276e-20      5.858976e-05
#> 4: 8.258122e-13      1.736564e-20      3.927098e-05
#> 5: 8.895869e-13      1.616901e-20      4.894331e-05
#> 6: 8.022017e-13      1.785590e-20      3.604005e-05
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
#>   FALSE  2741 2388
#>   TRUE    383 3516
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
#>       gene                                                     pathway
#>     <char>                                                      <char>
#>  1:   TCP1                                       protein stabilization
#>  2: HNRNPU                                              lncRNA binding
#>  3:  PSMA1 proteasomal ubiquitin-independent protein catabolic process
#>  4:  PSMB2 proteasomal ubiquitin-independent protein catabolic process
#>  5:  PSMA2 proteasomal ubiquitin-independent protein catabolic process
#>  6:   MCM5                                 3'-5' DNA helicase activity
#>  7:  PSMB3 proteasomal ubiquitin-independent protein catabolic process
#>  8:  PSMA4 proteasomal ubiquitin-independent protein catabolic process
#>  9:  PSMB2                                  spermatoproteasome complex
#> 10:   TCP1                                    unfolded protein binding
#> 11:   PHB2                                                    membrane
#> 12:  H2AZ1                                                 euchromatin
#> 13:   MCM5                       single-stranded DNA helicase activity
#> 14:  PSMA4                                  spermatoproteasome complex
#> 15:  U2AF1                                               nuclear speck
#>     similarity     avg_pval pval_ci_lower pval_ci_upper avg_global_fdr
#>          <num>        <num>         <num>         <num>          <num>
#>  1:  0.9967965 0.000000e+00  0.000000e+00  0.000000e+00   0.000000e+00
#>  2:  0.9943855 1.082080e-13  1.219902e-19  9.598283e-08   1.124543e-12
#>  3:  0.9943355 1.082080e-13  1.219902e-19  9.598283e-08   1.124543e-12
#>  4:  0.9929860 1.082080e-13  1.219902e-19  9.598283e-08   1.124543e-12
#>  5:  0.9926530 1.082080e-13  1.219902e-19  9.598283e-08   1.124543e-12
#>  6:  0.9920835 1.561071e-13  8.585570e-20  2.838417e-07   1.333312e-12
#>  7:  0.9920716 1.561071e-13  8.585570e-20  2.838417e-07   1.333312e-12
#>  8:  0.9917253 1.691846e-10  1.314014e-16  2.178320e-04   9.446964e-09
#>  9:  0.9915933 1.691846e-10  1.314014e-16  2.178320e-04   9.446964e-09
#> 10:  0.9914230 1.862118e-10  1.306065e-16  2.654909e-04   9.446964e-09
#> 11:  0.9913022 2.005906e-10  1.298366e-16  3.099016e-04   9.446964e-09
#> 12:  0.9912334 2.005906e-10  1.298366e-16  3.099016e-04   9.446964e-09
#> 13:  0.9903481 2.308686e-07  7.160984e-08  7.443154e-07   6.692991e-05
#> 14:  0.9902831 2.308686e-07  7.160984e-08  7.443154e-07   6.692991e-05
#> 15:  0.9902802 2.308686e-07  7.160984e-08  7.443154e-07   6.692991e-05
#>     global_fdr_ci_lower global_fdr_ci_upper avg_gene_fdr gene_fdr_ci_lower
#>                   <num>               <num>        <num>             <num>
#>  1:        0.000000e+00        0.000000e+00 0.000000e+00      0.000000e+00
#>  2:        1.291182e-20        9.794106e-05 5.042912e-13      2.787771e-20
#>  3:        1.291182e-20        9.794106e-05 3.733694e-13      3.719794e-20
#>  4:        1.291182e-20        9.794106e-05 3.640287e-13      3.811325e-20
#>  5:        1.291182e-20        9.794106e-05 3.793459e-13      3.663546e-20
#>  6:        1.096483e-20        1.621293e-04 4.441322e-13      3.149173e-20
#>  7:        1.096483e-20        1.621293e-04 5.204167e-13      2.704825e-20
#>  8:        1.398350e-16        6.382174e-01 1.544285e-09      1.349056e-16
#>  9:        1.398350e-16        6.382174e-01 1.204657e-09      1.368889e-16
#> 10:        1.398350e-16        6.382174e-01 1.132578e-09      1.357355e-16
#> 11:        1.398350e-16        6.382174e-01 2.267252e-09      1.368577e-16
#> 12:        1.398350e-16        6.382174e-01 1.477965e-09      1.356039e-16
#> 13:        1.937908e-05        2.311571e-04 2.654989e-06      8.235132e-07
#> 14:        1.937908e-05        2.311571e-04 4.501938e-06      1.396392e-06
#> 15:        1.937908e-05        2.311571e-04 3.924766e-06      1.217367e-06
#>     gene_fdr_ci_upper
#>                 <num>
#>  1:      0.000000e+00
#>  2:      9.122328e-06
#>  3:      3.747645e-06
#>  4:      3.476924e-06
#>  5:      3.927979e-06
#>  6:      6.263658e-06
#>  7:      1.001298e-05
#>  8:      1.767767e-02
#>  9:      1.060129e-02
#> 10:      9.450243e-03
#> 11:      3.756040e-02
#> 12:      1.610853e-02
#> 13:      8.559627e-06
#> 14:      1.451415e-05
#> 15:      1.265336e-05
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
