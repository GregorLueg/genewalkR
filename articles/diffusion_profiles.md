# Diffusion profiles

## Diffusion profiles

### Setup

``` r

library(genewalkR)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
library(magrittr)
```

### Background

[Ruiz, Zitnik and
Leskovec](https://www.nature.com/articles/s41467-021-21770-8) built the
multiscale interactome: drugs and diseases hang off proteins, proteins
talk to each other and to a hierarchy of biological functions. To
explain how a drug treats a disease, they run one personalised PageRank
per drug and per disease and compare the resulting visitation vectors,
the diffusion profiles. Two drugs that hit the same biology end up with
similar profiles, and a drug whose profile looks like a disease’s is a
candidate treatment.

It’s a personalised PageRank with two twists:

- **Sinks.** The seed is a source. Every other node of a sink type
  absorbs the walker: mass flows in, nothing flows out. In the paper
  drugs and diseases are sinks, so a drug’s profile doesn’t leak through
  some other drug into its targets.
- **Type weights.** At each step the walker first picks one of the node
  types among its neighbours, proportional to the type weight, then a
  neighbour of that type. This lets you push the walk towards, say,
  biological functions, no matter how many protein neighbours a node
  has.

Leave `type_weights = NULL` and `sink_types` empty and you get the plain
(edge-weighted) personalised PageRank. The graph is built once in Rust
and the seeds run in parallel.

### A toy interactome

Two drugs, two diseases, six proteins in a chain and three functions.
`drug_A` and `disease_X` share `P2`; `drug_B` and `disease_Y` sit at the
other end of the chain.

``` r

toy_nodes <- data.table(
  id = c(
    "drug_A",
    "drug_B",
    paste0("P", 1:6),
    "F1",
    "F2",
    "F3",
    "disease_X",
    "disease_Y"
  ),
  type = c(
    rep("drug", 2),
    rep("protein", 6),
    rep("function", 3),
    rep("disease", 2)
  )
)

toy_edges <- data.table(
  from = c(
    # drugs and diseases to proteins
    "drug_A", "drug_A", "drug_B", "disease_X", "disease_X", "disease_Y",
    # protein chain
    "P1", "P2", "P3", "P4", "P5",
    # proteins to functions
    "P1", "P2", "P3", "P4", "P5", "P6"
  ),
  to = c(
    "P1", "P2", "P5", "P2", "P3", "P6",
    "P2", "P3", "P4", "P5", "P6",
    "F1", "F1", "F2", "F2", "F3", "F3"
  )
)

toy_seeds <- c("drug_A", "drug_B", "disease_X", "disease_Y")
```

First the plain personalised PageRank: no sinks, no type weights.

``` r

plain <- DiffusionProfiles(
  graph_dt = toy_edges,
  node_dt = toy_nodes,
  sink_types = character()
) %>%
  generate_profiles(seeds = toy_seeds, .verbose = FALSE)

round(get_profiles(plain), 3)
#>           drug_A drug_B disease_X disease_Y
#> drug_A     0.150  0.009     0.071     0.006
#> drug_B     0.005  0.150     0.009     0.048
#> P1         0.168  0.013     0.103     0.009
#> P2         0.243  0.032     0.198     0.021
#> P3         0.111  0.059     0.159     0.040
#> P4         0.044  0.101     0.083     0.068
#> P5         0.023  0.250     0.044     0.225
#> P6         0.011  0.158     0.021     0.239
#> F1         0.123  0.009     0.071     0.006
#> F2         0.036  0.041     0.069     0.028
#> F3         0.008  0.116     0.015     0.149
#> disease_X  0.075  0.018     0.150     0.012
#> disease_Y  0.003  0.045     0.006     0.150
```

Each column sums to 1. Look at `drug_A`: 0.075 of its mass sits on
`disease_X`, and from there it keeps walking into `P3` and beyond. Now
make drugs and diseases sinks and up-weight functions five-fold.

``` r

constrained <- DiffusionProfiles(
  graph_dt = toy_edges,
  node_dt = toy_nodes,
  sink_types = c("drug", "disease"),
  type_weights = c(drug = 1, protein = 1, `function` = 5, disease = 1)
) %>%
  generate_profiles(seeds = toy_seeds, .verbose = FALSE)

constrained
#> DiffusionProfiles
#>   Graph: 13 nodes | 17 edges
#>   Sink types: drug, disease 
#>   Type weights: drug = 1, protein = 1, function = 5, disease = 1 
#>   Profiles: 4 seed(s)

round(get_profiles(constrained), 3)
#>           drug_A drug_B disease_X disease_Y
#> drug_A     0.176  0.000     0.026     0.000
#> drug_B     0.000  0.167     0.001     0.021
#> P1         0.212  0.001     0.063     0.000
#> P2         0.229  0.002     0.148     0.001
#> P3         0.024  0.019     0.166     0.010
#> P4         0.011  0.037     0.089     0.019
#> P5         0.001  0.284     0.011     0.171
#> P6         0.001  0.142     0.005     0.283
#> F1         0.289  0.001     0.128     0.001
#> F2         0.023  0.038     0.180     0.020
#> F3         0.001  0.288     0.009     0.304
#> disease_X  0.031  0.002     0.173     0.001
#> disease_Y  0.000  0.017     0.001     0.169
```

`F1` goes from 0.123 to 0.289 in the `drug_A` profile, and `F3` more
than doubles in the `drug_B` one. The profiles are now dominated by the
functions each drug touches, which is exactly what you want to compare
on.

Distances between profiles are computed in Rust as well. Lower means
more similar. Ruiz et al. use the correlation distance to rank
drug-disease pairs, which is the default here; `"canberra"`, `"l1"`,
`"l2"` and `"cosine"` are also available.

``` r

calculate_profile_distances(
  constrained,
  from = c("drug_A", "drug_B"),
  to = c("disease_X", "disease_Y")
) %>%
  round(3)
#>        disease_X disease_Y
#> drug_A     0.718     1.511
#> drug_B     1.594     0.265
```

`drug_A` pairs with `disease_X`, `drug_B` with `disease_Y`. Good.

### Reactome and a MYC signature

The DuckDB has no drugs or diseases, but Reactome gives us genes,
pathways and the pathway hierarchy. We play the drug-disease game with
gene signatures instead: a node for the MYC target genes that ship with
the package, and a node for a random gene set of the same size as the
control. Both get type `signature`, and signatures are sinks.

``` r

gene_gene <- get_interactions_reactome()[, .(
  from = as.character(from),
  to = as.character(to)
)]
#> Downloading database...
#> Download complete
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpFtsWvP/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
gene_pathway <- get_gene_to_reactome()[, .(
  from = as.character(from),
  to = as.character(to)
)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpFtsWvP/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
pathway_pathway <- get_reactome_hierarchy("child_of")[, .(
  from = as.character(from),
  to = as.character(to)
)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpFtsWvP/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

all_genes <- unique(c(gene_gene$from, gene_gene$to, gene_pathway$from))

set.seed(42L)
myc_targets <- intersect(myc_genes$ensembl_gene, all_genes)
random_genes <- sample(all_genes, length(myc_targets))

signature_edges <- rbind(
  data.table(from = "MYC_targets", to = myc_targets),
  data.table(from = "random", to = random_genes)
)

nodes <- rbind(
  data.table(id = all_genes, type = "gene"),
  data.table(
    id = unique(c(gene_pathway$to, pathway_pathway$from, pathway_pathway$to)),
    type = "pathway"
  ),
  data.table(id = c("MYC_targets", "random"), type = "signature")
) %>%
  unique(by = "id")

nodes[, .N, by = type]
#>         type     N
#>       <char> <int>
#> 1:      gene 11490
#> 2:   pathway  2883
#> 3: signature     2
```

The PageRank keeps parallel edges, so an edge listed twice gets twice
the transition weight. [`unique()`](https://rdrr.io/r/base/unique.html)
won’t catch it here: Reactome lists many gene pairs as both `A-B` and
`B-A`, which is the same edge in an undirected graph.
[`dedup_edges()`](https://gregorlueg.github.io/genewalkR/reference/dedup_edges.md)
handles that and tells you where the duplicates were.

``` r

edges <- rbind(gene_gene, gene_pathway, pathway_pathway, signature_edges) %>%
  .[from != to] %>%
  dedup_edges(node_dt = nodes)
#> Removed 27112 duplicate edge(s):
#>   gene-gene: 27112

nrow(edges)
#> [1] 78268
```

For readable output, the pathway names. A few pathway ids in the
hierarchy have no name in the database; those keep their id.

``` r

pathway_labels <- get_reactome_info()[, .(
  id = as.character(reactome_id),
  label = as.character(reactome_name)
)]
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpFtsWvP/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

label_pathways <- function(ids) {
  fcoalesce(pathway_labels[match(ids, id), label], ids)
}
```

We generate a profile for both signatures and for every single pathway.
That’s close to 3,000 PageRanks over more than 14k nodes.

``` r

pathways <- nodes[type == "pathway", id]

reactome_profiles <- DiffusionProfiles(
  graph_dt = edges,
  node_dt = nodes,
  sink_types = "signature"
)

tictoc <- system.time(
  reactome_profiles <- generate_profiles(
    reactome_profiles,
    seeds = c("MYC_targets", "random", pathways)
  )
)
#> Generated 2885 diffusion profile(s) over 14375 nodes.

tictoc
#>    user  system elapsed 
#>  84.403   0.332  21.897
```

A few seconds. Blazingly fast.

#### Where does the mass go?

The obvious first look: which pathways soak up the most mass from each
signature?

``` r

profiles <- get_profiles(reactome_profiles)

top_mass <- function(seed, n = 8L) {
  mass <- sort(profiles[pathways, seed], decreasing = TRUE)[seq_len(n)]
  data.table(pathway = label_pathways(names(mass)), mass = signif(mass, 3))
}

top_mass("MYC_targets")
#>                                                          pathway    mass
#>                                                           <char>   <num>
#> 1:                                 mRNA Splicing - Major Pathway 0.00884
#> 2:                                Dengue Virus-Host Interactions 0.00574
#> 3:                                          mRNA Polyadenylation 0.00521
#> 4: Major pathway of rRNA processing in the nucleolus and cytosol 0.00482
#> 5:           Interconversion of nucleotide di- and triphosphates 0.00433
#> 6:                                      Neutrophil degranulation 0.00359
#> 7:                  rRNA modification in the nucleus and cytosol 0.00259
#> 8:                   Regulation of expression of SLITs and ROBOs 0.00250
top_mass("random")
#>                                                        pathway    mass
#>                                                         <char>   <num>
#> 1:                               Generic Transcription Pathway 0.01180
#> 2:         Expression and translocation of olfactory receptors 0.01010
#> 3:                                              Keratinization 0.00504
#> 4:                                    Neutrophil degranulation 0.00412
#> 5: Antigen processing: Ubiquitination & Proteasome degradation 0.00387
#> 6:                         Formation of the cornified envelope 0.00334
#> 7:                                 Olfactory Signaling Pathway 0.00229
#> 8:                                    Stimuli-sensing channels 0.00199
```

MYC targets land on mRNA splicing, rRNA processing and nucleotide
metabolism. MYC drives ribosome biogenesis and nucleotide synthesis, so
that checks out. The random set lands on whatever is big and well
connected: olfactory receptors, the generic transcription pathway,
neutrophil degranulation. The latter also shows up for MYC. Raw mass is
biased towards hubs.

#### Comparing profiles

This is where the approach earns its keep. Instead of asking where the
signature’s walker ends up, ask which pathway’s own profile looks most
like the signature’s.

``` r

dists <- calculate_profile_distances(
  reactome_profiles,
  from = c("MYC_targets", "random"),
  to = pathways
)

closest <- function(seed, n = 8L) {
  d <- sort(dists[seed, ])[seq_len(n)]
  data.table(pathway = label_pathways(names(d)), distance = round(d, 3))
}

closest("MYC_targets")
#>                                                      pathway distance
#>                                                       <char>    <num>
#> 1:                                   Pyrimidine biosynthesis    0.826
#> 2:                5-Phosphoribose 1-diphosphate biosynthesis    0.835
#> 3:       Interconversion of nucleotide di- and triphosphates    0.886
#> 4:                             mRNA Splicing - Major Pathway    0.886
#> 5:                                      mRNA Polyadenylation    0.896
#> 6:                                  Metabolism of polyamines    0.898
#> 7: Defective HPRT1 disrupts guanine and hypoxanthine salvage    0.902
#> 8:                              Folding of actin by CCT/TriC    0.911
closest("random")
#>                                                                     pathway
#>                                                                      <char>
#> 1: Defective SLC35A2 causes congenital disorder of glycosylation 2M (CDG2M)
#> 2:                                           Defective ABCA12 causes ARCI4B
#> 3:                                             Intestinal hexose absorption
#> 4:                                          Defective PAPSS2 causes SEMD-PA
#> 5:                            Organic anion transport by SLC22 transporters
#> 6:                                            Generic Transcription Pathway
#> 7:                                    Synthesis of UDP-N-acetyl-glucosamine
#> 8:                                             Sulfide oxidation to sulfate
#>    distance
#>       <num>
#> 1:    0.853
#> 2:    0.853
#> 3:    0.865
#> 4:    0.884
#> 5:    0.884
#> 6:    0.907
#> 7:    0.913
#> 8:    0.914
```

The hubs are gone. For MYC you get pyrimidine and PRPP biosynthesis,
nucleotide interconversion and polyamine metabolism (ODC1 is one of the
classic MYC targets), i.e. the anabolic programme MYC is known for. The
random set gives a grab bag of small pathways without a common theme.

Don’t overread the absolute numbers though. On a graph this size even
the closest hits sit between 0.8 and 0.9, and the random set’s best hit
(0.842) is barely further away than MYC’s (0.826). The ranking within a
signature carries the biology; whether a single distance is surprising
needs a null, see below.

#### Choosing a metric

Which metric then? For a rough check we need something to score against.
Take the Reactome pathways where the MYC targets are over-represented
(hypergeometric test, BH \< 0.05) and ask how high each metric ranks
them, for MYC and for the random set.

``` r

gene_to_pathway <- unique(gene_pathway)
n_genes <- uniqueN(gene_to_pathway$from)
ora <- gene_to_pathway[,
  .(k = sum(from %in% myc_targets), size = .N),
  by = .(pathway = to)
]
ora[,
  fdr := p.adjust(
    phyper(k - 1, size, n_genes - size, length(myc_targets), lower.tail = FALSE),
    method = "BH"
  )
]
myc_pathways <- intersect(ora[fdr < 0.05, pathway], pathways)
length(myc_pathways)
#> [1] 149

# AUROC of the reference pathways when sorted by distance, lower = closer
auroc <- function(d, pos) {
  is_pos <- names(d) %in% pos
  n_pos <- sum(is_pos)
  n_neg <- length(d) - n_pos
  1 - (sum(rank(d)[is_pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

metric_check <- lapply(
  c("correlation", "cosine", "l1", "l2", "canberra"),
  \(metric) {
    d <- calculate_profile_distances(
      reactome_profiles,
      from = c("MYC_targets", "random"),
      to = pathways,
      metric = metric
    )
    top_myc <- names(sort(d["MYC_targets", ]))[1:20]
    data.table(
      metric = metric,
      auroc_myc = round(auroc(d["MYC_targets", ], myc_pathways), 3),
      auroc_random = round(auroc(d["random", ], myc_pathways), 3),
      top20_myc = sum(top_myc %in% myc_pathways),
      top20_random = sum(names(sort(d["random", ]))[1:20] %in% myc_pathways),
      # parent pathways without direct member genes have no size
      top20_median_size = median(
        ora[match(top_myc, pathway), size],
        na.rm = TRUE
      )
    )
  }
) %>%
  rbindlist()

metric_check
#>         metric auroc_myc auroc_random top20_myc top20_random top20_median_size
#>         <char>     <num>        <num>     <int>        <int>             <num>
#> 1: correlation     0.974        0.723        12            0              23.5
#> 2:      cosine     0.975        0.730        12            0              29.0
#> 3:          l1     0.911        0.628        10            2              42.0
#> 4:          l2     0.853        0.765        15            4              78.0
#> 5:    canberra     0.772        0.533         5            2              11.0
median(ora[pathway %in% pathways, size])
#> [1] 11
```

Correlation and cosine are practically the same thing here: profiles sum
to 1 and are mostly near zero, so centring barely changes anything. Both
give the best separation, with none of the reference pathways in the
random set’s top 20. L1 comes next. L2 drifts towards big pathways (look
at the median size of its top 20 against the median pathway) and picks
up the most reference pathways for the random set, i.e. it rewards size
rather than specificity. Canberra weights every node’s relative
difference equally, so the long tail of near-zero entries drowns the
signal on a graph this size.

So: correlation as the default, L1 if you want a second opinion. Leave
L2 out for ranking. Canberra makes sense when you compare a profile
against something on a different scale, which is what Ruiz et al. use it
for.

Take the table with a pinch of salt:

- The reference is partly circular. The signature node hangs off exactly
  the genes the over-representation test counts, and the pathway nodes
  hang off their member genes, so this rewards metrics that pick up
  direct gene overlap.
- It’s one signature, one random draw and one graph. On a different
  graph or with type weights the order may change.
- AUROC over 2,845 pathways says nothing about whether a single distance
  is significant. That still needs a null.

### Caveats and next steps

- Profiles are dense `nodes x seeds` matrices. 2,847 seeds over 14,295
  nodes is about 325 MB of doubles, so thousands of seeds on a large
  graph get painful fast.
- There’s no permutation test yet. For a proper null you’d rewire the
  graph (the configuration model used for GeneWalk would do) and
  regenerate the profiles. Watch this space.
- Type weights are a free parameter. Ruiz et al. tuned theirs on known
  drug-disease pairs; without such a gold standard, leave them at `NULL`
  or keep them flat.
