# GeneWalk

Ietswaart et al., Genome Biol 2021. Given a gene list, which of each gene's
annotated GO terms or Reactome pathways are relevant in this context? The
network is genes plus pathways: gene-gene interactions among your genes, gene to
pathway annotations, and the pathway hierarchy. node2vec embeds it; a gene and
one of its annotated pathways count as relevant when their cosine similarity
beats a null built from degree-preserving permuted networks.

It answers a different question from enrichment. Enrichment ranks pathways for
the list; GeneWalk ranks, per gene, the annotations that fit the list.

## Building the network from the DB

`GeneWalkGenerator` loads the full network once, then subsets it per gene list.
Build once, call `$create_for_genes()` as often as you like.

```r
gen <- GeneWalkGenerator$new()
gen$add_pathways("go", go_namespace = "biological_process")
gen$add_ppi("combined")
gen$build()

data(myc_genes)
gw <- gen$create_for_genes(genes = myc_genes$ensembl_gene)
gw
check_degree_distribution(gw)
```

- `$add_pathways()` and `$add_ppi()` take one or several sources and
  accumulate across calls. Unknown sources error.
- Pathway sources: `"go"` (default), `"reactome"`. GO namespaces default to all
  three; pass `go_namespace` to restrict.
- PPI sources: `"combined"` (default), `"string"`, `"signor"`, `"reactome"`,
  `"intact"`.
  `string_threshold` and `intact_threshold` (0 to 1) and `intact_physical_only`
  apply to their source. Unthresholded STRING is ~13M edges; set a threshold.
- `$build()` errors without at least one pathway and one PPI source. The whole
  pathway hierarchy goes into every gene network, not just the part your genes
  touch, so expect tens of thousands of `hierarchy` edges with full GO BP.
- `$create_for_genes()` warns about genes with no edges and drops them. A few
  dropped is normal (17 of the 200 MYC targets on Reactome); all of them means
  you passed symbols instead of Ensembl IDs.
- `$reset_choices()` clears sources and the built network.
  `$return_full_network_dt()` gives the whole thing.

`check_degree_distribution()` prints degree summaries per edge type. Few or no
`interaction` edges means your genes barely talk to each other in the PPI and
GeneWalk has little to work with. A random gene set shows exactly that.

## Custom networks

Bring your own tables via the S3 `DataBuilder`:

```r
db <- new_data_builder(
  ppis = get_interactions_reactome(),
  gene_to_pathways = get_gene_to_reactome(),
  pathway_hierarchy = get_reactome_hierarchy(relationship = "parent_of")
)
sub <- get_gw_data_filtered(db, gene_ids = my_genes)

gw <- with(sub, GeneWalk(
  graph_dt = gwn,
  gene_to_pathway_dt = genes_to_pathways,
  gene_ids = represented_genes,
  pathway_ids = represented_pathways
))
```

Or call `GeneWalk()` directly on any `graph_dt` with `from`, `to` and `type`
(`"interaction"`, `"part_of"`, `"hierarchy"`). All three types have to be there.

## Running it

```r
gw <- generate_initial_emb(gw, embd_dim = 8L, n_graph = 3L,
  genewalk_params = params_genewalk())
gw <- generate_permuted_emb(gw, n_perm = 3L)
gw <- calculate_genewalk_stats(gw)
res <- get_stats(gw)
```

- `n_graph` independent embeddings of the real network; `n_perm` permuted
  networks for the null. 3 and 3 follow the paper.
- `generate_permuted_emb()` reuses the node2vec parameters, `embd_dim` and
  `directed` stored by `generate_initial_emb()`. To change them, rerun from the
  start.
- The vignette drops to `params_genewalk(walks_per_node = 25L)` on real
  networks for speed (the paper uses 100). `num_workers > 1L` is much faster
  but no longer reproducible.

## Reading the stats

`get_stats()` returns a `data.table`, one row per annotated gene-pathway pair
(unannotated pairs are not tested), sorted by `similarity` descending:

| Column | Meaning |
|---|---|
| `gene`, `pathway` | IDs as in the network |
| `similarity`, `sem_sim` | mean cosine across the `n_graph` embeddings, and its standard error |
| `avg_pval`, `pval_ci_lower`, `pval_ci_upper` | p-value against the pooled permutation null, aggregated across the `n_graph` replicates, with 95% CI |
| `avg_global_fdr`, `global_fdr_*` | BH across all tested pairs |
| `avg_gene_fdr`, `gene_fdr_*` | BH within each gene's annotations |

`avg_gene_fdr` is the GeneWalk headline: per gene, which annotations hold up.
`avg_global_fdr` corrects across every tested pair instead.

IDs are Ensembl and GO/Reactome IDs. Translate for display:

```r
go <- get_gene_ontology_info()
genes <- get_gene_info()
res[, `:=`(
  symbol = genes$symbol[match(gene, genes$ensembl_id)],
  term = go$go_name[match(pathway, go$go_id)]
)]
```

## Plots

- `plot_similarities(gw)`: cosine similarity of the real pairs next to each
  permutation null. A real distribution that looks like the null means no
  signal.
- `plot_gw_results(gw, fdr_threshold = 0.25)`: per gene, the share of
  annotations passing `avg_gene_fdr` against the number of annotations.

## Synthetic data

`synthetic_genewalk_data()` builds a network with signal genes (annotations
inside one ontology subtree, PPI among themselves) and degree-matched noise
genes (annotations spread across subtrees). Signal gene IDs contain `"anchor"`,
so `grepl("anchor", res$gene)` is the ground truth.
