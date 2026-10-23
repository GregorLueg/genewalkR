# API index

Every export, grouped. `?fn` for the signature. Generics dispatch on the class
in the section header.

## GeneWalk (`GeneWalk`, S7)

| Function | What |
|---|---|
| `GeneWalk()` | constructor from graph, gene-to-pathway table, gene and pathway IDs |
| `GeneWalkGenerator` | R6 factory: `$add_pathways()`, `$add_ppi()`, `$build()`, `$create_for_genes()`, `$reset_choices()`, `$return_full_network_dt()` |
| `check_degree_distribution()` | degree summary per edge type |
| `generate_initial_emb()` | node2vec embeddings of the real network, `n_graph` times |
| `generate_permuted_emb()` | null similarities from `n_perm` degree-preserving permutations |
| `calculate_genewalk_stats()` | p-values and FDRs per gene-pathway pair |
| `get_stats()` | results `data.table` (also on `EmbedDrift`) |
| `get_embedding()` | list of the `n_graph` embedding matrices |
| `get_graph_dt()` | the edge table |
| `get_gene_pathway_dt()` | the gene-to-pathway table |
| `plot_similarities()` | real vs permuted similarity histograms |
| `plot_gw_results()` | per-gene share of significant annotations |

## Custom GeneWalk networks (`DataBuilder`, S3)

| Function | What |
|---|---|
| `new_data_builder()` | wrap PPI, gene-to-pathway and hierarchy tables |
| `get_gw_data()` | the full combined edge table |
| `get_gw_data_filtered()` | subset to a gene list; list of `gwn`, `genes_to_pathways`, `represented_genes`, `represented_pathways` |
| `n_nodes()`, `dim()` | size of the builder |

## Embeddings

| Function | What |
|---|---|
| `node2vec()` | node2vec on an edge table, returns a matrix |
| `metapath2vec()` | metapath2vec(++) on a typed graph, returns a matrix with attributes |
| `EmbedDrift()` | constructor from two edge tables |
| `generate_initial_embeddings()` | node2vec on both graphs |
| `calculate_drift()` | Procrustes alignment and per-node cosine |
| `get_embeddings()` | one or both embeddings |
| `get_graph_dts()` | one or both edge tables |
| `get_aligned_embd()` | graph 1 aligned onto graph 2, shared nodes |

## Diffusion

| Function | What |
|---|---|
| `DiffusionScores()` | constructor from an igraph object |
| `build_kernel()` | Laplacian, eigendecomposition, kernel transform |
| `diffuse_scores()` | raw, z or Monte Carlo scores for one or more inputs |
| `get_scores()` | scores `data.table` |
| `DiffusionProfiles()` | constructor from edge and node tables, sink types, type weights |
| `generate_profiles()` | one constrained personalised PageRank per seed |
| `get_profiles()` | nodes x seeds matrix |
| `calculate_profile_distances()` | seeds x seeds distance matrix |
| `dedup_edges()` | collapse duplicate edges, aggregate weights |

## Parameter bundles

`params_genewalk()`, `params_node2vec()`, `params_metapath2vec()`,
`params_kernel()`, `params_diffusion_profiles()`, `params_genewalk_data()`.

## Data

| Function | What |
|---|---|
| `download_database()`, `reload_db()` | fetch or re-fetch the DuckDB |
| `get_db_connection()` | read-only DBI connection |
| `get_gene_info()` | Ensembl ID, symbol, biotype |
| `get_gene_ontology_info()`, `get_gene_ontology_hierarchy()`, `get_gene_to_go()` | GO terms, DAG, annotations |
| `get_reactome_info()`, `get_reactome_hierarchy()`, `get_gene_to_reactome()` | Reactome pathways, hierarchy, annotations |
| `get_interactions_combined()` | union of all interaction sources with per-source flags |
| `get_interactions_string()`, `get_interactions_intact()`, `get_interactions_pc()`, `get_interactions_reactome()`, `get_interactions_signor()` | single sources |
| `myc_genes` | Hallmark MYC targets, the example gene list |

## Synthetic data and testing

| Function | What |
|---|---|
| `synthetic_genewalk_data()` | GeneWalk network with known signal and noise genes |
| `differential_graph_test_data()` | graph pair with known differential nodes |
| `node2vec_test_data()` | barbell, caveman or SBM graphs with cluster labels |
| `evaluate_node2vec_test()` | within- vs between-cluster cosine |
| `diffusion_test_lattice()` | square lattice as igraph plus edge table |

## Agent skill

`install_agent_skill()` copies this skill to where your agent looks for it.

## Rust bindings

`rs_build_synthetic_genewalk()`, `rs_cosine_sim()`,
`rs_differential_graph_data()`, `rs_diffuse()`, `rs_diffusion_profiles()`,
`rs_gene_walk()`, `rs_gene_walk_perm()`, `rs_gene_walk_test()`,
`rs_kernel_node_names()`, `rs_metapath2vec()`, `rs_node2vec()`,
`rs_node2vec_synthetic_data()`, `rs_procrustes_align()`,
`rs_profile_distances()`. Unvalidated, interfaces may change. Use the wrappers.
