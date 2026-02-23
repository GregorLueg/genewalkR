# Package index

## Core workflow

Main classes and methods to run GeneWalk

- [`GeneWalk()`](https://gregorlueg.github.io/genewalkR/reference/GeneWalk.md)
  : GeneWalk
- [`GeneWalkGenerator`](https://gregorlueg.github.io/genewalkR/reference/GeneWalkGenerator.md)
  : Gene Walk Network Generator
- [`generate_initial_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_initial_emb.md)
  : Generate the initial graph embedding
- [`generate_permuted_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_permuted_emb.md)
  : Generate the random, permuted embedding
- [`calculate_genewalk_stats()`](https://gregorlueg.github.io/genewalkR/reference/calculate_genewalk_stats.md)
  : Calculate the pathway and gene statistics
- [`plot_gw_results()`](https://gregorlueg.github.io/genewalkR/reference/plot_gw_results.md)
  : Plot the gene walk results

## Node2Vec

Run Node2Vec

- [`node2vec()`](https://gregorlueg.github.io/genewalkR/reference/node2vec.md)
  : Node2Vec

## Helpers and getters

Various getters and helper classes

- [`get_stats()`](https://gregorlueg.github.io/genewalkR/reference/get_stats.md)
  : Get the statistical results
- [`get_graph_dt()`](https://gregorlueg.github.io/genewalkR/reference/get_graph_dt.md)
  : Get the edge data.table
- [`get_gw_data()`](https://gregorlueg.github.io/genewalkR/reference/get_gw_data.md)
  : Get the full GW data stored in the class
- [`get_gw_data_filtered()`](https://gregorlueg.github.io/genewalkR/reference/get_gw_data_filtered.md)
  : Get the filtered GW data for a list of genes
- [`get_gene_pathway_dt()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_pathway_dt.md)
  : Get the gene to pathway data.table
- [`get_embedding()`](https://gregorlueg.github.io/genewalkR/reference/get_embedding.md)
  : Get the embedding
- [`get_gene_info()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_info.md)
  : Get the gene main information
- [`n_nodes()`](https://gregorlueg.github.io/genewalkR/reference/n_nodes.md)
  : Number of nodes represented
- [`new_data_builder()`](https://gregorlueg.github.io/genewalkR/reference/new_data_builder.md)
  : Generate a new DataBuilder helper class

## Data sources

How to download and interact with the internal DB.

- [`download_database()`](https://gregorlueg.github.io/genewalkR/reference/download_database.md)
  : Download database from GitHub release
- [`get_db_connection()`](https://gregorlueg.github.io/genewalkR/reference/get_db_connection.md)
  : Internal connector to the DB
- [`get_interactions_combined()`](https://gregorlueg.github.io/genewalkR/reference/get_interactions_combined.md)
  : Get the Combined network
- [`get_interactions_intact()`](https://gregorlueg.github.io/genewalkR/reference/get_interactions_intact.md)
  : Get the Intact network
- [`get_interactions_reactome()`](https://gregorlueg.github.io/genewalkR/reference/get_interactions_reactome.md)
  : Get the Reactome network
- [`get_interactions_signor()`](https://gregorlueg.github.io/genewalkR/reference/get_interactions_signor.md)
  : Get the SIGNOR network
- [`get_interactions_string()`](https://gregorlueg.github.io/genewalkR/reference/get_interactions_string.md)
  : Get the STRING network
- [`get_gene_ontology_info()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_ontology_info.md)
  : Get the Gene Ontology main information
- [`get_gene_ontology_hierarchy()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_ontology_hierarchy.md)
  : Get the Gene Ontology
- [`get_reactome_info()`](https://gregorlueg.github.io/genewalkR/reference/get_reactome_info.md)
  : Get the Reactome main information
- [`get_reactome_hierarchy()`](https://gregorlueg.github.io/genewalkR/reference/get_reactome_hierarchy.md)
  : Get the Reactome pathway hierarchy
- [`get_gene_to_go()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_to_go.md)
  : Get the Gene Ontology to gene connections
- [`get_gene_to_reactome()`](https://gregorlueg.github.io/genewalkR/reference/get_gene_to_reactome.md)
  : Get the Reactome to gene connections

## (Synthetic) Data & testing

- [`synthetic_genewalk_data()`](https://gregorlueg.github.io/genewalkR/reference/synthetic_genewalk_data.md)
  : Generate synthetic GeneWalk data
- [`node2vec_test_data()`](https://gregorlueg.github.io/genewalkR/reference/node2vec_test_data.md)
  : Generate Simple Node2Vec Test Data
- [`evaluate_node2vec_test()`](https://gregorlueg.github.io/genewalkR/reference/evaluate_node2vec_test.md)
  : Evaluate node2vec test results
- [`myc_genes`](https://gregorlueg.github.io/genewalkR/reference/myc_genes.md)
  : MYC target genes

## Params

Wrapper functions around core parameters

- [`params_genewalk_data()`](https://gregorlueg.github.io/genewalkR/reference/params_genewalk_data.md)
  : GeneWalk synthetic data params
- [`params_node2vec()`](https://gregorlueg.github.io/genewalkR/reference/params_node2vec.md)
  : Wrapper function for the node2vec parameters

## Rust wrappers

Everything rusty - only use this if you know what you are doing…

- [`rs_build_synthetic_genewalk()`](https://gregorlueg.github.io/genewalkR/reference/rs_build_synthetic_genewalk.md)
  : Generate synthetic GeneWalk data with controlled signal structure
- [`rs_cosine_sim()`](https://gregorlueg.github.io/genewalkR/reference/rs_cosine_sim.md)
  : Cosine similarity between two vectors
- [`rs_gene_walk()`](https://gregorlueg.github.io/genewalkR/reference/rs_gene_walk.md)
  : Generate GeneWalk node embeddings
- [`rs_gene_walk_perm()`](https://gregorlueg.github.io/genewalkR/reference/rs_gene_walk_perm.md)
  : Generate permuted embeddings for null distribution
- [`rs_gene_walk_test()`](https://gregorlueg.github.io/genewalkR/reference/rs_gene_walk_test.md)
  : Calculate the test statistics
- [`rs_generate_pathway_structure()`](https://gregorlueg.github.io/genewalkR/reference/rs_generate_pathway_structure.md)
  : Generate pathway structure and gene-pathway associations
- [`rs_node2vec_synthetic_data()`](https://gregorlueg.github.io/genewalkR/reference/rs_node2vec_synthetic_data.md)
  : Generate synthetic data for node2vec

## Internal

Internal helper functions

- [`assertGeneWalkDataParams()`](https://gregorlueg.github.io/genewalkR/reference/assertGeneWalkDataParams.md)
  : Assert GeneWalk data parameters
- [`assertGeneWalkDataTable()`](https://gregorlueg.github.io/genewalkR/reference/assertGeneWalkDataTable.md)
  : Assert GeneWalk input data.table
- [`assertGeneWalkGraphDt()`](https://gregorlueg.github.io/genewalkR/reference/assertGeneWalkGraphDt.md)
  : Assert GeneWalk graph data.table
- [`assertNode2VecParam()`](https://gregorlueg.github.io/genewalkR/reference/assertNode2VecParam.md)
  : Assert node2vec parameters
- [`checkGeneWalkDataParams()`](https://gregorlueg.github.io/genewalkR/reference/checkGeneWalkDataParams.md)
  : Check GeneWalk data parameters
- [`checkGeneWalkDataTable()`](https://gregorlueg.github.io/genewalkR/reference/checkGeneWalkDataTable.md)
  : Check GeneWalk input data.table
- [`checkGeneWalkGraphDt()`](https://gregorlueg.github.io/genewalkR/reference/checkGeneWalkGraphDt.md)
  : Check GeneWalk graph data.table
- [`checkNode2VecParams()`](https://gregorlueg.github.io/genewalkR/reference/checkNode2VecParams.md)
  : Check node2vec parameters
- [`gene_ontology_source`](https://gregorlueg.github.io/genewalkR/reference/gene_ontology_source.md)
  : Gene Ontology Data Source
- [`opentargets_source`](https://gregorlueg.github.io/genewalkR/reference/opentargets_source.md)
  : OpenTargets Data Source
