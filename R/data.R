#' G2G Interactions Network
#'
#' A combined dataset of gene-to-gene interactions from OpenTargets.
#' This dataset includes any interaction present in any of the source resources:
#' IntAct, STRING, Signor, or Reactome.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{graph}{A `data.table` with columns: `from`, `to`, `from_intact`, `from_string`, `from_signor`, `from_reactome`, `edge_type`, `source`}
#'   \item{node_id}{A `data.table` with columns: `id` (ENSEMBL gene ID), `name` (gene name)}
#' }
#' @source OpenTargets
"g2g_interactions"

#' Reactome Pathway Network
#'
#' A dataset containing Reactome pathways and their relationships to other pathways and genes.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{graph}{A `data.table` of pathway relationships: `from`, `to`, `edge_type`, `source`}
#'   \item{node_id}{A `data.table` with columns: `id` (Reactome pathway ID), `name` (pathway name)}
#' }
#' @source OpenTargets
"reactome"

#' Gene Ontology Network
#'
#' A dataset containing Gene Ontology terms and their relationships.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{graph}{A `data.table` of GO term relationships: `from`, `to`, `edge_type`, `source`}
#'   \item{node_id}{A `data.table` of GO term metadata: `id` (GO term ID), `name` (GO term name), `namespace` (biological_process, molecular_function, or cellular_component)}
#' }
#' @source OpenTargets
"gene_ontology"
