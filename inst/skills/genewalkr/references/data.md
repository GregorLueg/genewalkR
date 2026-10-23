# The bundled database

A DuckDB file fetched on first use (see `references/install.md`). Every getter
opens a read-only connection, runs one query, closes it and returns a
`data.table` with character columns. Every getter takes
`head_only = TRUE` to return five rows, which is the cheap way to see the
columns.

Genes are Ensembl IDs throughout (`ENSG...`), human only.

## Genes and annotations

| Getter | Columns | Rows |
|---|---|---|
| `get_gene_info()` | `ensembl_id`, `symbol`, `biotype` | ~79k |
| `get_gene_ontology_info()` | `go_id`, `go_name`, `name_space` | ~48k |
| `get_gene_ontology_hierarchy()` | `from`, `to`, `relationship` | ~74k |
| `get_gene_to_go()` | `from` (gene), `to` (GO term) | ~289k |
| `get_reactome_info()` | `reactome_id`, `reactome_name` | ~2.9k |
| `get_reactome_hierarchy(relationship)` | `from`, `to`, `relationship` | ~26k total |
| `get_gene_to_reactome()` | `from` (gene), `to` (pathway) | ~49k |

`get_reactome_hierarchy()` filters on `relationship`, any of `"child_of"`,
`"descendant_of"`, `"parent_of"`, `"ancestor_of"`. The default returns all four:
every direct edge twice (`child_of` and `parent_of` are the same 2.9k pairs
reversed) plus the 9.9k ancestor-descendant pairs, also twice. For an
undirected graph take one direct relationship. `GeneWalkGenerator` uses `"parent_of"`, the
metapath2vec and diffusion profile vignettes `"child_of"`.

`get_gene_ontology_hierarchy()` has no filter and returns every GO relationship:
`is_a` (~58k) plus `part_of`, `regulates`, `positively_regulates`,
`negatively_regulates`, `has_part` and a few rarer ones. Subset on
`relationship` if you only want the `is_a` tree.

## Interactions

| Getter | Extra columns | Rows | Filters |
|---|---|---|---|
| `get_interactions_combined()` | logical `intact`, `pc`, `reactome`, `signor`, `string` per source | ~2.0M | none |
| `get_interactions_string(threshold)` | `string_score`, `string_high_quality` | ~13M | `string_score >= threshold` |
| `get_interactions_intact(threshold, physical_interactions)` | `intact_score`, `physical_association` | ~1.2M | score threshold, physical only |
| `get_interactions_pc()` | `type` | ~1.8M | none |
| `get_interactions_reactome()` | `interaction_type` | ~54k | none |
| `get_interactions_signor()` | `interaction_type` | ~38k | none |

`pc` is PathwayCommons. Pull full STRING only if you mean it; use a threshold.

## Bundled data

`data(myc_genes)`: the 200 MSigDB Hallmark MYC targets with `ensembl_gene`,
`gene_symbol`, `gene_set`. The vignettes use it as the example gene list.

## Raw SQL

`get_db_connection()` hands you the connection. Tables: `gene_main`,
`gene_ontology_main`, `gene_ontology_dag`, `gene_ontology_genes`,
`reactome_main`, `reactome_ontology`, `reactome_genes`,
`interactions_combined`, `interactions_string`, `interactions_intact`,
`interactions_pc`, `interactions_reactome`, `interactions_signor`.

```r
con <- get_db_connection()
on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
DBI::dbGetQuery(
  con,
  "SELECT * FROM interactions_string WHERE string_score >= 0.9"
)
```

Raw SQL gives you a plain `data.frame`; `data.table::setDT()` it if you want
to carry on in data.table.
