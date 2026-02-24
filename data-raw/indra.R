# pull in the INDRA statements -------------------------------------------------

gene_info <- get_gene_info()

protein_coding_genes <- gene_info[biotype == "protein_coding", symbol]


gene <- protein_coding_genes[1]

library(httr)
library(igraph)

get_statements <- function(gene, role = "agent") {
  params <- list(format = "json", agent = gene, max_stmts = 1000)
  response <- GET("https://db.indra.bio/statements/from_agents", query = params)
  stop_for_status(response)
  content(response, as = "parsed")$statements
}

extract_agents <- function(stmt) {
  agents <- c()
  if (!is.null(stmt[["subj"]])) {
    agents <- c(agents, stmt[["subj"]][["name"]])
  }
  if (!is.null(stmt[["obj"]])) {
    agents <- c(agents, stmt[["obj"]][["name"]])
  }
  if (!is.null(stmt[["members"]])) {
    agents <- c(agents, sapply(stmt[["members"]], function(m) m[["name"]]))
  }
  agents[!is.na(agents)]
}

stmts <- get_statements("NMT1")

str(stmts[[1]])

library(mirai)

length(stmts)

stmt_i <- stmts[[2]]


extract_agents(stmt)
