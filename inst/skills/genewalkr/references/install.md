# Installing genewalkR

## The quick route

Binaries from r-universe, no Rust needed:

```r
install.packages(
  "genewalkR",
  repos = c("https://gregorlueg.r-universe.dev", "https://cloud.r-project.org")
)
```

## From source

Source installs compile the Rust crate, so you need a toolchain:

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

```r
install.packages("rextendr")
remotes::install_github("GregorLueg/genewalkR")
```

`configure` checks that `cargo` and `rustc` are on the `PATH` and fails early if
not. A cold build compiles the full dependency graph with LTO, so it takes a
while.

## The database

The `get_*` data accessors and `GeneWalkGenerator` read a DuckDB file (~200 MB
on disk) that is downloaded from a GitHub release on first use into
`tools::R_user_dir("genewalkR", "cache")`. No network, no DB. The methods
themselves (GeneWalk on your own tables, node2vec, diffusion) do not need it.

- `download_database()` fetches it explicitly.
- `reload_db()` wipes the cache and re-fetches. Use it after a package upgrade
  that points at a new data release, or if the file is corrupt.
- `get_db_connection()` opens it read-only. Close it yourself with
  `DBI::dbDisconnect(con, shutdown = TRUE)`.

## Checking an install

Run the smoke test in `SKILL.md`. Synthetic data, no network. If all three
stages print "yes" and `get_stats()` returns rows, the Rust side compiled and
linked.

## Related packages

| Package | Reach for it when |
|---|---|
| [bixverse](https://github.com/GregorLueg/bixverse) | gene set enrichment, single cell, bulk RNAseq, PageRank-style network methods and kNN graphs. Ships its own agent skill |
| [igraph](https://r.igraph.org/) | graph construction and plain graph statistics. `DiffusionScores()` takes an igraph object |

## Keeping this skill current

The skill is versioned with the package. After upgrading genewalkR:

```r
genewalkR::install_agent_skill(overwrite = TRUE)
```

Other agents work too. `agent = "codex"` (or `"generic"` with an explicit
`dest`) writes the entry file as `AGENTS.md` with the Claude-specific
frontmatter stripped. Only Claude Code auto-discovers the skill, so elsewhere
you have to point the agent at the directory yourself.
