# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

An R package over a Rust crate (extendr/rextendr). Graph methods for
computational biology: heat diffusion kernels, diffusion profiles (Ruiz
et al.), node2vec and metapath2vec embeddings, GeneWalk (Ietswaart et
al.) and embedding-drift detection. All heavy numerics live in
`src/rust`; the R layer validates, marshals and stores results.

## Build and test

``` r

Sys.setenv(DEV_BUILD = "1")   # fast iteration: keeps target dir + shared cargo registry, LTO off
devtools::load_all()          # runs configure -> cargo build -> regenerates R/extendr-wrappers.R
devtools::document()          # roxygen only; Rust wrappers come from the cargo build

tinytest::test_package("genewalkR")                    # all tests
tinytest::run_test_file("inst/tinytest/test_genewalk.R")  # single file
```

``` sh
cargo test --manifest-path src/rust/Cargo.toml   # Rust-side tests, opt-level 3
cargo clippy --manifest-path src/rust/Cargo.toml
R CMD build . && R CMD check genewalkR_*.tar.gz  # vignettes need quarto
```

Without `DEV_BUILD` every install rebuilds the full dependency graph
with LTO from a cold cargo registry. `DEBUG` is accepted but ignored:
this package is always a release build.

`configure` runs `tools/config.R`, which checks the MSRV from
`DESCRIPTION` (`SystemRequirements`) and expands `src/Makevars.in` into
`src/Makevars`. Do not edit `src/Makevars` or `src/Makevars.win`; they
are generated and gitignored.

## Architecture

### Rust and R boundary

`src/rust/src/lib.rs` holds the `extendr_module!` block and every
`#[extendr]` function. The roxygen for those functions lives in the Rust
doc comments; `src/rust/document.rs` (the `document` bin) regenerates
`R/extendr-wrappers.R` during the cargo build, so never hand-edit that
file. Rust modules:

- `graph.rs`: edge prep, CSR/CSC sparse types, SpMM, connected
  components, configuration-model rewiring for permutations.
- `diffusion.rs`: Laplacians, spectral decomposition (full or
  truncated), kernel families, and the raw/z/Monte Carlo diffusion
  methods.
- `embedding.rs`: node2vec and metapath2vec training via `node2vec-rs`.
- `data.rs`: synthetic graph generators (barbell, caveman, SBM,
  GeneWalk, differential graphs).
- `utils.rs`: cosine similarity, cross-cosine matrices, p-values, FDR.

`bixverse-rs` and `ann-search-rs` are shared crates of Gregor’s; `faer`
does the dense linear algebra. The constrained PageRank behind the
diffusion profiles comes from `bixverse_r::graph::page_rank`.

The diffusion kernel is the one stateful object:
[`rs_diffusion_kernel()`](https://gregorlueg.github.io/genewalkR/reference/rs_diffusion_kernel.md)
returns an `ExternalPtr<DiffusionKernel>` (spectral decomposition plus
transformed eigenvalues) that the R object holds and reuses across
[`diffuse_scores()`](https://gregorlueg.github.io/genewalkR/reference/diffuse_scores.md)
calls. It does not survive a session restart.

### R layer

Four S7 classes carry state through a staged pipeline, each generic
returning a new object with more properties filled in:

- `GeneWalk` (`R/genewalk.R`):
  [`GeneWalk()`](https://gregorlueg.github.io/genewalkR/reference/GeneWalk.md)
  -\>
  [`generate_initial_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_initial_emb.md)
  -\>
  [`generate_permuted_emb()`](https://gregorlueg.github.io/genewalkR/reference/generate_permuted_emb.md)
  -\>
  [`calculate_genewalk_stats()`](https://gregorlueg.github.io/genewalkR/reference/calculate_genewalk_stats.md).
- `EmbedDrift` (`R/embed_drift.R`):
  [`EmbedDrift()`](https://gregorlueg.github.io/genewalkR/reference/EmbedDrift.md)
  -\>
  [`generate_initial_embeddings()`](https://gregorlueg.github.io/genewalkR/reference/generate_initial_embeddings.md)
  -\>
  [`calculate_drift()`](https://gregorlueg.github.io/genewalkR/reference/calculate_drift.md)
  (Procrustes alignment in Rust).
- `DiffusionScores` (`R/diffusion.R`):
  [`DiffusionScores()`](https://gregorlueg.github.io/genewalkR/reference/DiffusionScores.md)
  -\>
  [`build_kernel()`](https://gregorlueg.github.io/genewalkR/reference/build_kernel.md)
  -\>
  [`diffuse_scores()`](https://gregorlueg.github.io/genewalkR/reference/diffuse_scores.md).
- `DiffusionProfiles` (`R/diffusion_profiles.R`):
  [`DiffusionProfiles()`](https://gregorlueg.github.io/genewalkR/reference/DiffusionProfiles.md)
  -\>
  [`generate_profiles()`](https://gregorlueg.github.io/genewalkR/reference/generate_profiles.md)
  -\>
  [`calculate_profile_distances()`](https://gregorlueg.github.io/genewalkR/reference/calculate_profile_distances.md).

Class definitions and their getters sit in `R/classes.R`.
[`S7::methods_register()`](https://rconsortium.github.io/S7/reference/methods_register.html)
runs in `R/zzz.R`.

`GeneWalkGenerator` (R6, `R/classes.R`) is the network factory: add
pathway and PPI sources, `$build()` pulls them from DuckDB into an S3
`DataBuilder`, then `$create_for_genes()` subsets to a gene bag and
returns a `GeneWalk`. The point is iterating over many gene bags against
one loaded network.

[`node2vec()`](https://gregorlueg.github.io/genewalkR/reference/node2vec.md)
(`R/node2vec.R`) and
[`metapath2vec()`](https://gregorlueg.github.io/genewalkR/reference/metapath2vec.md)
(`R/metapath2vec.R`) are the plain functional entry points: no class,
edge table in, embedding matrix out.

### Parameters and validation

Every tunable set goes through a `params_*()` wrapper returning a named
list, validated by a `check*`/`assert*` pair. Both are generated by
devforge: the specs live in `inst/params/params.R`, and
`devforge::forge_params()` writes `R/params-generated.R`,
`R/checkmate-params-generated.R` and `R/params-prelude-generated.R`.
Never hand-edit the generated files. Hand-written data checks (edge
tables and the like) sit in `R/checkmate_ext.R`. Rust receives the lists
and parses them into its own config structs. Adding a parameter means
touching the spec, regenerating, and the Rust-side parser.

### Data

Runtime data is a DuckDB file downloaded on first use from a GitHub
release into `tools::R_user_dir("genewalkR", "cache")`.
[`get_db_connection()`](https://gregorlueg.github.io/genewalkR/reference/get_db_connection.md)
opens it read-only,
[`reload_db()`](https://gregorlueg.github.io/genewalkR/reference/reload_db.md)
wipes and re-fetches. The `get_*` accessors in `R/data.R` are raw SQL
through
[`DBI::dbGetQuery()`](https://dbi.r-dbi.org/reference/dbGetQuery.html),
coerced to `data.table` with factor columns turned into character.

The `ingest_*`,
[`obo_parser()`](https://gregorlueg.github.io/genewalkR/reference/obo_parser.md)
and `download_*` functions in `R/data_processing.R` are not exported and
never run at analysis time. They build the DuckDB file, and
`data-raw/pull_data_to_db.R` is the script that drives them.

### Agent skill

`inst/skills/genewalkr/` is a skill that teaches coding agents to *use*
the package: `SKILL.md` (entry point, YAML frontmatter for Claude Code)
plus `references/*.md`.
[`install_agent_skill()`](https://gregorlueg.github.io/genewalkR/reference/install_agent_skill.md)
(`R/agent_skill.R`) copies it to `~/.claude/skills` or, for other
agents, writes `AGENTS.md` without the frontmatter. The directory is
lowercase because Claude skill names must be.

It ships with the package, so it goes stale with the code. When an
export, signature, default, output column or trap changes, update the
matching reference file and `references/api-index.md` in the same
change. Facts in it are checked against the code or the DB, not written
from memory.

## Conventions

- R is formatted with air (`air.toml`, format-on-save via the Posit VS
  Code extension). Vignettes are quarto (`.qmd`).
- Version bumps touch `DESCRIPTION`, `src/rust/Cargo.toml` and the
  package entry in `src/rust/Cargo.lock`. The README version badge reads
  `DESCRIPTION`. A green R-CMD-check on `main` triggers auto-tagging, so
  a bump is a release.
- The repo tracks a copier template (`.copier-answers.yml`);
  `configure`, `tools/`, `src/Makevars*.in` and the workflows come from
  it. Fix template bugs upstream where you can.
