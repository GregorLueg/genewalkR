# genewalkR 0.2.3

- Agent skill: `install_agent_skill()` copies a skill into Claude Code, Codex
  or any other agent so it knows how to use the package. Lives in
  `inst/skills/genewalkr` and is versioned with the package.
- `node2vec()` now uses the `weight` column; before it checked for `weight` but
  read `weights`, so every run was unweighted.
- `GeneWalk` reads edge weights from `weight` instead of `weights`, in line with
  the rest of the package. A `weights` column is now ignored.
- `GeneWalk()` errors on a `graph_dt` missing any of the `interaction`,
  `part_of` or `hierarchy` edge types. It ran the check but never acted on it.
- `GeneWalkGenerator$add_pathways()` and `$add_ppi()` take several sources in
  one call. `add_pathways(c("go", "reactome"))` used to keep GO only. Defaults
  stay `"go"` and `"combined"`.
- `plot_gw_results()`: `fdr_treshold` renamed to `fdr_threshold`.
- `calculate_drift()` docs: graph 1 is aligned onto graph 2, not the reverse.
- `build_kernel()` uses the `weight` edge attribute of the igraph if present
  (non-negative, validated). Before, every edge counted as 1. Unweighted graphs
  give identical results.

# genewalkR 0.2.2

- Rust side moved onto the split bixverse-rs 0.8.0 crates, pulled in via
  `bixverse-r`. Along with it: faer 0.24, ann-search-rs 0.11 and node2vec-rs
  0.4, so the build carries a single faer again.
- `rs_procrustes_align()` returns an R error on mismatched dimensions instead
  of panicking.

# genewalkR 0.2.1

- Diffusion profiles (constrained personalised PageRank) from Ruiz et al.:
  `DiffusionProfiles()`, `generate_profiles()` and
  `calculate_profile_distances()`, plus a vignette.

# genewalkR 0.2.0

- metapath2vec wired into the package.

# genewalkR 0.1.11

## Features

- Version bumps to `ann-search-rs` to take advantage of avx2 and avx512 
  instructions

# genewalkR 0.1.10

## Features

- Rust version updates
- Nicer website and updades to `README.md` to point to R-universe for install.

# genewalkR 0.1.9

## Features

- Rust version updates

# genewalkR 0.1.8

## Features

- New data release
- Various updates to Rust crates to be up-to-date

# genewalkR 0.1.7

## Features

- Vignettes and documentations

# genewalkR 0.1.6

## Features

- Various version bumps

# genewalkR 0.1.5

Shifting versioning to `x.x.x` pattern.

## Features

- Various diffusion methods for "guilty-by-association" methods implemented.

# genewalkR 0.0.1.4

## Features

- Gene drift implemented

# genewalkR 0.0.1.3

## Features

- Added functionality to redownload the DB.
- Added Pathway Commons to the DuckDB
- Aligned the statistics of GeneWalk with the original paper and only keeping
  the largest connected component.

# genewalkR 0.0.1.2

## Features

- Moved DB from git lfs to package version.

# genewalkR 0.0.1.1

## Features

- Added plotting functions and removed unnecessary dependencies in CI/CD.

# genewalkR 0.0.1.0

## Features

- Package release.
