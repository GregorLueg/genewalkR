# genewalkR package

![r_package](https://img.shields.io/badge/R_package-0.0.1.0-orange) 
[![CI](https://github.com/GregorLueg/genewalkR/actions/workflows/test.yml/badge.svg)](https://github.com/GregorLueg/genewalkR/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description

This package is an attempt to implement the GeneWalk approach from 
[Ietswaart et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02264-8)
(see GitHub [here](https://github.com/churchmanlab/genewalk)) into R
leveraging Rust under the hood to generate rapidly the random walks + use the
[Burn tensor library](https://burn.dev) to fit the SkipGram model for the 
generation of the embeddings.

## Installation

You will need Rust on your system to have the package working. An installation
guide is provided [here](https://www.rust-lang.org/tools/install). There is a
bunch of further help written [here](https://extendr.github.io/rextendr/index.html)
by the rextendr guys in terms of Rust set up. (bixverse uses rextendr to interface
with Rust.)

Steps for installation:

1. In the terminal, install [Rust](https://www.rust-lang.org/tools/install)

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. In R, install [rextendr](https://extendr.github.io/rextendr/index.html):

```
install.packages("rextendr")
```

3. Finally install genewalkR:

```
devtools::install_github("https://github.com/GregorLueg/genewalkR")
```

### Install genewalkR

```r
# From GitHub
remotes::install_github("GregorLueg/genewalkR")
```

**Note**: If LibTorch is installed in a non-standard location, set 
`LIBTORCH_PATH` before installing the package.