# genewalkR package

![r_package](https://img.shields.io/badge/R_package-0.0.1.0-orange) 


## Description

This package is an attempt to implement the GeneWalk approach from 
[Ietswaart et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02264-8)
(see GitHub [here](https://github.com/churchmanlab/genewalk)) into R
leveraging Rust under the hood to generate rapidly the random walks + use the
[Burn tensor library](https://burn.dev) to fit the SkipGram model for the 
generation of the embeddings.

## Installation

### Prerequisites

1. **Rust**: Install from [rust-lang.org](https://www.rust-lang.org/tools/install)
```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. **LibTorch 2.9.0** 

Download the CPU version for your platform:

**Linux:**

```bash

cd ~
wget https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-2.9.1%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-2.9.0+cpu.zip
rm libtorch-cxx11-abi-shared-with-deps-2.9.0+cpu.zip
```

**macOS:**

```bash

cd ~
curl -L https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-2.9.0.zip -o libtorch.zip
unzip libtorch.zip
rm libtorch.zip
```

   **Windows (PowerShell):**
```powershell

cd ~
Invoke-WebRequest -Uri "https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-2.9.0%2Bcpu.zip" -OutFile "libtorch.zip"
Expand-Archive -Path libtorch.zip -DestinationPath .
Remove-Item libtorch.zip
```

3. **Set environment variable**: 

Add to your `.bashrc`, `.zshrc`, or `.Renviron`:

```bash

export LIBTORCH_PATH="$HOME/libtorch"
```
   
For R specifically, add to `~/.Renviron`:
```
LIBTORCH_PATH=/path/to/your/home/libtorch
```

4. **Install rextendr** in R:
```r
install.packages("rextendr")
```

### Install genewalkR

```r
# From GitHub
remotes::install_github("GregorLueg/genewalkR")
```

**Note**: If LibTorch is installed in a non-standard location, set 
`LIBTORCH_PATH` before installing the package.