
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MolEvolvR 📦

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The *MolEvolvR* R-package provides a framework for characterizing
proteins using molecular evolution and phylogeny.

## Installation

You can install the development version of MolEvolvR from
[GitHub](https://github.com/) with:

``` r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.19")

# Install Development Version of molevolvR
BiocManager::install("JRaviLab/MolEvolvR", auth_token = '<PersonalAccessToken>')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(molevolvR)
## basic example code
## TBA
```

## Companion repos & projects

### Related MolEvolvR codebases

- MolEvolvR web-app 1.0 \| [🔒
  Repo](https://github.com/jravilab/molevol1.0/) \| [Live
  web-app](//jravilab.org/molevolvr) \|
  [Preprint](https://doi.org/10.1101/2022.02.18.461833) \| [Case
  studies](https://jravilab.cuanschutz.edu/molevolvr/?r=&p=help)
- MolEvolvR web-app 2.0 \| [🔒
  Repo](https://github.com/jravilab/molevolvr2.0/)

### Specific use cases

- PSP evolution across the tree of life \| [🔒
  Repo](https://github.com/jravilab/psp_app/) \| [Live
  web-app](//jravilab.org/psp)
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/38809013)
- Bacterial phage defense system, avcDI \|
  [Repo](https://github.com/JRaviLab/phage_defense_avcd) \|
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/35817890/)
- DciA evolution across bacteria \|
  [Repo](https://github.com/JRaviLab/dcia_evolution) \|
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/35880876/)
- Internalins in Listeria \|
  [Repo](https://github.com/JRaviLab/inlp_listeria) \|
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/35904424/)

## Current contributors

- [David Mayer](//github.com/the-mayer) \| R-package, back-end
- [Faisal Alquaddoomi](//github.com/falquaddoomi) \| Front-end 1.0,
  back-end
- [Evan Brenner](//github.com/epbrenner) \| MolEvolvR functionality
- [Vince Rubinetti](//github.com/vincerubinetti) \| Front-end 2.0
- [Dave Bunten](//github.com/d33bs) \| Data management, engineering
- [Janani Ravi](//github.com/jananiravi) \| PI

📜 [License](https://github.com/JRaviLab/MolEvolvR/blob/main/LICENSE.md)
