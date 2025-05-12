# PHP2530 Final Project: Bayesian Dirichlet–Multinomial Modeling of Gut Microbiota

This repository contains the final project for the course *PHP2530: Bayesian Statistical Methods* at Brown University. The project investigates gut microbiota composition in treatment-naïve Crohn’s disease (CD) patients using a hierarchical **Bayesian Dirichlet–Multinomial regression model** which allows parameters vary with covariates, comparing results with frequentist benchmarks such as PERMANOVA.

## Overview

The human gut microbiome plays a critical role in immune regulation. Crohn's disease has been associated with microbial dysbiosis, particularly in pediatric patients. Using 16S rRNA sequencing data from the RISK cohort, we:

- Implemented a hierarchical Dirichlet–Multinomial regression in Stan.
- Adjusted for covariates (diagnosis, age, sex) at the sample level.
- Quantified differential abundance using posterior summaries.
- Validated results using frequentist approaches.

## Repository Structure

```

├── R/                  # R scripts for analysis
│   ├── QualityControl.R          # Data preprocessing and filtering
│   ├── modeling.R                # Stan model fitting using rstan
│   ├── PosteriorInference.R     # Posterior summaries and inference
│   ├── ConvergencePlot.R        # Diagnostic plots (trace, ACF)
│   ├── FreqBenchmark.R          # PERMANOVA, Wilcoxon
├── data/
│   └── RISK\_CCFA.rds            # Processed 16S data from Qiita study 550
├── stan/
│   └── dm\_regression.stan       # Stan code for Dirichlet–Multinomial regression

````

## Key Methods

### Bayesian Model
- **Model**: Dirichlet–Multinomial with log-linear regression link
- **Priors**:
  - $\beta_{pk} \sim \mathcal{N}(0, 4)$
  - $\phi_i \sim \text{Gamma}(2, 0.1)$
- **Posterior inference**: MCMC via NUTS (No-U-Turn Sampler)

### Frequentist Benchmarks
- PERMANOVA (Bray–Curtis dissimilarity)
- Wilcoxon rank-sum tests
- DESeq2 log-fold change estimation

## Results Summary

- Significant compositional shifts detected between CD and non-IBD controls.
- Taxa such as *Sutterella*, *Bacteroides*, and *Haemophilus* showed strong differential abundance.
- Frequentist methods corroborated a subset of Bayesian findings but lacked covariate adjustment or full uncertainty quantification.

## Usage

1. Clone the repository
```bash
git clone https://github.com/wjyin01/PHP2530.git
cd PHP2530
````

2. Open R and install required packages

```r
install.packages(c("rstan", "phyloseq", "vegan", "DESeq2", "tidyverse", "bayesplot"))
```

3. Run analysis scripts

Project by **Weijing Yin**, May 2025.
Course: PHP2530 - Bayesian Statistical Methods, Brown University.
