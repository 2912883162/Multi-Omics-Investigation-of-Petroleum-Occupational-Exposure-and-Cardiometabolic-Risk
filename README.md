# Multi-Omics Investigation of Petroleum Occupational Exposure and Cardiometabolic Risk

## About the Study

Environmental chemical exposures represent a critical yet underappreciated determinant of residual cardiometabolic risk. This project investigates the biological links connecting petroleum occupational exposure to host cardiovascular health through an integrated multi-omics framework. 

By combining population-level epidemiological profiling with deeply phenotyped metagenomic, targeted exposome, and host plasma proteomic cohorts, this work establishes low-density lipoprotein cholesterol (LDL-C) as a key intermediate phenotype and delineates an exposome-microbiome-lipid axis linking petroleum-derived chemical mixtures to selective gut microbial guild remodeling, altered endotoxigenesis, and downstream host thrombo-inflammatory alterations.

## Overview of Scripts

This repository provides the core analytical scripts and statistical modeling routines utilized across the multi-omics investigation:

- **Ecological Network Robustness**: Topological stability and attack simulations for prioritized microbial guilds (`01_network_robustness_simulation.R`).
- **Exposome Feature Selection**: Constrained elastic net regression and Bayesian kernel machine regression (BKMR) for mixture deconvolution and exposure index construction (`02_exposome_screening_and_pmds.R`).
- **Mediation Modeling**: Dose-response curve fitting via restricted cubic splines and bootstrap causal mediation analysis connecting chemical exposure, microbial dysbiosis, and host lipid levels (`03_mediation_and_rcs_analysis.R`).
- **Functional Pathway Deconstruction**: Stratified metagenomic analysis quantifying guild-level contributions to core bacterial biosynthetic pathways (`04_metagenomic_pathway_contribution.R`).
- **Host Proteomic Integration**: Weighted gene co-expression network analysis (WGCNA) and objective topological extraction of host protein-protein interaction subnetworks (`05_proteomics_wgcna_and_ppi.R`).

## Software Requirements

The analysis scripts are developed in R (>= 4.3.0) and utilize standard packages for multi-omics integration and biostatistics, including `igraph`, `glmnet`, `bkmr`, `mediation`, `rms`, and `WGCNA`.
