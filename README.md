# SASPF Distribution Analysis

This repository contains the R scripts and supplementary material used in the manuscript:

**"A New SAS-Based Power Function Distribution: Properties and Applications"**

## Description

This repository provides reproducible code for:

- Simulation studies for parameter estimation.
- Fitting the SASPF distribution and competing models to real datasets.
- Goodness-of-fit analysis using AIC, BIC, KS statistics, etc.
- Generating PDF, CDF, Survival, and Hazard function plots.
- Producing Q–Q and P–P plots.
- Constructing Total Time on Test (TTT) plots.

## Repository Structure

``` id="cz8a9l"
SASPF-Distribution-Analysis/
│
├── simulation/
│   └── simulation.r
│
├── fitting/
│   ├── fitting_1_2_data.r
│   └── fitting_third_data.r
│
├── plots/
│   ├── pdf_cdf_sf_hzf.r
│   ├── PP_QQ_BC.r
│   ├── PP_QQ_IR.r
│   ├── PP_QQ_Mye.r
│   └── TTTplot.r
│
└── README.md
