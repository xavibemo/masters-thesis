# MASTER'S THESIS REPOSITORY

Welcome to the repository for Xavier Benedicto's Master's Thesis! This project was undergone at the Barcelona Supercomputing Center (BSC) and the Universidad Autónoma de Madrid (UAM) during the years 2023-2024.


## Overview

This repository contains all the resources, code, data, and documentation related to my Master's thesis titled "**INVESTIGATING METABOLIC CHANGES IN GASTRIC ADENOCARCINOMA CELL LINE UNDER DIFFERENT DRUG TREATMENT: A CONSTRAINT-BASED MODELING APPROACH**". This thesis delves into the complicated metabolism shifts in phenotype in gastric cancer after different drug treatments utilizing diverse Bioinformatics methodologies, offering valuable insights into crucial metabolic alterations that may occur. 


## Table of Contents

1. [Abstract](#abstract)
2. [Scripts](#scripts)
3. [Expression Data Visualizer]()


## Abstract

Gastric carcinoma (GC) is one of the leading cause of cancer death globally. The management of advanced stages of GC necessitates a multifaceted approach, encompassing not only surgical interventions but also comprehensive multidisciplinary strategies such as combinatorial drug treatments, recently emerging as promising avenues. Recent advances in computational biology have facilitated the development of genome-scale metabolic models (GEMs), which offer a comprehensive framework for understanding the intricate metabolic dynamics. Integrating these GEMs with high-throughput omics data, GEMs provide a representation of cellular metabolism though metabolic tasks, allowing for the elucidation of complex metabolic phenotypes associated with the effects of combinatorial treatments. In this thesis, an approach to take advantage of metabolic tasks and high-throughput omics data using the previously described Tasks Inferred from Differential Expression (TIDEs) approach in conjunction with vital genes to metabolic tasks (anchor genes) is presented as ag-TIDEs. The aim is to uncover aberrant metabolic phenotypes applying the ag-TIDEs approach to try to offer a mechanistic explanation to the observed synergies in the combinatorial drug treatments of TAK1, MEK, and PI3K inhibitors on cells from the AGS cell line. Preliminary findings reveal distinct metabolic alterations induced by the drug treatments, characterized by an impairment in nucleotide synthesis and oxidative stress protection in MEK and PI3K treatments, and propose that a combination of metabolic phenotypes is key to explain the observed synergistic effects in AGS cell growth. 


## Scripts

All scripts used in the making of this master's thesis can be found at the directory `scripts/`. 


## Expression Data Visualizer 

A small web app was also developed to visualize the differential expression results gene by gene. This resource can be found by following the link https://expression-data-visualizer.shinyapps.io/shiny/ or by locally launching the app using R.

```R
library(shiny)
runApp("expression_visualizer_web_app/")
```