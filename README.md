This repository contains R scripts used to generate the main figures for the manuscript:

### **Transcriptional programs of cell identity and p53-induced stress responses are associated with distinctive features of spatial genome organization**  
Gony Shanel, Tsung-Han S. Hsieh, Claudia Cattoglio, Hadar Amira Haham, Hsin-Jung Chou, Jack Z. Li, Ron Shamir, Xavier Darzacq, and Ran Elkon.


---

## Overview

This repository provides original R scripts developed for specific analyses in Figures 2–6 of the manuscript.  
The scripts are organized by figure and cover tasks such as A/B compartment assignment, TAD boundary analysis, chromatin loop quantification, enhancer-promoter loop identification, and integration of gene expression and chromatin architecture in response to p53 activation.

Each script operates on example input files provided in the `data/` folder and can be run independently.  
Supporting functions are collected in `_functions.r` and `_functions_hic.r`.

> **Note:** This repository contains only the custom R scripts created for the study.  
> Operation preformed with external command-line tools (deeptools, FAN-C, ect.) are not included here.

Scripts are organized into folders by figure number:
- `Figure_2/`: A/B compartment analysis
- `Figure_3/`: TAD detection and characterization
- `Figure_4/`: Chromatin loop detection and quantification
- `Figure_5/`: Identification of functional enhancer-promoter loops
- `Figure_6/`: Analysis of p53-induced transcriptional responses

Supporting R functions used across analyses are provided in:
- `_functions.r`
- `_functions_hic.r`

Small example data files needed to run the scripts are provided in the `data/` folder.

---

## Requirements

**R packages**:
- data.table
- ggplot2
- bedtoolsr
- parallel
- DESeq2

**External tools**:
- bedtools

---

## Usage

1. Example input files are placed in the `data/` folder.
2. Scripts can be run independently within each `Figure_X/` folder.
3. Paths to data files must be manually set at the beginning of each script if needed.
4. Supporting functions are sourced automatically from `_functions.r` and `_functions_hic.r`.

> ⚠️ These scripts are designed for reproducibility and demonstration purposes using small example files, not for full-scale reprocessing of raw datasets.

---

## Citation

If you use these scripts, please cite:

Gony Shanel et al.,  
**"Transcriptional programs of cell identity and p53-induced stress responses are associated with distinctive features of spatial genome organization."** (2025)

---
