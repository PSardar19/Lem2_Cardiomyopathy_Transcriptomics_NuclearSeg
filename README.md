# Investigating the Role of Lem2 in Dilated Cardiomyopathy  
### Transcriptomic Analysis and Nuclear Morphometry

This repository contains the full analysis pipeline and interactive visualisation tools developed for the MSc Research Project titled:  
**"Investigating the Role of Lem2 in Dilated Cardiomyopathy: Transcriptomic Analysis and Nuclear Morphometry"**

The project aims to explore how the **Lem2 p.L13R mutation** perturbs gene expression programs and nuclear architecture in a murine model of neonatal dilated cardiomyopathy (DCM). 
It comprises:
1. Longitudinal direct RNA-seq analysis of mouse hearts (E18.5, P1, P3)
2. Machine learning-based nuclear segmentation of light-sheet microscopy images
3. Interactive data exploration via R Shiny dashboard



---

## 🧪 Key Technologies & Tools

| Category                 | Tool / Technology           | Version          |
|-------------------------|-----------------------------|------------------|
| **Transcriptomics**     | FLAMES                      | GitHub commit cf8f221               |
|                         |NanoSeq                      | v3.1.0               |
| **Differential Analysis**| DESeq2                     | v1.42.1               |
|                         | IsoformSwitchAnalyzeR       | v2.6.0               |
| **Gene Set Enrichment** | clusterProfiler             | v4.12.0               |
|                         | fgsea                     | v1.26.0                |
|                         | msigdbr                     | v2023.2                |
| **Image Analysis**      | FIJI                        | v2.16.0/1.54p                |
|                         | StarDist                    | v0.9.1                |
|                         | Tensorflow                  |v2.14.0              |
| **Interactive Visualization**| R                     | v4.4.1                |
|                         | Shiny                       | v1.10.0                |

# 🚀 Getting Started

## 1. Clone this repository

```bash
git clone https://github.com/PSardar19/Lem2_Cardiomyopathy_Transcriptomics_NuclearSeg.git
cd Lem2_Cardiomyopathy_Transcriptomics_NuclearSeg
```
## 2. Explore Each Module

### DirectRNASeqAnalysis:
Quantification pipelines using **FLAMES** and **NanoSeq**, followed by differential analyses:
- **DGE** – Differential gene expression  
- **DTE** – Differential transcript expression  
- **DTU** – Differential transcript usage  
Also includes pathway enrichment methods:  
- **GSEA** – Gene Set Enrichment Analysis  
- **ORA** – Over-Representation Analysis  

### NuclearSegmentation:
Benchmarks for 2D and 3D segmentation, including:
- Model training with **StarDist2D**  
- Morphometric analysis workflows for 2D nuclear features

### RShinyDashboard:
A modular **Shiny** application for interactive visualisation of:
- Stagewise transcriptomic effects  
- Temporal shifts in cardiac gene regulation
- Pathway Enrichment analysis (GSEA) 
Designed for intuitive exploration and accessibility.

NOTE: Each folder contains its own README with specific instructions and setup steps.

# 📂 Data Availability

Raw sequencing data, processed expression matrices, and results from Differential Expression analysis are **not included** in this repository due to:
- 📦 **Storage limitations**
- 🛡️ **Supervisory restrictions on data sharing**


