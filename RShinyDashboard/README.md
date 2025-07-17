# RShiny Dashboard for Transcriptomic Exploration

This interactive dashboard provides a streamlined interface to explore gene expression and pathway enrichment dynamics in **Lem2-mediated Dilated Cardiomyopathy**. Built using **modularised RShiny components**, the dashboard is scalable, maintainable, and easy to extend.

---

## 🧩 Modules Included

The dashboard is divided into distinct, self-contained modules — each designed to visualise a key aspect of the transcriptomic analysis. All modules support outputs from both **NanoSeq** and **FLAMES** pipelines.

### 🔹 Stagewise DGE Module
- Summary statistics with interactive volcano plots  
- Filterable & downloadable DGE result tables  
- Comparison across all developmental stages and both pipelines  

### 🔹 Temporal DGE Module
- Volcano plots & MA plots for developmental transitions  
- Temporal expression trajectories for selected genes  
- Interactive filtering and export of results  

### 🔹 GSEA Module
- Gene Set Enrichment Analysis (GO BP, Hallmark)  
- Stagewise and temporal enrichment results  
- Cross-condition comparisons (WT vs Mut)  
- Visualisation of pathway enrichment across stages and transitions

---

## 📥 Data Input & Flow

- The dashboard visualises results generated from the **DirectRNASeqAnalysis/Analysis/** directory.
- Before being loaded into the dashboard, all analysis outputs are **restructured into dashboard-friendly formats** using scripts provided in the `DataPrep/` folder.

---

# 🖥️ Running the Dashboard

### 1️⃣ Install Required Packages

To run the dashboard locally, the following R packages must be installed:

| Package         | Version   |
|----------------|-----------|
| bslib          | 0.9.0     |
| tibble         | 3.2.1     |
| purrr          | 1.0.4     |
| plotly         | 4.10.4    |
| ggplot2        | 3.5.2     |
| DT             | 0.33      |
| tidyr          | 1.3.1     |
| dplyr          | 1.1.4     |
| shinydashboard | 0.7.3     |
| shiny          | 1.10.0    |

### 2️⃣ Launch the Dashboard
In R or RStudio, run the following:

```bash
shiny::runApp("RShinyDashboard/Modularised_Dashboard")
```
This will open the dashboard in the web browser with all modules loaded.






  
