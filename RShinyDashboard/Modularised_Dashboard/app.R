# Load necessary libraries
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(purrr)
library(tibble)
library(bslib)

# Load pre-processed data
# These Rdata files were ontained after running the respective scripts in DataPrep/
load("../data/dge_data2.RData")
load("../data/temporal_data.RData")
gsea_data <- readRDS("../data/gsea_full_data.rds")

# Source modules
source("../Dashboard_modules/mod_dge2.R")
source("../Dashboard_modules/mod_TemporalDGE2.R")
source("../Dashboard_modules/mod_fgsea.R")

# UI for the dashboard -> Load module UI
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  style = "background-color: #f8f9fa; padding: 20px;",
  titlePanel("Direct RNA-seq DGE Analysis Dashboard"),
  div(
    style = "background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1);",
    tabsetPanel(
      tabPanel("Stage-wise DGE", mod_dge_ui("dge")),
      tabPanel("Temporal DGE", temporalDGE_ui("temporal_dge")),
      tabPanel("GSEA Analysis", gseaModuleUI("gsea"))
    )
  )
)

server <- function(input, output, session) {
  mod_dge_server("dge", dge_data)
  temporalDGE_server("temporal_dge", temporal_dge_data, expression_data)
  gseaModuleServer("gsea", gsea_data)
}

shinyApp(ui, server)