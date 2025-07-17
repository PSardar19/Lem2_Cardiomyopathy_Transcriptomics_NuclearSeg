library(shiny)
library(dplyr)
library(DT)
library(plotly)
library(shinydashboard)

mod_dge_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # --- Inject CSS ---
    tags$style(HTML("
  /* ---------- Plot Container Styling ---------- */
  .plot-container {
    padding: 15px;
    border-radius: 10px;
    background-color: #fff;
    margin-bottom: 20px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    border: 1px solid #ddd !important;  /* ✅ Keep only outer border */
    box-sizing: border-box;
  }

  /* ✅ Remove all inner Plotly-injected borders */
  .plot-container .js-plotly-plot,
  .plot-container .plotly,
  .plot-container .main-svg,
  .plot-container .plot-container {
    border: none !important;
    box-shadow: none !important;
    background: transparent !important;
    margin: 0 !important;
    padding: 0 !important;
    outline: none !important;
    stroke: none !important;
  }

  /* ---------- Value Box Styling ---------- */
  .value-box-custom .small-box {
    height: 110px;
    border: 1px solid #ccc !important;
    box-shadow: 0 2px 6px rgba(0,0,0,0.1);
    border-radius: 10px;
    background-color: #ffffff;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    text-align: center;
    padding: 0px;
    margin: 0 10px 10px 10px;
  }

  .value-box-custom .small-box .icon-large {
    font-size: 30px;
    margin-bottom: 5px;
  }

  .value-box-custom .small-box .inner h3,
  .value-box-custom .small-box .inner p {
    margin: 0;
    padding: 0;
    font-size: 15px;
  }

  /* ---------- Layout Adjustments ---------- */
  .row.equal-height [class*='col-'] {
    display: flex;
    flex-direction: column;
  }

  .row.equal-height .value-box-custom {
    flex: 1;
  }

  .no-gutters {
    margin-left: 0 !important;
    margin-right: 0 !important;
  }

  .no-gutters > [class*='col-'] {
    padding-left: 5px !important;
    padding-right: 5px !important;
  }
  
.sidebar-panel {
  background-color: #1d2d44 !important;
  color: white !important;
  min-height: 100vh;
  padding: 20px;
  border-right: 1px solid #34495E;
}

.sidebar-panel label {
  color: white !important;
}

.sidebar-panel .form-control {
  background-color: #2a3b57 !important;
  color: white !important;
  border: 1px solid #34495E !important;
}

.sidebar-panel .form-control::placeholder {
  color: #bbb !important;
}

.sidebar-panel .selectize-input {
  background-color: #2a3b57 !important;
  color: white !important;
  border: 1px solid #34495E !important;
}

.sidebar-panel .selectize-input > input {
  color: white !important;
}

.sidebar-panel .selectize-dropdown {
  background-color: #34495E !important;
  color: white !important;
}

.sidebar-panel .irs-grid-text,
.sidebar-panel .irs-single,
.sidebar-panel .irs-bar,
.sidebar-panel .irs-from,
.sidebar-panel .irs-to {
  color: white !important;
}
  .irs-bar,
  .irs-bar-edge {
    background: #ffffff !important;
    border-color: #ffffff !important;
  }

  /* Keep text white but remove background box */
  .irs-single,
  .irs-from,
  .irs-to {
    background: transparent !important;
    color: #ffffff !important;
    border: none !important;
    box-shadow: none !important;
  }

  /* Customize track line */
  .irs-line {
    background: #99ccff !important;
    border-color: #99ccff !important;
  }

  /* Tick labels (min/max/grid) */
  .irs-grid-text,
  .irs-min,
  .irs-max {
    color: #ffffff !important;
  }

  /* Slider handle styling */
  .irs-slider {
    border: 2px solid #ffffff !important;
    background: #003366 !important;
  }
")),
    
    fluidRow(
      class = "no-gutters",
      style = "margin: 0; padding: 0;",  # remove extra space from the row
      
      # Sidebar: left panel with background and border
      column(
        width = 2,
        class = "sidebar-panel",
        style = "background-color:	#1d2d44; min-height: 100vh; padding: 20px; border-right: 1px solid #ddd;",
        selectInput(ns("pipeline"), "Select Pipeline:", choices = c("FLAMES", "NanoSeq")),
        selectInput(ns("stage"), "Select Stage:", choices = c("E18.5", "P1", "P3")),
        sliderInput(ns("padj"), "Adjusted p-value cutoff:", min = 0, max = 0.1, value = 0.05, step = 0.005),
        sliderInput(ns("logfc"), "Log2 Fold Change cutoff (absolute):", min = 0, max = 5, value = 0.263, step = 0.001)
      ),
      
      # Main content: right panel with padding only
      column(
        width = 10,
        tabsetPanel(
          tabPanel("Summary & Plots",
                   fluidRow(class = "equal-height no-gutters",
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("totalGenesBox"), width = 12))),
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("upGenesBox"), width = 12))),
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("downGenesBox"), width = 12)))
                   ),
                   fluidRow(class = "equal-height no-gutters",
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("sigGenesBox"), width = 12))),
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("topUpBox"), width = 12))),
                            column(4, div(class = "value-box-custom", valueBoxOutput(ns("topDownBox"), width = 12)))
                   ),
                   fluidRow(class = "no-gutters",
                            column(6, div(class = "plot-container", plotlyOutput(ns("volcano_plot"), height = "350px"))),
                            column(6, div(class = "plot-container", plotlyOutput(ns("ma_plot"), height = "350px")))
                   )
          ),
          tabPanel("Filtered Data Table",
                   br(),
                   downloadButton(ns("downloadFiltered"), "Export data to CSV"),
                   br(), br(),
                   DTOutput(ns("dge_table"))
          )
        )
      )
    )
  )
}

# Server function for the DGE module
mod_dge_server <- function(id, dge_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    selected_data <- reactive({
      req(input$pipeline, input$stage)
      dge_data[[input$pipeline]][[input$stage]]
    })
    
    filtered_data <- reactive({
      req(input$pipeline, input$stage)
      
      df <- selected_data() %>%
        filter(!is.na(padj)) %>%
        filter(padj <= input$padj, abs(log2FoldChange) >= input$logfc)
      
      # Round all numeric columns except padj and pvalue
      df <- df %>%
        mutate(across(where(is.numeric) & !matches("padj|pvalue"), ~ round(.x, 2))) %>%
        mutate(
          padj = formatC(padj, format = "e", digits = 2),
          pvalue = formatC(pvalue, format = "e", digits = 2)
        )
      
      # Move gene_name to the first column if it exists
      if ("gene_name" %in% colnames(df)) {
        df <- df %>% select(gene_name, everything())
      }
      
      df
    })
    
    output$dge_table <- renderDT({
      filtered_data()
    }, options = list(pageLength = 10))
    
    output$totalGenesBox <- renderValueBox({
      valueBox(
        value = nrow(selected_data()),
        subtitle = "Total Genes Tested",
        icon = icon("list")
      )
    })
    
    output$sigGenesBox <- renderValueBox({
      sig <- selected_data() %>% filter(!is.na(padj) & padj < input$padj & abs(log2FoldChange) >= input$logfc)
      valueBox(
        value = nrow(sig),
        subtitle = paste("Significant Genes (padj <=", input$padj, ")"),
        icon = icon("check-circle")
      )
    })
    
    output$upGenesBox <- renderValueBox({
      up <- filtered_data() %>% filter(log2FoldChange > input$logfc)
      valueBox(
        value = nrow(up),
        subtitle = paste("Significantly Upregulated Genes"),
        icon = icon("arrow-up")
      )
    })
    
    output$downGenesBox <- renderValueBox({
      down <- filtered_data() %>% filter(log2FoldChange < -input$logfc)
      valueBox(
        value = nrow(down),
        subtitle = paste("Significantly Downregulated Genes"),
        icon = icon("arrow-down")
      )
    })
    
    output$topUpBox <- renderValueBox({
      up <- filtered_data() %>% filter(log2FoldChange > 0)
      if (nrow(up) > 0) {
        gene <- up[which.max(up$log2FoldChange), ]
        valueBox(
          value = paste0(gene$gene_name, " (", round(gene$log2FoldChange, 2), ")"),
          subtitle = "Top Upregulated Gene",
          icon = icon("arrow-up")
        )
      }
    })
    
    output$topDownBox <- renderValueBox({
      down <- filtered_data() %>% filter(log2FoldChange < 0)
      if (nrow(down) > 0) {
        gene <- down[which.min(down$log2FoldChange), ]
        valueBox(
          value = paste0(gene$gene_name, " (", round(gene$log2FoldChange, 2), ")"),
          subtitle = "Top Downregulated Gene",
          icon = icon("arrow-down")
        )
      }
    })
    
    output$volcano_plot <- renderPlotly({
      data <- selected_data() %>% filter(!is.na(padj))
      data <- data %>% mutate(
        significance = case_when(
          padj <= input$padj & log2FoldChange >= input$logfc ~ "Up",
          padj <= input$padj & log2FoldChange <= -input$logfc ~ "Down",
          TRUE ~ "NS"
        )
      )
      
      plot_ly(
        data,
        x = ~log2FoldChange,
        y = ~-log10(padj),
        type = 'scatter',
        mode = 'markers',
        color = ~significance,
        colors = c("Down" = "blue", "NS" = "grey", "Up" = "red"),
        text = ~paste("Gene ID:", gene_id,
                      "<br>Gene Name:", gene_name,
                      "<br>Log2FC:", round(log2FoldChange, 2)),
        hoverinfo = "text",
        marker = list(size = 7, opacity = 0.6)
      ) %>%
        layout(
          title = "Volcano Plot",
          xaxis = list(title = "Log2 Fold Change"),
          yaxis = list(title = "-log10(padj)"),
          paper_bgcolor = 'rgba(0,0,0,0)',
          plot_bgcolor = 'rgba(0,0,0,0)',
          shapes = list(
            list(type = "line", x0 = input$logfc, x1 = input$logfc,
                 y0 = 0, y1 = max(-log10(data$padj), na.rm = TRUE),
                 line = list(color = "black", dash = "dot")),
            list(type = "line", x0 = -input$logfc, x1 = -input$logfc,
                 y0 = 0, y1 = max(-log10(data$padj), na.rm = TRUE),
                 line = list(color = "black", dash = "dot")),
            list(type = "line", x0 = min(data$log2FoldChange, na.rm = TRUE),
                 x1 = max(data$log2FoldChange, na.rm = TRUE),
                 y0 = -log10(input$padj), y1 = -log10(input$padj),
                 line = list(color = "black", dash = "dot"))
          )
        )
    })
    
    output$ma_plot <- renderPlotly({
      df <- selected_data()
      req(nrow(df) > 0)
      df$log2BaseMean <- log2(df$baseMean + 1)
      
      # Create 3 categories: Up, Down, Not Significant
      df <- df %>%
        mutate(
          significant = case_when(
            !is.na(padj) & padj <= input$padj & log2FoldChange >= input$logfc ~ "Up",
            !is.na(padj) & padj <= input$padj & log2FoldChange <= -input$logfc ~ "Down",
            TRUE ~ "Not Significant"
          )
        )
      
      plot_ly(
        data = df,
        x = ~log2BaseMean,
        y = ~log2FoldChange,
        type = "scatter",
        mode = "markers",
        color = ~significant,
        colors = c("Not Significant" = "grey", "Up" = "red", "Down" = "blue"),  # Added blue
        marker = list(size = 6),
        text = ~paste("Gene ID:", gene_id, "<br>Gene Name:", gene_name, "<br>log2FC:", round(log2FoldChange, 3)),
        hoverinfo = "text"
      ) %>%
        layout(
          xaxis = list(title = "Log2(Base Mean + 1)"),
          yaxis = list(title = "Log2 Fold Change"),
          title = "MA Plot",
          paper_bgcolor = 'rgba(0,0,0,0)',
          plot_bgcolor = 'rgba(0,0,0,0)'
        )
    })
    
    output$downloadFiltered <- downloadHandler(
      filename = function() {
        paste0("filtered_genes_", input$pipeline, "_", input$stage, ".csv")
      },
      content = function(file) {
        write.csv(filtered_data(), file, row.names = FALSE)
      }
    )
  })
}