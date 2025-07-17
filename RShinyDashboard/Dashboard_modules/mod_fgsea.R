# GSEA Shiny Module - Adapted for our data structure
library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(fgsea) 
library(stringr)

# UI Module
gseaModuleUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("GSEA Analysis Dashboard"),
    
    # Sidebar with filters
    sidebarLayout(
      sidebarPanel(
        h4("Filters"),
        
        # Pipeline selection
        selectInput(ns("pipeline"), 
                    "Pipeline:",
                    choices = c("FLAMES", "NanoSeq"),
                    selected = "FLAMES"),
        
        selectInput(ns("analysis_type"), 
                    "Analysis Type:",
                    choices = c("Stagewise", "Temporal"),
                    selected = "Stagewise"),
        
        # Stage selection
        conditionalPanel(
          condition = paste0("input['", ns("analysis_type"), "'] == 'Stagewise'"),
          selectInput(ns("stage"), 
                      "Stage:",
                      choices = c("E18.5", "P1", "P3"),
                      selected = "E18.5")
        ),
        
        # Transition selection for Temporal
        conditionalPanel(
          condition = paste0("input['", ns("analysis_type"), "'] == 'Temporal'"),
          selectInput(ns("transition"), 
                      "Transition:",
                      choices = c("E18.5 -> P1", "P1 -> P3", "E18.5 -> P3"))
        ),
        
        # Condition selection for Temporal
        conditionalPanel(
          condition = paste0("input['", ns("analysis_type"), "'] == 'Temporal'"),
          selectInput(ns("condition"), 
                      "Condition:",
                      choices = c("WT", "Mut"),
                      selected = "WT")
        ),
        
        # Collection selection
        selectInput(ns("collection"), 
                    "Pathway Database:",
                    choices = c("All", "HALLMARK", "KEGG", "REACTOME", "GO_BP", "GO_MF", "GO_CC"),
                    selected = "HALLMARK"),
        
        # Significance level
        selectInput(ns("significance"), 
                    "Significance Level:",
                    choices = c("All", "Highly Significant", "Significant", "Not Significant"),
                    selected = "Significant"),
        
        # Direction filter
        selectInput(ns("direction"), 
                    "Direction:",
                    choices = c("All", "Activated", "Suppressed"),
                    selected = "All"),
        
        # Selecting top N pathways
        numericInput(ns("top_n"), 
                     label = "Number of Top Pathways:",
                     value = 5,
                     min = 1,
                     max = 20),
        
        # P-value threshold
        numericInput(ns("pvalue_threshold"), 
                     "Adjusted P-value Threshold:",
                     value = 0.05,
                     min = 0,
                     max = 1,
                     step = 0.01),
        
        hr(),
        h5("Data Summary"),
        verbatimTextOutput(ns("data_summary")),
      ),
      
      # Main panel with tabs
      mainPanel(
        tabsetPanel(
          id = ns("main_tabs"),
          
          # Results Tab
          tabPanel(
            "Pathway Results",
            fluidRow(
              column(12,
                     h4("Enrichment Analysis Results"),
                     br()
              )
            ),
            
            # Results Table
            fluidRow(
              column(12,
                     h5("Filtered Results"),
                     DTOutput(ns("results_table"))
              )
            ),
            
            br(),
            
            # Plots
            fluidRow(
              column(12,
                     h5("Pathway Enrichment Dot Plot"),
                     plotlyOutput(ns("dot_plot"), height = "400px")
              )
            ),
            
            br(),
            
            fluidRow(
              column(12,
                     h5("Top Pathways by NES"),
                     plotOutput(ns("top_pathways_plot"), height = "500px")
              )
            ),
            
            br(),
            
            fluidRow(
              column(8,  # Adjust width as needed
                     h5("Enrichment Plot for Selected Pathway"),
                     selectInput(ns("selected_pathway"), 
                                 "Select Pathway:",
                                 choices = NULL),
                     plotOutput(ns("gsea_plot"), height = "400px")
              ),
              column(4,  # Adjust width as needed
                     h5("Leading Edge Genes in Selected Pathway"),
                     uiOutput(ns("gene_list"))
              )
            )
          )
        )
      )
    )
  )
}

gseaModuleServer <- function(id, gsea_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for debugging
    debug_info <- reactiveValues(
      filtered_count = 0,
      gene_key = "",
      ranked_key = ""
    )
    
    # FIXED filtered_data function with condition handling
    filtered_data <- reactive({
      req(gsea_data)
      data <- gsea_data$gsea_results
      
      # Create required columns
      data <- data %>%
        mutate(
          Direction = ifelse(NES > 0, "Activated", "Suppressed"),
          Significance = case_when(
            padj < 0.01 ~ "Highly Significant",
            padj < 0.05 ~ "Significant",
            TRUE ~ "Not Significant"
          )
        )
      
      # Start with all data and track each filtering step
      filtered <- data
      step_results <- list()
      step_results[["original"]] <- nrow(filtered)
      
      # Apply filters one by one and track results
      if (!is.null(input$pipeline) && input$pipeline != "All") {
        filtered <- filtered %>% filter(pipeline == input$pipeline)
        step_results[["after_pipeline"]] <- nrow(filtered)
      }
      
      if (!is.null(input$analysis_type) && input$analysis_type != "All") {
        filtered <- filtered %>% filter(analysis_type == input$analysis_type)
        step_results[["after_analysis_type"]] <- nrow(filtered)
      }
      
      # Special handling for stage/transition based on analysis type
      if (input$analysis_type == "Stagewise" && !is.null(input$stage) && input$stage != "All") {
        filtered <- filtered %>% filter(stage == input$stage)
        step_results[["after_stage"]] <- nrow(filtered)
      }
      
      if (input$analysis_type == "Temporal" && !is.null(input$transition) && input$transition != "All") {
        filtered <- filtered %>% filter(transition == input$transition)
        step_results[["after_transition"]] <- nrow(filtered)
      }
      
      # CONDITION FILTER FIX - handle different value representations
      if (input$analysis_type == "Temporal" && !is.null(input$condition)) {
        # Map UI condition values to data values
        condition_value <- case_when(
          input$condition == "WT" ~ "WT",
          input$condition == "Mut" ~ "Mut",
          TRUE ~ input$condition
        )
        
        filtered <- filtered %>% 
          filter(condition == condition_value | 
                   condition == paste0(condition_value, "_vs_WT") |
                   condition == paste0("WT_vs_", condition_value))
        
        step_results[["after_condition"]] <- nrow(filtered)
      }
      
      if (!is.null(input$collection) && input$collection != "All") {
        filtered <- filtered %>% filter(collection == input$collection)
        step_results[["after_collection"]] <- nrow(filtered)
      }
      
      if (!is.null(input$significance) && input$significance != "All") {
        if (input$significance == "Significant") {
          # Include both "Significant" and "Highly Significant"
          filtered <- filtered %>% filter(Significance %in% c("Significant", "Highly Significant"))
        } else {
          filtered <- filtered %>% filter(Significance == input$significance)
        }
        step_results[["after_significance"]] <- nrow(filtered)
      }
      
      if (!is.null(input$direction) && input$direction != "All") {
        filtered <- filtered %>% filter(Direction == input$direction)
        step_results[["after_direction"]] <- nrow(filtered)
      }
      
      if (!is.null(input$pvalue_threshold) && !is.na(input$pvalue_threshold)) {
        filtered <- filtered %>% filter(padj <= input$pvalue_threshold)
        step_results[["after_pvalue"]] <- nrow(filtered)
      }
      
      # Store step results for debugging
      attr(filtered, "step_results") <- step_results
      debug_info$filtered_count <- nrow(filtered)  # UPDATE DEBUG COUNT
      return(as.data.frame(filtered))  # Convert to data.frame for better compatibility
    })
    
    # TEMPORARY DIAGNOSTIC - Keep for debugging
    output$diagnostic <- renderText({
      req(gsea_data)
      data <- gsea_data$gsea_results
      filtered <- filtered_data()
      step_results <- attr(filtered, "step_results")
      
      # Show step-by-step filtering
      step_text <- ""
      if (!is.null(step_results)) {
        for (step_name in names(step_results)) {
          step_text <- paste0(step_text, step_name, ": ", step_results[[step_name]], " rows\n")
        }
      }
      
      # Get current condition value
      cond_value <- if (input$analysis_type == "Temporal") input$condition else "N/A"
      
      paste0(
        "=== DATA INSPECTION ===\n",
        "Original pathways: ", nrow(data), "\n",
        "Unique conditions: ", paste(unique(data$condition), collapse = ", "), "\n",
        "Current UI condition: ", cond_value, "\n",
        "=== FILTERING STEPS ===\n",
        step_text,
        "Filtered count: ", nrow(filtered), "\n",
        "=== CURRENT FILTER VALUES ===\n",
        "input$pipeline: '", input$pipeline, "'\n",
        "input$analysis_type: '", input$analysis_type, "'\n",
        "input$stage: '", ifelse(!is.null(input$stage), input$stage, "N/A"), "'\n",
        "input$transition: '", ifelse(!is.null(input$transition), input$transition, "N/A"), "'\n",
        "input$condition: '", cond_value, "'\n",
        "input$collection: '", input$collection, "'\n",
        "input$significance: '", input$significance, "'\n",
        "input$direction: '", input$direction, "'\n",
        "input$pvalue_threshold: ", input$pvalue_threshold, "\n"
      )
    })
    
    # Update pathway choices
    observe({
      data <- filtered_data()
      if (nrow(data) > 0) {
        choices <- setNames(
          paste(data$collection, data$pathway, sep = "::"),
          paste0(data$pathway, " [", data$collection, "]")
        )
        updateSelectInput(session, "selected_pathway", choices = choices)
      } else {
        updateSelectInput(session, "selected_pathway", choices = NULL)
      }
    })
    
    # Data summary
    output$data_summary <- renderText({
      data <- filtered_data()
      n_pathways <- nrow(data)
      n_activated <- sum(data$Direction == "Activated")
      n_suppressed <- sum(data$Direction == "Suppressed")
      
      paste0(
        "Pathways: ", n_pathways, "\n",
        "Activated: ", n_activated, "\n",
        "Suppressed: ", n_suppressed
      )
    })
    
    # Results Table - SIMPLIFIED FOR DEBUGGING
    output$results_table <- renderDT({
      data <- filtered_data()
      
      # If no data, return a message
      if (is.null(data) || nrow(data) == 0) {
        return(datatable(
          data.frame(Message = "No pathways match the current filters"),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }
      
      # Simplified version with key columns only
      display_data <- data %>%
        select(
          Pathway = pathway,
          Collection = collection,
          NES,
          `Adj P-value` = padj,
          Direction
        ) %>%
        arrange(`Adj P-value`)
      
      datatable(
        display_data,
        options = list(
          scrollX = TRUE, 
          pageLength = 10,
          dom = 'ft'
        ),
        rownames = FALSE,
        selection = 'single'
      ) %>%
        formatRound(c("NES", "Adj P-value"), 4)
    })
    
    # Dot Plot - USING PATHWAY ID
    output$dot_plot <- renderPlotly({
      data <- filtered_data()
      if (nrow(data) == 0) return(NULL)
      
      # Create hover text with both ID and description
      data$hover_text <- paste(
        "<b>", data$pathway, "</b><br>",
        "Description: ", data$pathway, "<br>",
        "Database: ", data$collection, "<br>",
        "NES: ", round(data$NES, 3), "<br>",
        "P.adj: ", formatC(data$padj, format = "e", digits = 2)
      )
      
      p <- ggplot(data, aes(x = NES, y = -log10(padj), 
                            size = size, 
                            color = Direction,
                            text = hover_text)) +
        geom_point(alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
        scale_color_manual(values = c("Activated" = "red", "Suppressed" = "blue")) +
        labs(x = "Normalized Enrichment Score (NES)", 
             y = "-log10(Adjusted P-value)",
             size = "Gene Count") +
        theme_bw() +
        theme(legend.position = "bottom")
      
      ggplotly(p, tooltip = "text") %>%
        layout(hoverlabel = list(bgcolor = "white"))
    })
    
    # Top Pathways by NES - USING PATHWAY ID
    output$top_pathways_plot <- renderPlot({
      data <- filtered_data()
      if (nrow(data) == 0) return(NULL)
      
      # Get user-selected N value with validation
      top_n <- ifelse(is.null(input$top_n) || input$top_n < 1, 
                      5, 
                      min(input$top_n, 20))
      
      # Get top pathways by absolute NES magnitude
      top_data <- data %>%
        arrange(desc(abs(NES))) %>%
        head(top_n * 2)  # Get 2x top_n to ensure we cover both directions
      
      # If we have both directions, take top_n from each
      if (n_distinct(top_data$Direction) > 1) {
        top_activated <- top_data %>%
          filter(Direction == "Activated") %>%
          head(top_n)
        
        top_suppressed <- top_data %>%
          filter(Direction == "Suppressed") %>%
          head(top_n)
        
        top_data <- bind_rows(top_activated, top_suppressed)
      }
      
      # Order by NES (activated first, suppressed last)
      top_data <- top_data %>%
        arrange(desc(NES)) %>%
        mutate(
          Type = ifelse(NES > 0, "Activated", "Suppressed"),
          Type = factor(Type, levels = c("Activated", "Suppressed"))
        )
      
      if (nrow(top_data) == 0) return(NULL)
      
      # Calculate optimal wrap width
      max_char <- max(nchar(top_data$pathway), na.rm = TRUE)
      wrap_width <- case_when(
        max_char > 60 ~ 25,
        max_char > 40 ~ 30,
        TRUE ~ 40
      )
      
      # Create a single plot with consistent bar width
      p <- ggplot(top_data, aes(x = reorder(pathway, NES), y = NES, fill = Type)) +
        geom_bar(stat = "identity", width = 0.85) +
        scale_fill_manual(values = c("Activated" = "red", "Suppressed" = "blue")) +
        coord_flip() +
        labs(x = "", y = "Normalized Enrichment Score (NES)",
             title = paste("Top Pathways by NES Magnitude"),
             fill = "Direction") +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 12, lineheight = 0.8),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom"
        ) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = wrap_width))
      
      # Add a subtle separator between activated and suppressed
      if (n_distinct(top_data$Type) > 1) {
        # Find the index where activated meets suppressed
        sep_index <- sum(top_data$Type == "Activated")
      }
      
      return(p)
    })
    
    # Gene list for selected pathway
    output$gene_list <- renderUI({
      req(input$selected_pathway)
      
      # Split the selection into collection and pathway ID
      selected <- str_split(input$selected_pathway, "::", simplify = TRUE)
      collection <- selected[1]
      pathway_id <- selected[2]
      
      # Get the current filtered data
      data <- filtered_data()
      selected_data <- data[data$collection == collection & data$pathway == pathway_id, ]
      
      if (nrow(selected_data) == 0) {
        return(tags$p("Selected pathway not found in filtered data"))
      }
      
      # Construct key for leading_edge_genes
      gene_key <- paste(
        selected_data$file_key[1],
        collection,
        pathway_id,
        sep = "::"
      )
      
      debug_info$gene_key <- gene_key
      
      # Get genes from leading_edge_genes
      all_genes <- gsea_data$leading_edge_genes[[gene_key]]
      
      if (is.null(all_genes)) {
        return(tags$div(
          class = "alert alert-warning",
          HTML(paste(
            "Gene information not available for:<br>",
            "<b>Key:</b> ", gene_key
          ))
        ))
      }
      
      # Display in a scrollable box
      tags$div(
        style = "max-height: 300px; overflow-y: auto; border: 1px solid #ddd; padding: 10px;",
        tags$h5(paste("Leading Edge Genes for ", pathway_id, " (", length(all_genes), ")", sep = "")),
        tags$ul(
          lapply(all_genes, function(gene) tags$li(gene))
        )
      )
    })
    
    # GSEA Plot for selected pathway
    output$gsea_plot <- renderPlot({
      req(input$selected_pathway)
      
      # Split the selection into collection and pathway ID
      selected <- str_split(input$selected_pathway, "::", simplify = TRUE)
      collection <- selected[1]
      pathway_id <- selected[2]
      
      # Get the current filtered data
      data <- filtered_data()
      selected_data <- data[data$collection == collection & data$pathway == pathway_id, ]
      
      if (nrow(selected_data) == 0) {
        return(tags$div(
          class = "alert alert-danger",
          "Selected pathway not found in filtered data"
        ))
      }
      
      # Get the ranked genes for this analysis
      key <- selected_data$file_key[1]
      debug_info$ranked_key <- key
      ranked_genes <- gsea_data$ranked_genes_list[[key]]
      
      if (is.null(ranked_genes)) {
        return(tags$div(
          class = "alert alert-danger",
          HTML(paste(
            "Ranked genes not available for:<br>",
            "<b>Key:</b> ", key
          ))
        ))
      }
      
      # Construct key for leading_edge_genes
      gene_key <- paste(key, collection, pathway_id, sep = "::")
      pathway_genes <- gsea_data$leading_edge_genes[[gene_key]]
      
      if (is.null(pathway_genes)) {
        return(tags$div(
          class = "alert alert-danger",
          HTML(paste(
            "Leading edge genes not available for:<br>",
            "<b>Key:</b> ", gene_key
          ))
        ))
      }
      
      # Generate enrichment plot
      plotEnrichment(pathway = pathway_genes, stats = ranked_genes) +
        labs(title = pathway_id,
             subtitle = paste("NES:", round(selected_data$NES[1], 3), 
                              "| P.adj:", formatC(selected_data$padj[1], format = "e", digits = 2))) +
        theme(plot.title = element_text(size = 12, hjust = 0.5))
    })
    
    # Debug information output
    output$debug_info <- renderUI({
      req(gsea_data)
      
      tags$div(
        class = "well",
        tags$h4("Debug Information"),
        tags$p(tags$b("Filtered pathways:"), debug_info$filtered_count),
        tags$p(tags$b("Current gene key:"), debug_info$gene_key),
        tags$p(tags$b("Current ranked genes key:"), debug_info$ranked_key)
      )
    })
    
    # Observe table selection to update pathway selection
    observeEvent(input$results_table_rows_selected, {
      data <- filtered_data()
      if (!is.null(input$results_table_rows_selected)) {
        row_index <- input$results_table_rows_selected
        selected_row <- data[row_index, ]
        
        # Create composite key
        composite_key <- paste(selected_row$collection, selected_row$pathway, sep = "::")
        
        # Update the pathway selector
        updateSelectInput(session, "selected_pathway", selected = composite_key)
      }
    })
  })
}