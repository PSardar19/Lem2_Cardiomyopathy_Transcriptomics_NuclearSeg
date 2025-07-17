library(shiny)
library(DT)
library(plotly)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

temporalDGE_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    # Sidebar
    column(
      width = 3,  # You can adjust width as needed
      class = "sidebar-panel",
      style = "background-color: #1d2d44; color: white; min-height: 100vh; padding: 20px; border-right: 1px solid #34495E;",
      
      radioButtons(ns("tab_selector"), "Select Analysis:",
                   choices = c("Summary", "Temporal Trajectories", "Expression Heatmaps"),
                   selected = "Summary"),
      
      # Filters for Summary tab
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Summary'", ns("tab_selector")),
        selectInput(ns("pipeline"), "Select Pipeline:", choices = NULL),
        radioButtons(ns("condition"), "Select Condition:", choices = c("WT", "Mut")),
        selectInput(ns("transition"), "Select Transition:", choices = NULL),
        sliderInput(ns("padj"), "Adjusted p-value cutoff:", min = 0, max = 0.1, value = 0.05, step = 0.005),
        sliderInput(ns("logfc"), "Log2 Fold Change cutoff (absolute):", min = 0, max = 5, value = 1.2, step = 0.1)
      ),
      
      # Temporal Trajectories
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Temporal Trajectories'", ns("tab_selector")),
        selectInput(ns("pipeline_traj"), "Select Pipeline:", choices = NULL, selected = "NanoSeq"),
        selectInput(ns("gene_name"), "Select Gene:", choices = NULL, selected = "Myh6")
      ),
      
      # Heatmap tab filters
      conditionalPanel(
        condition = sprintf("input['%s'] == 'Expression Heatmaps'", ns("tab_selector")),
        selectInput(ns("pipeline_heat"), "Select Pipeline:", choices = NULL),
        selectInput(ns("transition_heat"), "Order genes by significance in:", choices = NULL),
        numericInput(ns("num_genes"), "Number of genes per heatmap:", value = 25, min = 1, max = 100, step = 1),
        sliderInput(ns("padj_heat"), "Adjusted p-value cutoff:", min = 0, max = 0.1, value = 0.05, step = 0.005),
        sliderInput(ns("logfc_heat"), "Log2 Fold Change cutoff (absolute):", min = 0, max = 5, value = 0.263, step = 0.1),
        radioButtons(ns("scale_option"), "Expression scaling:",
                     choices = c("Z-score per gene", "Raw normalized values"),
                     selected = "Z-score per gene"),
        radioButtons(ns("cluster_genes"), "Cluster genes:", choices = c("Yes", "No"), selected = "No")
      )
    ),
    
    # Main content
    column(
      width = 9,
      tabsetPanel(
        id = ns("main_tabs"),
        tabPanel("Summary",
                 tabsetPanel(
                   tabPanel("Plots",
                            fluidRow(
                              valueBoxOutput(ns("totalGenesBox")),
                              valueBoxOutput(ns("sigGenesBox")),
                              valueBoxOutput(ns("upGenesBox")),
                              valueBoxOutput(ns("downGenesBox")),
                              valueBoxOutput(ns("topUpBox")),
                              valueBoxOutput(ns("topDownBox"))
                            ),
                            br(),
                            fluidRow(
                              column(6, plotlyOutput(ns("volcano_plot"), height = "350px")),
                              column(6, plotlyOutput(ns("ma_plot"), height = "350px"))
                            )
                   ),
                   tabPanel("Filtered Table",
                            DTOutput(ns("temporal_dge_table"))
                   )
                 )
        ),
        tabPanel("Temporal Trajectories",
                 fluidRow(
                   column(6, plotlyOutput(ns("expression_trajectory_plot"), height = "500px")),
                   column(6, plotlyOutput(ns("temporal_trajectory_plot"), height = "500px"))
                 )
        ),
        tabPanel("Expression Heatmaps",
                 h3(textOutput(ns("heatmap_title"))),
                 fluidRow(
                   column(6,
                          h4("WT Condition - Top Upregulated Genes"),
                          plotlyOutput(ns("heatmap_wt_up"), height = "450px")
                   ),
                   column(6,
                          h4("WT Condition - Top Downregulated Genes"),
                          plotlyOutput(ns("heatmap_wt_down"), height = "450px")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(6,
                          h4("Mut Condition - Top Upregulated Genes"),
                          plotlyOutput(ns("heatmap_mut_up"), height = "450px")
                   ),
                   column(6,
                          h4("Mut Condition - Top Downregulated Genes"),
                          plotlyOutput(ns("heatmap_mut_down"), height = "450px")
                   )
                 ),
                 br(),
                 hr(),
                 fluidRow(
                   column(4,
                          radioButtons(ns("fc_corr_condition"), "Condition:",
                                       choices = c("WT", "Mut"),
                                       selected = "WT", inline = TRUE)
                   )
                 ),
                 plotlyOutput(ns("fc_correlation_plot"))
        )
      )
    )
  )
}

temporalDGE_server <- function(id, temporal_dge_data, expression_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Helper function to convert internal names to display names
    to_display_name <- function(original) {
      case_when(
        original == "P1_vs_E18.5." ~ "E18.5 → P1",
        original == "P3_vs_P1." ~ "P1 → P3",
        original == "P3_vs_E18.5." ~ "E18.5 → P3",
        TRUE ~ original
      )
    }
    
    # Helper function to convert display names back to internal names
    to_internal_name <- function(display) {
      case_when(
        display == "E18.5 → P1" ~ "P1_vs_E18.5.",
        display == "P1 → P3" ~ "P3_vs_P1.",
        display == "E18.5 → P3" ~ "P3_vs_E18.5.",
        TRUE ~ display
      )
    }
    
    observe({
      updateSelectInput(session, "pipeline", choices = names(temporal_dge_data))
      updateSelectInput(session, "pipeline_traj", choices = names(temporal_dge_data))
    })
    
    # For the summary tab
    observeEvent(input$pipeline, {
      req(input$pipeline)
      pipeline_list <- temporal_dge_data[[input$pipeline]]
      transitions <- unique(unlist(lapply(pipeline_list, names)))
      
      # Convert to display names
      display_names <- sapply(transitions, to_display_name)
      
      updateSelectInput(session, "transition", 
                        choices = setNames(transitions, display_names))
    })
    
    observe({
      selected_pipeline <- input$pipeline_traj %||% "NanoSeq"
      
      all_genes <- unique(unlist(lapply(temporal_dge_data[[selected_pipeline]], function(cond_list) {
        unlist(lapply(cond_list, function(df) df$gene_name))
      })))
      
      updateSelectInput(session, "pipeline_traj",
                        choices = names(temporal_dge_data),
                        selected = "NanoSeq")
      
      updateSelectInput(session, "gene_name",
                        choices = sort(all_genes),
                        selected = "Myh6")
    })
    
    observe({
      updateSelectInput(session, "pipeline_heat", choices = names(temporal_dge_data))
    })
    
    selected_data <- reactive({
      req(input$pipeline, input$condition, input$transition)
      
      # Use internal name for data access
      internal_transition <- input$transition
      
      temporal_dge_data[[input$pipeline]][[input$condition]][[internal_transition]]
    })
    
    filtered_summary_data <- reactive({
      req(input$pipeline, input$condition, input$transition)
      df <- temporal_dge_data[[input$pipeline]][[input$condition]][[input$transition]]  # Fixed variable name
      
      # Filter first
      df_filtered <- df %>% 
        filter(padj <= input$padj, abs(log2FoldChange) >= input$logfc)
      
      # Then format
      df_filtered <- df_filtered %>%
        mutate(across(where(is.numeric) & !matches("padj|pvalue"), ~ round(.x, 2))) %>%
        mutate(
          padj = formatC(padj, format = "e", digits = 2),
          pvalue = formatC(pvalue, format = "e", digits = 2)
        )
      
      if ("gene_name" %in% colnames(df_filtered)) {
        df_filtered <- df_filtered %>% select(gene_name, everything())
      }
      
      df_filtered      
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
        colors = c("Not Significant" = "grey", "Up" = "red", "Down" = "blue"),
        marker = list(size = 6),
        text = ~paste("Gene ID:", gene_id, "<br>Gene Name:", gene_name, "<br>log2FC:", round(log2FoldChange, 3)),
        hoverinfo = "text"
      ) %>%
        layout(
          xaxis = list(title = "Log2(Base Mean + 1)"),
          yaxis = list(title = "Log2 Fold Change"),
          title = "MA Plot"
        )
    })
    
    output$temporal_dge_table <- renderDT({
      df <- filtered_summary_data()
      datatable(df)
    })
    
    output$totalGenesBox <- renderValueBox({
      valueBox(length(unique(temporal_dge_data[[input$pipeline]][[input$condition]][[input$transition]]$gene_name)), "Total Genes")  # Fixed
    })
    
    output$sigGenesBox <- renderValueBox({
      df <- filtered_summary_data()
      valueBox(nrow(df), "Significant Genes")
    })
    
    output$upGenesBox <- renderValueBox({
      df <- filtered_summary_data()
      valueBox(sum(df$log2FoldChange > 0), "Upregulated Genes")
    })
    
    output$downGenesBox <- renderValueBox({
      df <- filtered_summary_data()
      valueBox(sum(df$log2FoldChange < 0), "Downregulated Genes")
    })
    
    output$topUpBox <- renderValueBox({
      df <- filtered_summary_data() %>% arrange(desc(log2FoldChange))
      valueBox(ifelse(nrow(df) > 0, df$gene_name[1], "NA"), "Top Upregulated")
    })
    
    output$topDownBox <- renderValueBox({
      df <- filtered_summary_data() %>% arrange(log2FoldChange)
      valueBox(ifelse(nrow(df) > 0, df$gene_name[1], "NA"), "Top Downregulated")
    })
    
    trajectory_data <- reactive({
      req(input$pipeline_traj, input$gene_name)
      pipeline_list <- temporal_dge_data[[input$pipeline_traj]]  # Fixed
      
      bind_rows(
        lapply(names(pipeline_list), function(cond) {
          cond_list <- pipeline_list[[cond]]
          bind_rows(lapply(names(cond_list), function(trans) {
            df <- cond_list[[trans]]
            gene_row <- df %>% filter(gene_name == input$gene_name)
            if (nrow(gene_row) == 0) return(NULL)
            data.frame(
              gene_name = input$gene_name,
              condition = cond,
              transition = trans,
              log2FoldChange = gene_row$log2FoldChange,
              padj = gene_row$padj
            )
          }))
        })
      )
    })
    
    output$temporal_trajectory_plot <- renderPlotly({
      df <- trajectory_data()
      req(nrow(df) > 0)
      
      transition_labels <- c(
        "P1_vs_E18.5." = "E18.5 → P1",
        "P3_vs_P1." = "P1 → P3",
        "P3_vs_E18.5." = "E18.5 → P3"
      )
      
      df <- df %>%
        mutate(
          pretty_transition = sapply(transition, to_display_name),
          pretty_transition = factor(pretty_transition, 
                                     levels = c("E18.5 → P1", "P1 → P3", "E18.5 → P3")),
          significance = case_when(
            padj < 0.001 ~ "***",
            padj < 0.01  ~ "**",
            padj < 0.05  ~ "*",
            TRUE         ~ ""
          )
        )
      
      plot_ly(
        data = df,
        x = ~pretty_transition,
        y = ~log2FoldChange,
        color = ~condition,
        colors = c("WT" = "#1f77b4", "Mut" = "salmon"),
        type = 'bar',
        text = ~significance,
        textposition = 'outside',
        hoverinfo = 'text',
        hovertext = ~paste("Transition:", pretty_transition,
                           "<br>Condition:", condition,
                           "<br>Log2FC:", round(log2FoldChange, 2),
                           "<br>Adj.P.Val:", formatC(padj, format = "e", digits = 2)),
        barmode = 'group'
      ) %>%
        layout(
          title = paste("Stepwise Transitions for", input$gene_name),
          xaxis = list(title = "Transition"),
          yaxis = list(title = "Log2 Fold Change"),
          margin = list(t = 60),
          legend = list(title = list(text = "Condition"))
        )
    })
    
    # Expression trajectory plot
    output$expression_trajectory_plot <- renderPlotly({
      req(input$pipeline_traj, input$gene_name)
      
      exp_df <- expression_data %>%
        filter(
          pipeline == input$pipeline_traj,
          external_gene_name == input$gene_name
        )
      
      validate(need(nrow(exp_df) > 0, "No expression data for selected gene/pipeline"))
      
      plot_ly(
        data = exp_df,
        x = ~stage,
        y = ~mean_expr,
        color = ~condition,
        colors = c("WT" = "#1f77b4", "Mut" = "salmon"),
        type = 'scatter',
        mode = 'lines+markers',
        linetype = ~condition,
        text = ~paste(
          "Stage:", stage,
          "<br>Condition:", condition,
          "<br>Mean:", round(mean_expr, 2)),
        hoverinfo = 'text',
        error_y = list(array = ~sd_expr, 
                       color = '#000000',
                       thickness = 1)
      ) %>%
        layout(
          title = paste(input$pipeline_traj, "Expression for", input$gene_name),
          xaxis = list(title = "Developmental Stage", 
                       categoryorder = "array",
                       categoryarray = c("E18.5", "P1", "P3")),
          yaxis = list(title = "Normalized Expression"),
          margin = list(t = 60),
          showlegend = TRUE
        )
    })
    
    output$fc_correlation_plot <- renderPlotly({
      req(input$pipeline_traj, input$fc_corr_condition)
      
      pipeline_list <- temporal_dge_data[[input$pipeline_traj]]
      condition <- input$fc_corr_condition
      
      transitions_to_compare <- c("P1_vs_E18.5.", "P3_vs_P1.")
      
      all_dfs <- lapply(transitions_to_compare, function(trans) {
        df <- pipeline_list[[condition]][[trans]]
        df %>%
          select(gene_name, log2FoldChange, padj) %>%
          rename(
            !!paste0("log2FC_", trans) := log2FoldChange,
            !!paste0("padj_", trans) := padj
          )
      })
      
      merged_df <- reduce(all_dfs, full_join, by = "gene_name") %>%
        drop_na()
      
      # Assign significance and set proper factor levels
      merged_df <- merged_df %>%
        mutate(
          significance = case_when(
            padj_P1_vs_E18.5. < 0.05 & padj_P3_vs_P1. < 0.05 ~ "Both",
            padj_P1_vs_E18.5. < 0.05 ~ "P1vsE18.5",
            padj_P3_vs_P1. < 0.05 ~ "P3vsP1",
            TRUE ~ "Not significant"
          ),
          significance = factor(significance, levels = c("Not significant", "P1vsE18.5", "P3vsP1", "Both"))
        )
      
      # Compute correlation
      r <- cor(merged_df$log2FC_P1_vs_E18.5., merged_df$log2FC_P3_vs_P1.)
      r_squared <- r^2
      
      # Create the plot
      plot_ly(
        data = merged_df,
        x = ~log2FC_P1_vs_E18.5.,
        y = ~log2FC_P3_vs_P1.,
        type = 'scatter',
        mode = 'markers',
        text = ~gene_name,
        color = ~significance,
        colors = c("grey", "blue", "red", "purple"),
        marker = list(size = 6, opacity = 0.7)
      ) %>%
        layout(
          title = paste("Correlation of Fold Changes Across Transitions -", condition),
          xaxis = list(title = "log2FC: P1 vs E18.5"),
          yaxis = list(title = "log2FC: P3 vs P1"),
          showlegend = TRUE,
          annotations = list(
            # Quadrant labels fixed to plot area corners (adjust x/y to suit)
            list(x = 0.1, y = 0.9, xref = "paper", yref = "paper",
                 text = "↓ P1vsE18.5<br>↑ P3vsP1", showarrow = FALSE),
            list(x = 0.9, y = 0.9, xref = "paper", yref = "paper",
                 text = "↑ in both", showarrow = FALSE),
            list(x = 0.1, y = 0.1, xref = "paper", yref = "paper",
                 text = "↓ in both", showarrow = FALSE),
            list(x = 0.9, y = 0.1, xref = "paper", yref = "paper",
                 text = "↑ P1vsE18.5<br>↓ P3vsP1", showarrow = FALSE),
            # R & R² below title (centered horizontally)
            list(x = 0.5, y = 1.05, xref = "paper", yref = "paper",
                 text = paste0("R = ", round(r, 2), " &nbsp;&nbsp; R² = ", round(r_squared, 2)),
                 showarrow = FALSE, font = list(size = 14), align = 'center')
          )
        )
      
    })
    # Update transition choices for heatmap tab
    observeEvent(input$pipeline_heat, {
      req(input$pipeline_heat)
      pipeline_list <- temporal_dge_data[[input$pipeline_heat]]
      transitions <- unique(unlist(lapply(pipeline_list, names)))
      
      # Convert to display names
      display_names <- sapply(transitions, to_display_name)
      
      updateSelectInput(session, "transition_heat", 
                        choices = setNames(transitions, display_names))
    })
    
    # Dynamic title for heatmaps
    output$heatmap_title <- renderText({
      req(input$transition_heat)
      paste("Expression Patterns for Genes Significant in:", 
            to_display_name(input$transition_heat))
    })
    
    # Get top genes ordered by significance in selected transition
    top_genes <- reactive({
      req(input$pipeline_heat, input$transition_heat, input$num_genes,
          input$padj_heat, input$logfc_heat)
      
      # Use internal name for data access
      internal_transition <- input$transition_heat
      
      pipeline_data <- temporal_dge_data[[input$pipeline_heat]]
      
      # Function to get top genes for a condition
      get_top_genes <- function(condition) {
        trans_df <- pipeline_data[[condition]][[input$transition_heat]]
        if (is.null(trans_df)) return(list(up = character(0), down = character(0)))
        
        # Filter significant genes
        sig_genes <- trans_df %>% 
          filter(padj <= input$padj_heat, 
                 abs(log2FoldChange) >= input$logfc_heat)
        
        # Get top upregulated
        up_genes <- sig_genes %>%
          filter(log2FoldChange > 0) %>%
          arrange(padj, desc(log2FoldChange)) %>%
          head(input$num_genes) %>%
          pull(gene_name)
        
        # Get top downregulated
        down_genes <- sig_genes %>%
          filter(log2FoldChange < 0) %>%
          arrange(padj, log2FoldChange) %>%  # Most negative first
          head(input$num_genes) %>%
          pull(gene_name)
        
        list(up = up_genes, down = down_genes)
      }
      
      # Get genes for both conditions
      wt_genes <- get_top_genes("WT")
      mut_genes <- get_top_genes("Mut")
      
      list(
        wt_up = wt_genes$up,
        wt_down = wt_genes$down,
        mut_up = mut_genes$up,
        mut_down = mut_genes$down
      )
    })
    
    # Create heatmap function - UPDATED FOR FLAT DF STRUCTURE
    create_heatmap <- function(genes, cond) {
      if (length(genes) == 0) {
        return(plot_ly() %>% 
                 layout(title = paste("No significant genes found for", cond)))
      }
      
      # Get expression data - simplified for flat structure
      exp_df <- expression_data %>%
        filter(pipeline == input$pipeline_heat,
               external_gene_name %in% genes,
               condition == cond) %>%
        mutate(stage = factor(stage, levels = c("E18.5", "P1", "P3")))
      
      # Apply scaling
      if (input$scale_option == "Z-score per gene") {
        exp_df <- exp_df %>%
          group_by(external_gene_name) %>%
          mutate(expr_scaled = scale(mean_expr)) %>%
          ungroup()
      } else {
        exp_df$expr_scaled <- exp_df$mean_expr
      }
      
      # Create matrix - no need for pivoting since we have direct mapping
      heat_mat <- exp_df %>%
        select(external_gene_name, stage, expr_scaled) %>%
        pivot_wider(names_from = stage, values_from = expr_scaled) %>%
        column_to_rownames("external_gene_name") %>%
        as.matrix()
      
      # Preserve significance order (unless clustering)
      if (input$cluster_genes == "No") {
        # Order by the original gene significance ranking
        heat_mat <- heat_mat[genes[genes %in% rownames(heat_mat)], ]
      } else if (nrow(heat_mat) > 1) {
        # Cluster genes if requested
        hc <- hclust(dist(heat_mat))
        heat_mat <- heat_mat[hc$order, ]
      }
      
      # Create heatmap
      plot_ly(
        x = colnames(heat_mat),
        y = rownames(heat_mat),
        z = heat_mat,
        type = "heatmap",
        colors = viridis::viridis(100),
        colorbar = list(title = ifelse(input$scale_option == "Z-score per gene", 
                                       "Z-score", "Expression")),
        hoverinfo = "x+y+z",
        hovertext = ~paste("Gene: ", rownames(heat_mat), "<br>",
                           "Stage: ", colnames(heat_mat), "<br>",
                           "Value: ", round(heat_mat, 2))
      ) %>%
        layout(
          xaxis = list(title = "Developmental Stage"),
          yaxis = list(title = "Genes", tickfont = list(size = 9)),
          margin = list(l = 150)
        )
    }
    
    # Render heatmaps - UNCHANGED
    output$heatmap_wt_up <- renderPlotly({
      create_heatmap(top_genes()$wt_up, "WT")
    })
    
    output$heatmap_wt_down <- renderPlotly({
      create_heatmap(top_genes()$wt_down, "WT")
    })
    
    output$heatmap_mut_up <- renderPlotly({
      create_heatmap(top_genes()$mut_up, "Mut")
    })
    
    output$heatmap_mut_down <- renderPlotly({
      create_heatmap(top_genes()$mut_down, "Mut")
    })
    
  })
}