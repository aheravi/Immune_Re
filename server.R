# Global resource path for other plot assets
addResourcePath(
  prefix = 'my_plots', 
  directoryPath = '/mnt/data/analysis/alireza/temp'
)

function(input, output, session) {
  
  values <- reactiveValues(authenticated = FALSE)
  
  dataModal <- function() {
    modalDialog(
      textInput("username2", "Username:"),
      passwordInput("password2", "Password:"),
      footer = tagList(actionButton("ok", "OK")),
      easyClose = FALSE,
      fade = FALSE,
      tags$script(HTML("
        $(document).on('shown.bs.modal', function() {
          $('#username2').focus();
        });
        $(document).off('keypress').on('keypress', function(e) {
          if (e.which == 13 && $('#shiny-modal').is(':visible')) {
            $('#ok').click();
          }
        });
      "))
    )
  }
  
  # Show login modal on launch
  observe({
    if (!values$authenticated) {
      showModal(dataModal())
    }
  })
  
  # Authentication logic

  # Authentication logic
  observeEvent(input$ok, {
    Username <- input$username2
    Password <- input$password2
    
    if (nchar(Username) > 0 && nchar(Password) > 0) {
      show_modal_spinner(text = "Authenticating ... Please wait.")
      
      if (Username == "Re" && Password == "AM") {
        values$authenticated <- TRUE
        removeModal()
        showNotification("Login successful!", type = "message")
        remove_modal_spinner()
        
        # --- SUCCESSFUL LOGIN: NOW SCAN AND POPULATE DROPDOWN ---
        base_qc_dir <- "/mnt/data/analysis/alireza/qualimap/"
        if (dir.exists(base_qc_dir)) {
          qc_folders <- list.dirs(base_qc_dir, full.names = FALSE, recursive = FALSE)
          
          # Only update if folders actually exist, otherwise show a placeholder
          if (length(qc_folders) > 0) {
            updateSelectInput(session, "qc_file", 
                              choices = qc_folders, 
                              selected = qc_folders[1])
          } else {
            updateSelectInput(session, "qc_file", 
                              choices = "No subfolders found in path")
          }
        } else {
          updateSelectInput(session, "qc_file", 
                            choices = "Directory does not exist")
        }
        # --------------------------------------------------------
        
      } else {
        remove_modal_spinner()
        showNotification("Invalid Username or Password.", type = "error")
        showModal(dataModal()) 
      }
    } else {
      showNotification("Please enter both Username and Password.", type = "warning")
      showModal(dataModal())  
    }
  })  

  
  # QC Datatable Output
  output$QC_metrics <- renderDT({
    req(exists("qc_data")) # Prevents rendering errors if qc_data hasn't loaded yet
    
    # Dynamic page length fallback calculation
    p_len <- if(exists("meta_data")) nrow(meta_data) else 10
    
    datatable(
      qc_data,
      filter = "top",
      rownames = FALSE,
      extensions = 'Buttons',
      options = list(
        scrollX = TRUE,          
        scrollY = "50vh",          
        autoWidth = TRUE,
        pageLength = p_len,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      )
    )
  })
  
  # Dynamic Qualimap HTML Report Viewer
  observeEvent(input$qc_file, {
    req(input$qc_file)
    
    base_qc_dir <- "/mnt/data/analysis/alireza/qualimap/"
    
    # 1. List all subdirectories inside the main qualimap folder
    all_dirs <- list.dirs(base_qc_dir, full.names = FALSE, recursive = FALSE)
    
    # 2. Match the folder name using your input selection pattern
    # e.g., if input$qc_file is "RV09_05_2100", match that folder
    pattern <- paste0("^", input$qc_file)
    matched_dir <- all_dirs[grepl(pattern, all_dirs)]
    
    # Safely handle missing directory matches
    if (length(matched_dir) == 0) {
      output$qc_html <- renderUI({
        h3(style = "color: red; padding: 20px;", "No matching QC results folder found.")
      })
      return(NULL)
    }
    
    # 3. Target the specific html report within that matched subfolder
    # Constructed path will be like: "qc/RV09_05_2100/qualimapReport.html"
    qc_url <- paste0("qc/", matched_dir[1], "/qualimapReport.html")
    
    # Verify the actual file exists on disk before attempting to serve it
    full_disk_path <- file.path(base_qc_dir, matched_dir[1], "qualimapReport.html")
    if (!file.exists(full_disk_path)) {
      output$qc_html <- renderUI({
        h3(style = "color: red; padding: 20px;", "Folder found, but qualimapReport.html is missing inside it.")
      })
      return(NULL)
    }
    
    # Render into iframe smoothly
    output$qc_html <- renderUI({
      tags$iframe(
        src = qc_url,
        style = "width:100%; height: calc(100vh - 160px); border:none; background: white;"
      )
    })
  })
  
  output$meta_data1 <- renderDT({
    
    datatable(meta_data, filter = "top", rownames = FALSE, extensions = 'Buttons',
              options = list (
                scrollX = TRUE,
                scrollY = TRUE,
                autoWidth = TRUE,
                pageLength = nrow(meta_data),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              ))
    
  })
  
  # Inside your server.R function(input, output, session) block:
  
  output$three_d_plot <- renderPlotly({
    req(input$plot_type, input$color_by)
    
    if (input$plot_type == "PCA") {
      
      plot_ly(pca_df, 
              x = ~PC1, y = ~PC2, z = ~PC3, 
              color = pca_df[[input$color_by]], 
              colors = "Set1", # <--- FORCES HIGH CONTRAST CATEGORICAL COLORS
              text = ~Sample, hoverinfo = "text",
              type = "scatter3d", mode = "markers", size = I(50)) %>%
        layout(scene = list(
          xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
          yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
          zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
        ))
      
    } else {
      
      plot_ly(umap_df, 
              x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, 
              color = umap_df[[input$color_by]], 
              colors = "Set1", # <--- FORCES HIGH CONTRAST CATEGORICAL COLORS
              text = ~Sample, hoverinfo = "text",
              type = "scatter3d", mode = "markers", size = I(50))
      
    }
  })

  
  
  
  # ---- Populate sample dropdowns ----
  observe({
    updateSelectInput(session, "ir_gene_sample",
                      choices = names(ir_immdata$data),
                      selected = names(ir_immdata$data)[1])
    updateSelectInput(session, "ir_vj_sample",
                      choices = names(ir_immdata$data),
                      selected = names(ir_immdata$data)[1])
  })
  
  # ---- Overview / QC ----
  output$ir_chain_table <- renderDT({
    chain_counts <- bind_rows(lapply(names(ir_immdata$data), function(s) {
      as.data.frame(table(ir_immdata$data[[s]]$chain)) %>%
        rename(Chain = Var1, Count = Freq) %>%
        mutate(Sample = s)
    })) %>%
      pivot_wider(names_from = Chain, values_from = Count, values_fill = 0)
    
    datatable(chain_counts, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$ir_meta_table <- renderDT({
    datatable(ir_immdata$meta, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  library(plotly)
  
  # ---- Diversity & Clonality ----
  output$ir_diversity_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_div)
    div <- repDiversity(dat, .method = input$ir_div_method)
    p <- vis(div)
    ggplotly(p)
  })
  
  output$ir_clonality_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_div)
    clon <- repClonality(dat, .method = "top", .head = c(10, 100, 1000))
    p <- vis(clon)
    ggplotly(p)
  })
  
  # ---- Gene Usage ----
  output$ir_geneusage_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_gene)
    seg <- input$ir_gene_segment
    
    if (input$ir_facet) {
      gene_usage_all <- bind_rows(lapply(names(dat), function(s) {
        dat[[s]] %>%
          count(.data[[seg]], name = "n") %>%
          mutate(Sample = s)
      }))
      names(gene_usage_all)[1] <- "Gene"
      
      p <- ggplot(gene_usage_all, aes(x = reorder(Gene, -n), y = n,
                                      text = paste0("Gene: ", Gene, "<br>Count: ", n,
                                                    "<br>Sample: ", Sample))) +
        geom_bar(stat = "identity") +
        facet_wrap(~Sample, scales = "free_y") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
        labs(x = seg, y = "Count")
    } else {
      df <- dat[[input$ir_gene_sample]]
      gene_usage <- df %>% count(.data[[seg]], sort = TRUE)
      names(gene_usage)[1] <- "Gene"
      top_genes <- gene_usage %>% slice_max(n, n = input$ir_topN)
      
      p <- ggplot(top_genes, aes(x = reorder(Gene, -n), y = n,
                                 text = paste0("Gene: ", Gene, "<br>Count: ", n))) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = seg, y = "Count",
             title = paste0("Top ", input$ir_topN, " ", seg, " (", input$ir_gene_sample, ")"))
    }
    
    ggplotly(p, tooltip = "text")
  })
  
  # ---- V-J Pairing ----
  output$ir_vj_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_vj)
    df <- dat[[input$ir_vj_sample]]
    req(nrow(df) > 0)
    
    vj_table <- df %>%
      count(V.name, J.name) %>%
      complete(V.name, J.name, fill = list(n = 0))
    
    p <- ggplot(vj_table, aes(x = J.name, y = V.name, fill = n,
                              text = paste0("V: ", V.name, "<br>J: ", J.name, "<br>Count: ", n))) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkred") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6)) +
      labs(title = paste("V-J pairing -", input$ir_vj_sample,
                         "(", input$ir_chain_filter_vj, ")"))
    
    ggplotly(p, tooltip = "text")
  })
  
  # ---- CDR3 Length ----
  output$ir_cdr3len_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_cdr3)
    
    cdr3_len <- bind_rows(lapply(names(dat), function(s) {
      dat[[s]] %>% mutate(len = nchar(CDR3.aa), Sample = s)
    }))
    
    p <- ggplot(cdr3_len, aes(x = len, fill = Sample)) +
      geom_density(alpha = 0.4) +
      labs(x = "CDR3 length (aa)", y = "Density",
           title = paste("CDR3 length distribution -", input$ir_chain_filter_cdr3))
    
    ggplotly(p)
  })
  
  # ---- Clonotype Tracking ----
  output$ir_tracking_plot <- renderPlotly({
    dat <- ir_filter_chain(ir_immdata$data, input$ir_chain_filter_track)
    
    if (input$ir_track_mode == "topn") {
      tc <- trackClonotypes(dat, list(1, input$ir_track_topn), .col = "aa")
    } else {
      seqs <- strsplit(input$ir_track_custom, "\n")[[1]]
      seqs <- trimws(seqs)
      seqs <- seqs[seqs != ""]
      req(length(seqs) > 0)
      
      # Check which sequences actually exist in the (chain-filtered) data
      all_cdr3 <- unique(unlist(lapply(dat, function(df) df$CDR3.aa)))
      found <- seqs[seqs %in% all_cdr3]
      missing <- seqs[!seqs %in% all_cdr3]
      
      validate(
        need(length(found) > 0,
             paste0("None of the entered CDR3 sequences were found in the '",
                    input$ir_chain_filter_track, "' chain data.\n",
                    "Missing: ", paste(missing, collapse = ", ")))
      )
      
      tc <- trackClonotypes(dat, found, .col = "aa")
    }
    
    p <- vis(tc)
    ggplotly(p)
  })
  
  
  ## debugging CDR3 inputs:
  #output$debug_mode <- renderText({ paste("Mode:", input$ir_track_mode) })
  }