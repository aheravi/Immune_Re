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
  
}