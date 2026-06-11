dashboardPage(
  dashboardHeader(title = "Gene Expression Analysis "),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem(text = "About", tabName = "about", icon = icon("dashboard")),
      menuItem(text = "QC", icon = icon("chart-line"), startExpanded = FALSE,
               menuSubItem(text = "QC metrics", tabName = "qc_metrics", icon = icon("chart-line")),
               # Using UI container elements smoothly in sidebars instead of messy nested menuItems
               div(style = "padding: 10px 15px;",
                   selectInput("qc_file", label = "Select a Sample QC Report:", choices = NULL)
               ),
               menuSubItem(text = "HTML Summary", tabName = "qc", icon = icon("file-code"))
      ),
      
      menuItem(text = "Dimensional Reduction", tabName = "dim_red", icon = icon("cube")),
      menuItem(text = "Immune Repertoire", icon = icon("dna"), startExpanded = FALSE,
               menuSubItem(text = "Overview / QC", tabName = "ir_overview", icon = icon("table")),
               menuSubItem(text = "Diversity & Clonality", tabName = "ir_diversity", icon = icon("chart-pie")),
               menuSubItem(text = "Gene Usage", tabName = "ir_geneusage", icon = icon("chart-bar")),
               menuSubItem(text = "V-J Pairing", tabName = "ir_vj", icon = icon("th")),
               menuSubItem(text = "CDR3 Length", tabName = "ir_cdr3len", icon = icon("ruler-horizontal")),
               menuSubItem(text = "Clonotype Tracking", tabName = "ir_tracking", icon = icon("route"))
      ),
      menuItem(text = "Meta Data", tabName = "meta_data", icon = icon("table"))
    )
  ), 
  
  dashboardBody(
    useShinyjs(), 
    tags$head(
      tags$style(HTML("
        #Clonal_Sharing, #Clonal_Sharing_Venn {
          height: calc(100vh - 120px) !important;
        }
        #fullHeightImage {
          height: calc(100vh - 60px);
          width: auto;
          display: block;
          margin-left: auto;
          margin-right: auto;
        }
      "))
    ), 
    
    tabItems(
      tabItem(tabName = "about",
              fluidPage(
                h5(style = "line-height: 1.8;",
                   tags$br(), "A Shiny app for visualization of RNA-Seq analysis.",
                   tags$br(), tags$br(), tags$strong("Scope of Work:"),
                   tags$br(), "i. Quality control (QC)",
                   tags$br(), "ii. Differential gene expression analysis",
                   tags$br(), tags$br(), "Feedback welcome: bb@gmail.com"
                )
              )
      ),
      tabItem(tabName = "qc_metrics",
              withSpinner(DTOutput("QC_metrics"), type = 6, color = "#0dc5c1", color.background = "white")
      ),
      tabItem(tabName = "qc",
              withSpinner(uiOutput("qc_html"), type = 6, color = "#0dc5c1", color.background = "white")
      ),
      
      # FIX: Clean automated 2-column layout panel (No buttons)
      tabItem(tabName = "dim_red",
              fluidPage(
                fluidRow(
                  # Left Column: Compact footprint options
                  box(title = "Controls", width = 2, status = "warning", solidHeader = TRUE,
                      selectInput("plot_type", "Framework:", choices = c("PCA", "UMAP")),
                      selectInput("color_by", "Color Matrix:", choices = c("Group1", "Group2", "Group3", "Group4"))
                  ),
                  
                  # Right Column: Main high-real estate viewport panel
                  box(title = "3D Interactive Space", width = 10, status = "primary", solidHeader = TRUE,
                      withSpinner(
                        plotlyOutput("three_d_plot", height = "750px"),
                        type = 6, color = "#0dc5c1", color.background = "white"
                      )
                  )
                )
              )
      ),
      # --- Immune Repertoire: Overview / QC ---
      tabItem(tabName = "ir_overview",
              fluidPage(
                fluidRow(
                  box(title = "Chain Composition per Sample", width = 12, status = "primary", solidHeader = TRUE,
                      withSpinner(DTOutput("ir_chain_table"), type = 6, color = "#0dc5c1", color.background = "white")
                  )
                ),
                fluidRow(
                  box(title = "Sample Metadata", width = 12, status = "primary", solidHeader = TRUE,
                      withSpinner(DTOutput("ir_meta_table"), type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      
      # --- Immune Repertoire: Diversity & Clonality ---
      tabItem(tabName = "ir_diversity",
              fluidPage(
                fluidRow(
                  box(title = "Controls", width = 3, status = "warning", solidHeader = TRUE,
                      selectInput("ir_chain_filter_div", "Chain:",
                                  choices = c("All", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"),
                                  selected = "All"),
                      selectInput("ir_div_method", "Diversity Method:",
                                  choices = c("div" = "div", "chao1" = "chao1",
                                              "hill" = "hill", "gini.simp" = "gini.simp",
                                              "inv.simp" = "inv.simp", "gini" = "gini", "d50" = "d50"))
                  ),
                  box(title = "Diversity", width = 9, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_diversity_plot", height = "400px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                ),
                fluidRow(
                  box(title = "Clonal Proportions (Top Clones)", width = 12, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_clonality_plot", height = "400px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      
      # --- Immune Repertoire: Gene Usage ---
      tabItem(tabName = "ir_geneusage",
              fluidPage(
                fluidRow(
                  box(title = "Controls", width = 3, status = "warning", solidHeader = TRUE,
                      selectInput("ir_chain_filter_gene", "Chain:",
                                  choices = c("All", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"),
                                  selected = "IGK"),
                      selectInput("ir_gene_segment", "Gene Segment:",
                                  choices = c("V.name", "J.name", "D.name"), selected = "V.name"),
                      selectInput("ir_gene_sample", "Sample (single-sample plot):",
                                  choices = NULL),
                      numericInput("ir_topN", "Top N genes:", value = 20, min = 5, max = 100),
                      checkboxInput("ir_facet", "Facet by sample (all samples)", value = FALSE)
                  ),
                  box(title = "Gene Usage", width = 9, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_geneusage_plot", height = "500px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      
      # --- Immune Repertoire: V-J Pairing ---
      tabItem(tabName = "ir_vj",
              fluidPage(
                fluidRow(
                  box(title = "Controls", width = 3, status = "warning", solidHeader = TRUE,
                      selectInput("ir_vj_sample", "Sample:", choices = NULL),
                      selectInput("ir_chain_filter_vj", "Chain:",
                                  choices = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"),
                                  selected = "IGK")
                  ),
                  box(title = "V-J Pairing Heatmap", width = 9, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_vj_plot", height = "600px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      
      # --- Immune Repertoire: CDR3 Length ---
      tabItem(tabName = "ir_cdr3len",
              fluidPage(
                fluidRow(
                  box(title = "Controls", width = 3, status = "warning", solidHeader = TRUE,
                      selectInput("ir_chain_filter_cdr3", "Chain:",
                                  choices = c("All", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"),
                                  selected = "IGK")
                  ),
                  box(title = "CDR3 Length Distribution", width = 9, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_cdr3len_plot", height = "450px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      
      # --- Immune Repertoire: Clonotype Tracking ---
      tabItem(tabName = "ir_tracking",
              fluidPage(
                fluidRow(
                  box(title = "Controls", width = 3, status = "warning", solidHeader = TRUE,
                      selectInput("ir_chain_filter_track", "Chain:",
                                  choices = c("All", "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"),
                                  selected = "IGK"),
                      #verbatimTextOutput("debug_mode"),
                      radioButtons("ir_track_mode", "Track by:",
                                   choices = c("Top N from first sample" = "topn",
                                               "Specific CDR3 sequences" = "custom"),
                                   selected = "topn"),
                      conditionalPanel(
                        condition = "input.ir_track_mode == 'topn'",
                        numericInput("ir_track_topn", "Top N clonotypes:", value = 10, min = 1, max = 50)
                      ),
                      # add these to the box to see the results:
                      #YCGQATHLPPTF
                      #YCGQGTHFPPTF
                      conditionalPanel(
                        condition = "input.ir_track_mode == 'custom'",
                        textAreaInput("ir_track_custom", "CDR3.aa sequences (one per line):",
                                      value = "", rows = 5, width = "100%")
                      )
                      
                  ),
                  box(title = "Clonotype Tracking Across Samples", width = 9, status = "primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("ir_tracking_plot", height = "500px"),
                                  type = 6, color = "#0dc5c1", color.background = "white")
                  )
                )
              )
      ),
      tabItem(tabName = "meta_data",
              shiny_busy()
      )
    )
  ) 
)