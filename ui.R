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
      tabItem(tabName = "meta_data",
              shiny_busy()
      )
    )
  ) 
)