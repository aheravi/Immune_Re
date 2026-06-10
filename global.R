suppressMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinyjs)
  library(shinyBS)
  library(shinycssloaders)
  library(shinyauthr)
  library(shinydisconnect)
  library(shinybusy)
  library(plotly)
  library(DT)
  library(tidyverse) # Loads ggplot2, magrittr, and tidyr automatically
  library(ggalluvial)
  library(vegan)
  library(janitor)
  library(pheatmap)
  library(units)
  library(sf)
})

theme_set(theme_bw(18))
  
  setwd("~/ShinyApps/ShinyReverb/")

  #meta_data = data.table::fread("/mnt/data/analysis/alireza/sample.csv")
  
# FIX: Map the 'qc' web prefix directly to the directory where the actual qualimap HTMLs sit
addResourcePath("qc", "/mnt/data/analysis/alireza/qualimap")

   
  pca_df  <- readRDS("pca_df_3d.rds")
  umap_df <- readRDS("umap_df_3d.rds")
  var_exp <- readRDS("var_exp.rds")


#meta_Choices <- c("Sequencer", "Year", "Protocol", "Name")
m <- list(l = 50, r = 50, b = 100, t = 100, pad = 4)

# Optimized CSS busy spinner
shiny_busy <- function() {
  HTML('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span data-display-if="$(\'html\').attr(\'class\')==\'shiny-busy\'"><i class="fa fa-spinner fa-pulse fa-fw" style="color:orange; font-size:100px; position:fixed; left:500px; top:300px; z-index:9999;"></i></span>')
}