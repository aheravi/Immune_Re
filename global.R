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
  library(immunarch)
  library(dplyr)
  #library(ggplot2)
  library(stringr)
  #library(tidyr)
})

theme_set(theme_bw(18))
  
  setwd("~/ShinyApps/ShinyReverb/")

  #meta_data = data.table::fread("/mnt/data/analysis/alireza/sample.csv")
  
# FIX: Map the 'qc' web prefix directly to the directory where the actual qualimap HTMLs sit
addResourcePath("qc", "/mnt/data/analysis/alireza/qualimap")

   
  pca_df  <- readRDS("pca_df_3d.rds")
  umap_df <- readRDS("umap_df_3d.rds")
  var_exp <- readRDS("var_exp.rds")

  
  #immdata <- repLoad("/mnt/data/analysis/alireza/trust4/immunarch_input/", .mode = "AIRR")
  immdata <- readRDS("immdata.RDS")
  
  # Strip alleles
  ir_data <- lapply(immdata$data, function(df) {
    df$V.name <- sub("\\*.*", "", df$V.name)
    df$J.name <- sub("\\*.*", "", df$J.name)
    df$D.name <- sub("\\*.*", "", df$D.name)
    df$chain  <- substr(df$V.name, 1, 3)
    df$Proportion <- df$Clones / sum(df$Clones, na.rm = TRUE)
    df
  })
  
  # Build metadata
  ir_sample_names <- names(ir_data)
  ir_parsed <- str_match(ir_sample_names, "^(RV\\d+)_(\\d+)_(\\d+)D(.*)$")
  
  ir_meta <- data.frame(
    Sample  = ir_sample_names,
    Group1  = ir_parsed[,2],
    Group2  = ir_parsed[,3],
    Group3  = ir_parsed[,4],
    Timepoint_raw = ir_parsed[,5],
    stringsAsFactors = FALSE
  )
  ir_meta$Timepoint <- recode(ir_meta$Timepoint_raw,
                              "_1" = "01",
                              "_9" = "09",
                              .default = ir_meta$Timepoint_raw)
  
  ir_immdata <- list(data = ir_data, meta = ir_meta)
  
  # Helper: filter by chain
  ir_filter_chain <- function(data_list, chain) {
    if (chain == "All") return(data_list)
    lapply(data_list, function(df) df %>% filter(chain == !!chain))
  }
  
  
  
  
  

#meta_Choices <- c("Sequencer", "Year", "Protocol", "Name")
m <- list(l = 50, r = 50, b = 100, t = 100, pad = 4)

# Optimized CSS busy spinner
shiny_busy <- function() {
  HTML('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span data-display-if="$(\'html\').attr(\'class\')==\'shiny-busy\'"><i class="fa fa-spinner fa-pulse fa-fw" style="color:orange; font-size:100px; position:fixed; left:500px; top:300px; z-index:9999;"></i></span>')
}