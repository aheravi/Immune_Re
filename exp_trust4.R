

  # correlate expression with trust4
  
  library(dplyr)
  library(immunarch)
  library(stringr)
  
  ir_meta = readRDS("~/ShinyApps/ShinyReverb/ir_meta")
  ir_immdata = readRDS("~/ShinyApps/ShinyReverb/ir_immdata")
  ir_filter_chain = readRDS("~/ShinyApps/ShinyReverb/ir_filter_chain.RDS")
  

  # Building the repertoire summary table
  repertoire_summary <- bind_rows(lapply(names(ir_immdata$data), function(s) {
  df <- ir_immdata$data[[s]]
  
  data.frame(
    Sample = s,
    
    # Total clones per chain
    IGH_clones = sum(df$Clones[df$chain == "IGH"], na.rm = TRUE),
    IGK_clones = sum(df$Clones[df$chain == "IGK"], na.rm = TRUE),
    IGL_clones = sum(df$Clones[df$chain == "IGL"], na.rm = TRUE),
    TRA_clones = sum(df$Clones[df$chain == "TRA"], na.rm = TRUE),
    TRB_clones = sum(df$Clones[df$chain == "TRB"], na.rm = TRUE),
    
    # Total unique clonotypes (repertoire size) per chain
    IGH_unique = sum(df$chain == "IGH"),
    IGK_unique = sum(df$chain == "IGK"),
    TRB_unique = sum(df$chain == "TRB"),
    
    stringsAsFactors = FALSE
  )
}))

  # Add diversity metrics
  div_chao <- repDiversity(ir_immdata$data, .method = "chao1")
  div_shannon <- repDiversity(ir_immdata$data, .method = "div")
  
  # These are matrices, not data frames - convert first
  div_chao_df <- as.data.frame(div_chao)
  div_chao_df$Sample <- rownames(div_chao_df)
  
  div_shannon_df <- as.data.frame(div_shannon)
  div_shannon_df$Sample <- rownames(div_shannon_df)
  
  # Check column names
  head(div_chao_df)
  head(div_shannon_df)
  
  repertoire_summary <- repertoire_summary %>%
    left_join(div_chao_df %>% select(Sample, Chao1 = Estimator), by = "Sample") %>%
    left_join(div_shannon_df %>% select(Sample, Shannon = Value), by = "Sample")
  
  
  # Add top-clone clonality (e.g. top 10 proportion)
  clon <- repClonality(ir_immdata$data, .method = "top", .head = c(10))
  clon_df <- as.data.frame(clon)
  clon_df$Sample <- rownames(clon_df)
  
  head(clon_df)
  repertoire_summary <- repertoire_summary %>%
    left_join(clon_df, by = "Sample")
  
  head(repertoire_summary)
  
  # fixing smaple names like D_0 to D00
  if (FALSE){ # fixing smaple names like D_0 to D00
  # fixing smaple names like D_0 to D00
  fix_sample_name <- function(x) {
    m <- stringr::str_match(x, "^(.*D)_(\\d+)$")
    ifelse(!is.na(m[,1]),
           paste0(m[,2], sprintf("%02d", as.integer(m[,3]))),
           x)
  }
  repertoire_summary$Sample <- fix_sample_name(repertoire_summary$Sample)
  ir_immdata$meta$Sample    <- fix_sample_name(ir_immdata$meta$Sample)
  
  head(repertoire_summary)
  table(repertoire_summary$Sample)[table(repertoire_summary$Sample) > 1]
  }
  
  
  
  
  
  # preparing expression data
  # tpm_matrix: rows = genes (gene symbols), columns = sample names
  #tpm_matrix <- read.delim("/mnt/data/analysis/alireza/StringTie/final_tpm_matrix.txt", row.names = 1, check.names = FALSE)
  tpm_matrix = readRDS("/mnt/data/analysis/alireza/StringTie/final_tpm_matrix.RDS")
  
  library(stringr)
  
  orig_names <- colnames(tpm_matrix)  # original, with D_1/D_9 style names still present
  
  # For each duplicated new name, find the original column names that map to it
  for (nm in unique(dup_names)) {
    # Reverse-match: which original names become this new name?
    matches <- orig_names[fix_sample_name(orig_names) == nm & orig_names != "Gene_ID"]
    
    if (length(matches) == 2) {
      x <- tpm_matrix[[matches[1]]]
      y <- tpm_matrix[[matches[2]]]
      r <- cor(x, y, method = "spearman")
      cat(sprintf("%-20s : %-20s vs %-20s | r = %.4f\n", nm, matches[1], matches[2], r))
    } else {
      cat(sprintf("%-20s : found %d matches: %s\n", nm, length(matches), paste(matches, collapse=", ")))
    }
  }
  
  
  colnames(tpm_matrix) <- ifelse(colnames(tpm_matrix) == "Gene_ID",
                                 "Gene_ID",
                                 fix_sample_name(colnames(tpm_matrix)))
  colnames(tpm_matrix)
  new_names <- ifelse(colnames(tpm_matrix) == "Gene_ID",
                      "Gene_ID",
                      fix_sample_name(colnames(tpm_matrix)))
  
  dup_names <- new_names[duplicated(new_names)]
  dup_names
  
  # Find which positions have the duplicated name
  which(colnames(tpm_matrix) == "RV09_05_21084D01")
  
  
  
  # Pick immune marker genes of interest
  genes_of_interest <- c("IGHM", "IGHG1", "IGKC", "IGLC1", "MZB1", "CD79A", "CD19",
                         "CD3D", "CD3E", "TRBC1", "TRAC", "CD8A", "CD4", "MKI67")
  
  # Subset and transpose to samples × genes
  expr_subset <- tpm_matrix[rownames(tpm_matrix) %in% genes_of_interest, ]
  expr_subset_t <- as.data.frame(t(expr_subset))
  expr_subset_t$Sample <- rownames(expr_subset_t)
  
  # Example if expression matrix uses different naming
  sample_map <- data.frame(
    trust4_name = c("RV09_05_21084D00", "RV09_05_21084D01", "RV09_05_21084D09"),
    expr_name   = c("D00_sample", "D01_sample", "D09_sample")  # adjust to your actual names
  )
  
  expr_subset_t <- expr_subset_t %>%
    left_join(sample_map, by = c("Sample" = "expr_name")) %>%
    mutate(Sample = trust4_name) %>%
    select(-trust4_name)
  
  # Merge and correlate
  merged <- repertoire_summary %>%
    left_join(expr_subset_t, by = "Sample")
  
  head(merged)
  
  #Pairwise correlation (e.g. IGK clones vs IGKC expression):
  cor.test(merged$IGK_clones, merged$IGKC, method = "spearman")
  
  # scatter plot:
  library(ggplot2)
  
  ggplot(merged, aes(x = IGKC, y = IGK_clones, label = Sample)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    geom_text(vjust = -1, size = 3) +
    labs(x = "IGKC expression (TPM)", y = "IGK clone count",
         title = "IGKC expression vs IGK repertoire size")
  
  #Correlation matrix across multiple metrics:
  library(corrplot)
  
  cor_data <- merged %>%
    select(IGH_clones, IGK_clones, TRB_clones, Chao1, Shannon,
           IGHM, IGKC, CD3D, CD79A, MZB1) %>%
    na.omit()
  
  cor_mat <- cor(cor_data, method = "spearman")
  
  corrplot(cor_mat, method = "circle", type = "upper",
           tl.cex = 0.7, tl.col = "black")
  
  #With only 3–5 samples, statistical correlation is limited
  #With current sample size (3 timepoints × few subjects), formal correlation tests (cor.test) won't have meaningful p-values. This is better treated as:
  #Exploratory/descriptive — does the trend make biological sense (e.g., does B-cell clonal expansion at D09 coincide with increased IGHG1/MZB1 expression)?
  #Per-subject longitudinal trajectories — plot expression and repertoire metrics over timepoints for the same subject, side by side, rather than computing a single correlation coefficient:

  merged_long <- merged %>%
    left_join(ir_immdata$meta, by = "Sample") %>%
    select(Sample, Timepoint, Group1, IGK_clones, IGKC, Shannon) %>%
    tidyr::pivot_longer(cols = c(IGK_clones, IGKC, Shannon),
                        names_to = "Metric", values_to = "Value")
  
  ggplot(merged_long, aes(x = Timepoint, y = Value, color = Metric, group = Metric)) +
    geom_line() + geom_point() +
    facet_wrap(~Group1, scales = "free_y") +
    labs(title = "Repertoire metrics vs gene expression over time")
  
  
  