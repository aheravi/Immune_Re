
# stringTie :
# stringTie :

  
  library(data.table)
  library(matrixStats)
  library(umap)
  library(tidyr)
  #library(ggplot2)
  library(plotly)
  

  setwd("/mnt/data/analysis/alireza/StringTie/")
  
  tpm_matrix = readRDS("final_tpm_matrix.RDS")
  
  #3D plots:
  
  # --- 1. Prepare Matrix (From your previous step) ---
  mat <- as.matrix(tpm_matrix[, -1, with = FALSE])
  rownames(mat) <- tpm_matrix$Gene_ID
  log_mat <- log2(mat + 1)
  
  # Get top 1000 variable genes
  gene_vars <- rowVars(log_mat)
  variable_mat <- log_mat[order(gene_vars, decreasing = TRUE)[1:2000], ]
  transposed_mat <- t(variable_mat)
  
  # --- 2. Generate 3D PCA Data ---
  pca_3d <- prcomp(transposed_mat, center = TRUE, scale. = FALSE)
  pca_df_3d <- as.data.frame(pca_3d$x[, 1:3]) # Keep PC1, PC2, and PC3
  pca_df_3d$Sample <- rownames(pca_df_3d)
  
  # Extract your metadata groupings
  #pca_df_3d <- separate(pca_df_3d, Sample, into = c("Group1", "Group2"), 
  #                      sep = "_", remove = FALSE, extra = "drop")
  
  pca_df_3d <- tidyr::extract(
    pca_df_3d, 
    col = Sample, 
    into = c("Group1", "Group2", "Group3", "Group4"), 
    # Slices at the first two underscores, then separates numbers from "D..."
    regex = "^([^_]+)_([^_]+)_([0-9]+)(D.*)$", 
    remove = FALSE
  )
  
  # Calculate % variance explained for axis labels
  var_exp <- round(100 * (pca_3d$sdev^2 / sum(pca_3d$sdev^2)), 1)
  
  
  # --- 3. Generate 3D UMAP Data ---
  set.seed(42)
  umap_config_3d <- umap.defaults
  umap_config_3d$n_components <- 3 # Force UMAP to calculate 3 dimensions
  
  umap_3d <- umap(transposed_mat, config = umap_config_3d)
  umap_df_3d <- as.data.frame(umap_3d$layout)
  colnames(umap_df_3d) <- c("UMAP1", "UMAP2", "UMAP3")
  umap_df_3d$Sample <- rownames(umap_df_3d)
  
  # Extract groupings
  #umap_df_3d <- separate(umap_df_3d, Sample, into = c("Group1", "Group2"), 
  #                       sep = "_", remove = FALSE, extra = "drop")
  
  umap_df_3d <- tidyr::extract(
    umap_df_3d, 
    col = Sample, 
    into = c("Group1", "Group2", "Group3", "Group4"), 
    # Slices at the first two underscores, then separates numbers from "D..."
    regex = "^([^_]+)_([^_]+)_([0-9]+)(D.*)$", 
    remove = FALSE
  )
  
  
  # Create the 3D PCA Plot
  pca_plot_3d <- plot_ly(pca_df_3d, 
                         x = ~PC1, y = ~PC2, z = ~PC3, 
                         color = ~Group1, symbol = ~Group2,
                         text = ~Sample, hoverinfo = "text",
                         type = "scatter3d", mode = "markers", size = I(50)) %>%
    layout(scene = list(
      xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
      zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
    ), title = "3D PCA Decomposition")
  
  # new pca plot based on the new groups
  pca_plot_3d <- plot_ly(
    pca_df_3d, 
    x = ~PC1, y = ~PC2, z = ~PC3, 
    color = ~Group1, 
    symbol = ~Group2,
    # Creates a beautifully formatted custom hover label
    text = ~paste0(
      "<b>Sample ID:</b> ", Sample, "<br>",
      "<b>Batch/Group 1:</b> ", Group1, "<br>",
      "<b>Time/Group 2:</b> ", Group2, "<br>",
      "<b>ID Num:</b> ", Group3, "<br>",
      "<b>Suffix:</b> ", Group4
    ), 
    hoverinfo = "text",
    type = "scatter3d", 
    mode = "markers", 
    size = I(50)
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
        yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
        zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
      ), 
      title = "3D PCA Decomposition"
    )
  
  # Create the 3D UMAP Plot
  umap_plot_3d <- plot_ly(umap_df_3d, 
                          x = ~UMAP1, y = ~UMAP2, z = ~UMAP3, 
                          color = ~Group1, symbol = ~Group2,
                          text = ~Sample, hoverinfo = "text",
                          type = "scatter3d", mode = "markers", size = I(50)) %>%
    layout(title = "3D UMAP Manifold Projection")
  
  # Preview them
  pca_plot_3d
  umap_plot_3d
  
  # Save the clean data frames containing the coordinates and metadata
  saveRDS(pca_df_3d, "~/ShinyApps/ShinyReverb//pca_df_3d.rds")
  saveRDS(umap_df_3d, "~/ShinyApps/ShinyReverb//umap_df_3d.rds")
  saveRDS(var_exp, "~/ShinyApps/ShinyReverb//var_exp.rds")
  
  
  a
  
  
  
  
  # 2D plots
  # --- 1. Prepare the Data Table for Matrix Math ---
  mat <- as.matrix(tpm_matrix[, -1, with = FALSE])
  rownames(mat) <- tpm_matrix$Gene_ID
  
  # Log-transform TPM data (Crucial for PCA/UMAP to prevent extreme values from dominating)
  log_mat <- log2(mat + 1)
  
  
  # --- 2. Find the Most Variable Genes ---
  # Calculate variance for each gene across all samples
  gene_vars <- rowVars(log_mat)
  
  # Alternative: Use rowMads(log_mat) if you want a metric less sensitive to single outliers
  
  # Select the top N most variable genes (e.g., top 1000 or top 2000)
  top_n <- 1000
  variable_genes_idx <- order(gene_vars, decreasing = TRUE)[1:top_n]
  
  # Filter our matrix down to just those highly variable genes
  variable_mat <- log_mat[variable_genes_idx, ]
  
  
  # --- 3. Run Principal Component Analysis (PCA) ---
  # For PCA, dimensions are usually calculated across samples, so we transpose the matrix (t)
  # We center the data to mean 0; scaling is optional depending on your normalization
  pca_res <- prcomp(t(variable_mat), center = TRUE, scale. = FALSE)
  
  # Extract PCA coordinates into a data frame for plotting
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Sample <- rownames(pca_df)
  
  
  # --- 4. Run Uniform Manifold Approximation and Projection (UMAP) ---
  # Set a seed to make the UMAP layout reproducible
  set.seed(42)
  
  # Run UMAP on the transposed matrix (samples as rows, genes as columns)
  umap_config <- umap.defaults
  # If you have very few samples, you may need to reduce n_neighbors (default is 15)
  # umap_config$n_neighbors <- 5 
  
  umap_res <- umap(t(variable_mat), config = umap_config)
  
  # Extract UMAP coordinates into a data frame
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Sample <- rownames(umap_df)
  
  
  # --- 5. Optional: Plot the Results ---
  # Calculate Variance Explained for PCA Axis labels
  var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  
  # Plot PCA
  ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 3, color = "darkblue") +
    geom_text(vjust = -1, size = 3) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"), 
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "PCA of Top Highly Variable Genes") +
    theme_minimal()
  
  # Plot UMAP
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = Sample)) +
    geom_point(size = 3, color = "darkred") +
    geom_text(vjust = -1, size = 3) +
    labs(title = "UMAP of Top Highly Variable Genes") +
    theme_minimal()
  
  
  # color based on identifiers:
  library(tidyr)
  library(ggplot2)
  
  # ==========================================
  # 1. Extract Groups from Sample Names (PCA)
  # ==========================================
  pca_df$Sample <- rownames(pca_df)
  
  # Use separate() to split the Sample column by the underscores "_"
  # extra = "drop" ensures that the 3rd part (e.g., 21084D00) is ignored
  pca_df <- tidyr::separate(pca_df, Sample, 
                     into = c("Group1", "Group2"), 
                     sep = "_", 
                     remove = FALSE, 
                     extra = "drop")
  
  # ==========================================
  # 2. Extract Groups from Sample Names (UMAP)
  # ==========================================
  umap_df$Sample <- rownames(umap_df)
  
  umap_df <- tidyr::separate(umap_df, Sample, 
                      into = c("Group1", "Group2"), 
                      sep = "_", 
                      remove = FALSE, 
                      extra = "drop")
  
  
  # ==========================================
  # 3. Plot PCA Color-Coded by Groups
  # ==========================================
  # Option A: Color by the first part (e.g., RV09)
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group1)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"), 
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "PCA colored by Group 1 (First Part)",
         color = "Group 1") +
    theme_minimal()
  
  # Option B: Color by the first part AND change shape by the second part (e.g., 05)
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group1, shape = Group2)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"), 
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "PCA Dual-coded (Color & Shape)",
         color = "Group 1", 
         shape = "Group 2") +
    theme_minimal()
  
  
  # ==========================================
  # 4. Plot UMAP Color-Coded by Groups
  # ==========================================
  # Plotting UMAP using Group 1 for coloring
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group1)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(title = "UMAP colored by Group 1 (First Part)",
         color = "Group 1") +
    theme_minimal()
  
  # Plotting UMAP using Group 2 for coloring
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group2)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(title = "UMAP colored by Group 2 (Second Part)",
         color = "Group 2") +
    theme_minimal()
  
  
  
  
  
  
  # generating tpm matrix:
  # generating tpm matrix:
  # generating tpm matrix:
  
  #if(!require(data.table)) install.packages("data.table")
  library(data.table)
  
  setwd("/mnt/data/analysis/alireza/StringTie/")
  
  # 1. Read your sample list
  samples <- read.table("sample_list.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(samples) <- c("SampleID", "GTF_Path")
  
  tpm_list <- list()
  
  # 2. Process each file
  for(i in 1:nrow(samples)) {
    sample_id <- samples$SampleID[i]
    abund_file <- file.path("stringtie_output", sample_id, paste0(sample_id, "_gene_abund.txt"))
    
    if(file.exists(abund_file)) {
      # Read the file
      dt <- fread(abund_file, select = c("Gene Name", "TPM"))
      setnames(dt, c("Gene Name", "TPM"), c("Gene_ID", "TPM"))
      
      # CRITICAL FIX: Aggregate duplicate Gene IDs by summing their TPM values
      dt_clean <- dt[, .(TPM = sum(TPM)), by = Gene_ID]
      
      # Rename TPM column to the unique Sample ID
      setnames(dt_clean, "TPM", sample_id)
      
      tpm_list[[sample_id]] <- dt_clean
    } else {
      warning(paste("Missing file for sample:", sample_id))
    }
  }
  
  # 3. Merge the cleaned datasets (No duplicates = perfectly clean outer join)
  cat("Merging 66 cleaned samples...\n")
  final_tpm_matrix <- Reduce(function(x, y) merge(x, y, by="Gene_ID", all=TRUE), tpm_list)
  
  # 4. Fill structural missing values with 0
  final_tpm_matrix[is.na(final_tpm_matrix)] <- 0
  
  # Keep only genes that have a TPM > 0.5 in at least 5 or more samples
  keep <- rowSums(final_tpm_matrix[, -1] > 0.5) >= 5
  filtered_tpm_matrix <- final_tpm_matrix[keep, ]
  
  # 5. Save the clean matrix
  fwrite(filtered_tpm_matrix, "gene_tpm_matrix.csv")
  saveRDS(filtered_tpm_matrix, "final_tpm_matrix.RDS")
  cat("Success! Master TPM matrix written to gene_tpm_matrix.csv\n")
  
  tpm_matrix = readRDS("final_tpm_matrix.RDS")
  
  #QC
  all.equal(filtered_tpm_matrix, tpm_matrix) 
  all.equal(filtered_tpm_matrix, tpm_matrix, tolerance = 0)  
  
  
  
  
  
  
  
  
  


    # stringTie strand-specific parameter:
    setwd("/home/alireza/data/analysis/alireza/StringTie/temp")
    unstranded = data.table::fread("unstranded_test.txt")
    fr = data.table::fread("fr_test.txt")
    rf = data.table::fread("rf_test.txt")
    
    summary(unstranded$TPM)
    summary(fr$TPM)
    summary(rf$TPM)
    
    
    
    ### featureCounts
    # convert featureCounts output to TPM
    
    setwd("/mnt/data/analysis/alireza/featureCounts/")
    
    # 1. Load the featureCounts output
    counts_file <- read.table("final_counts.txt", header=TRUE, skip=1, row.names=1, check.names=FALSE)
    
    # 2. Extract gene lengths (convert from base pairs to kilobases)
    gene_lengths_kb <- counts_file$Length / 1000
    
    # 3. Extract just the raw count columns (columns 6 onwards)
    # (Columns 1-5 contain Chr, Start, End, Strand, Length)
    raw_counts <- counts_file[, 6:ncol(counts_file)]
    
    # 4. Calculate RPK (Reads Per Kilobase)
    rpk <- raw_counts / gene_lengths_kb
    
    # 5. Calculate per-sample scaling factors (sum of RPK / 1,000,000)
    scaling_factors <- colSums(rpk) / 1e6
    
    # 6. Calculate TPM by dividing RPK by scaling factors
    tpm <- t_tpm <- t(t(rpk) / scaling_factors)
    
    # 7. Save the TPM matrix to a file
    write.table(tpm, "tpm_counts.txt", sep="\t", quote=FALSE, col.names=NA)
    