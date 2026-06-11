

   # immunarch
    
#mkdir -p /mnt/data/analysis/trust4/immunarch_input
#cd /mnt/data/analysis/trust4/
#  for f in *_airr.tsv; do
#     SAMPLE=$(basename "$f" _trust4_airr.tsv)
#     cp "$f" immunarch_input/${SAMPLE}.tsv
#   done

library(immunarch)
library(ggpubr)
suppressWarnings(vis(tc))

  # Load all files in the directory
  immdata <- repLoad("/mnt/data/analysis/alireza/trust4/immunarch_input/", .mode = "AIRR")
    #If repLoad doesn't auto-detect AIRR format:
    #immdata <- repLoad("/mnt/data/analysis/trust4/immunarch_input/", .format = "airr")
  saveRDS(immdata, "~/ShinyApps/ShinyReverb/immdata.RDS")
  
  # Check what loaded
  names(immdata$data)
  head(immdata$data[[1]])
  
  
  
  #### QC
  # Chain composition per sample
  lapply(immdata$data, function(df) {
    table(substr(df$V.name, 1, 3))
  })
  
  # Strip alleles for cleaner gene usage
  immdata$data <- lapply(immdata$data, function(df) {
    df$V.name <- sub("\\*.*", "", df$V.name)
    df$J.name <- sub("\\*.*", "", df$J.name)
    df$D.name <- sub("\\*.*", "", df$D.name)
    df
  })
  
  # Split by chain if needed
  library(dplyr)
  igk_data <- lapply(immdata$data, function(df) {
    df %>% filter(substr(V.name, 1, 3) == "IGK")
  })
  
  
  # adding metadata
  # immdata$meta should have a "Sample" column matching your filenames
  library(stringr)
  library(dplyr)
  
  sample_names <- names(immdata$data)
  
  parsed <- stringr::str_match(sample_names, "^(RV\\d+)_(\\d+)_(\\d+)D(.*)$")
  
  immdata$meta <- data.frame(
    Sample  = sample_names,
    Group1  = parsed[,2],
    Group2  = parsed[,3],
    Group3  = parsed[,4],
    Timepoint_raw = parsed[,5],
    stringsAsFactors = FALSE
  )
  
  # Fix the known mislabeled entries
  immdata$meta$Timepoint <- recode(immdata$meta$Timepoint_raw,
                                   "_1" = "01",
                                   "_9" = "09",
                                   .default = immdata$meta$Timepoint_raw
  )
  
  #immdata$meta
  
  immdata$data <- lapply(immdata$data, function(df) {
    df$Proportion <- df$Clones / sum(df$Clones, na.rm = TRUE)
    df
  })
  
  # Verify
  head(immdata$data[[1]][, c("Clones", "Proportion")])
  
  # Repertoire diversity:
  div_chao <- repDiversity(immdata$data, .method = "chao1")
  vis(div_chao)
  
  div_shannon <- repDiversity(immdata$data, .method = "div")
  vis(div_shannon)
  
  #Clonal proportions (top clones):
  clon <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000))
  vis(clon)
  
  # V/J gene usage:
  #gene_usage <- geneUsage(immdata$data, "hs.trbv")  # this is good for human
  gene_usage <- immdata$data[[1]] %>%
    count(V.name, sort = TRUE)
  #vis(gene_usage)
  
  ggplot(gene_usage, aes(x = reorder(V.name, -n), y = n)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(x = "V gene", y = "Count", title = "V gene usage")
  
  #If there are many V genes and the x-axis is too crowded, show only the top N:
  top_genes <- gene_usage %>% slice_max(n, n = 20)
  
  ggplot(top_genes, aes(x = reorder(V.name, -n), y = n)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "V gene", y = "Count", title = "Top 20 V genes")
  
  
  #For comparing gene usage across all samples at once, build a combined data frame:
  gene_usage_all <- bind_rows(lapply(names(immdata$data), function(s) {
    immdata$data[[s]] %>%
      count(V.name, name = "n") %>%
      mutate(Sample = s)
  }))
  
  ggplot(gene_usage_all, aes(x = reorder(V.name, -n), y = n, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(x = "V gene", y = "Count", title = "V gene usage by sample")
  
  # This is messy with many genes — consider faceting instead of dodging:
  ggplot(gene_usage_all, aes(x = reorder(V.name, -n), y = n)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Sample, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
    labs(x = "V gene", y = "Count")
  
  
  
  # Clonotype tracking across samples:
  tc <- trackClonotypes(immdata$data, list(1, 10), .col = "aa")
  vis(tc)
  
  tc$CDR3.aa
  clonotypes_of_interest <- c("YCGQATHLPPTF", "YCGQGTHFPPTF")
  tc <- trackClonotypes(immdata$data, clonotypes_of_interest, .col = "aa")
  vis(tc)
  
  
  #1. Gene usage functions assume human gene names (hs.trbv, hs.ighv, etc.). Since you used IMGT macaque gene names (e.g. IGKV2-65), built-in geneUsage reference panels won't match well. You'll likely need to build gene usage manually:
  library(dplyr)
  
  gene_usage_manual <- immdata$data[[1]] %>%
    count(V.name, sort = TRUE)
  
  ggplot(gene_usage_manual, aes(x = reorder(V.name, -n), y = n)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  
  #Multiple chain types in one file — your TRUST4 output mixes IGH, IGK, IGL, TRA, TRB etc. You may want to split by chain before analysis:
  library(dplyr)
  immdata$data <- lapply(immdata$data, function(df) {
    df$chain <- substr(df$V.name, 1, 3)
    df
  })
  
  # Analyze IGK only, for example

  # Allele suffixes (*01, *02) — immunarch usually handles these fine, but if gene usage plots look fragmented, strip alleles:
  # Fix typo - second line should strip J.name, not V.name again
  immdata$data <- lapply(immdata$data, function(df) {
    df$V.name <- sub("\\*.*", "", df$V.name)
    df$J.name <- sub("\\*.*", "", df$J.name)
    df$D.name <- sub("\\*.*", "", df$D.name)
    df
  })
  
  # Note: do allele-stripping BEFORE subsetting by chain, 
  # or re-subset igk_data after stripping
  igk_data <- lapply(immdata$data, function(df) df %>% filter(chain == "IGK"))

  #Now visualize the IGK subset
  #Clonality within IGK only
  igk_clon <- repClonality(igk_data, .method = "top", .head = c(10, 100, 1000))
  vis(igk_clon)  
  
  #Diversity within IGK only
  igk_div <- repDiversity(igk_data, .method = "div")
  vis(igk_div)

  #3. IGK V gene usage per sample
  igk_vgene <- bind_rows(lapply(names(igk_data), function(s) {
    igk_data[[s]] %>%
      count(V.name, name = "n") %>%
      mutate(Sample = s, Freq = n / sum(n))
  }))
  
  ggplot(igk_vgene, aes(x = reorder(V.name, -Freq), y = Freq, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(x = "IGKV gene", y = "Frequency", title = "IGK V gene usage")  

  #4. IGK V-J pairing heatmap (per sample)
  #This shows which V and J genes are used together — useful for spotting biased recombination:
  for (s in names(igk_data)) {
    vj_table <- igk_data[[s]] %>%
      count(V.name, J.name) %>%
      tidyr::complete(V.name, J.name, fill = list(n = 0))
    
    p <- ggplot(vj_table, aes(x = J.name, y = V.name, fill = n)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkred") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6)) +
      labs(title = paste("IGK V-J pairing -", s))
    
    print(p)
  } 

  #5. CDR3 length distribution (IGK)
  #CDR3 length distribution is a classic repertoire QC/comparison metric:  
  cdr3_len <- bind_rows(lapply(names(igk_data), function(s) {
    igk_data[[s]] %>%
      mutate(len = nchar(CDR3.aa), Sample = s)
  }))
  
  ggplot(cdr3_len, aes(x = len, fill = Sample)) +
    geom_density(alpha = 0.4) +
    labs(x = "CDR3 length (aa)", y = "Density", title = "IGK CDR3 length distribution")

  #6. Track top IGK clonotypes across samples (e.g. timepoints)
  igk_tc <- trackClonotypes(igk_data, list(1, 10), .col = "aa")
  vis(igk_tc)

  #Interpretation notes
  #Clonality plots (repClonality) show what fraction of the repertoire is occupied by the top N clones — a shift toward higher clonality at later timepoints can indicate clonal expansion (e.g. in response to vaccination/infection).
  #Diversity (repDiversity) — declining diversity over time often pairs with increasing clonality, both pointing to a focused immune response.
  #V/J usage shifts between timepoints can indicate selection of particular gene segments during an immune response.
  #CDR3 length distribution shifts (e.g. shorter CDR3s becoming more common) can reflect selection pressure on the repertoire.
  #Clonotype tracking is the most direct way to see if specific clones expand from baseline (Day0) through later timepoints (Day1, Day9) — a hallmark of an antigen-specific response.
  
  
  