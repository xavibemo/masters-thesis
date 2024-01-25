## UTILS

library(ggsignif)


# Custom function to map p.values to "*"
map_pvalues_to_asterisks <- function(
    p_values,
    alpha_levels = c(0.05, 0.01, 0.001, 0.0001),
    symbols = c("*", "**", "***", "****")
) {
  
  # Managing NAs
  p_values[is.na(p_values)] <- 1
  
  # Checking if p_values are numeric and between (0, 1)
  if (!is.numeric(p_values) || any(p_values < 0 | p_values > 1)) {
    stop("Invalid p-values. Please provide a numeric vector of p-values between 0 and 1.")
  }
  # Checking if alpha levels are numeric and between (0, 1)
  if (!is.numeric(alpha_levels) || any(alpha_levels < 0 | alpha_levels > 1)) {
    stop("Invalid alpha levels. Please provide a numeric vector of alpha levels between 0 and 1.")
  }
  if (length(alpha_levels) != length(symbols)) {
    stop("The lengths of alpha_levels and symbols must be the same.")
  }
  
  symbols_to_use <- rep("", length(p_values))
  
  for (i in seq_along(alpha_levels)) {
    symbols_to_use[p_values <= alpha_levels[i]] <- symbols[i]
  }
  
  return(symbols_to_use)
}


# Custom function to plot TPMs as bars using ggplot
TPMplot <- function(longTPM, gene, LFCvalues, 
  gray_color = TRUE, 
  minimal_theme = TRUE, 
  show_comparisons = TRUE,
  pval = 0.05, 
  lfc = log2(1.5)
) {
  
  # Creating all pairwise comparisons
  comparisons <- data.frame(
    V1 = rep("Ctr", 5), 
    V2 = c("Oxo", "PD", "PI", "PI_Oxo", "PI_PD")
  ) %>% t() %>% as.data.frame() %>% as.list() 
  
  # Keeping only significative comparisons
  rownames(LFCvalues) <- LFCvalues$treatment.id
  LFCvalues$log2FoldChange[is.na(LFCvalues$log2FoldChange)] <- 0
  keep <- LFCvalues$padj < pval & abs(LFCvalues$log2FoldChange) > lfc
  comparisons <- comparisons[keep]
  LFCvalues <- LFCvalues[keep,]
  
  # Summary to display data in bars
  longTPM.summary <- longTPM %>%
    group_by(treatment.id) %>%
    summarise(
      sd = sd(tpm),
      tpm = median(tpm)
    )
  
  # Setting color scale
  if (gray_color) color.scale <- gray.colors(6)
  else color.scale <- c("#ff5349" ,"#f79256", "#fbd1a2", "#7dcfb6", "#00b2ca", "#1d4e89")
  
  # Plot
  gg <- ggplot(longTPM, aes(x = treatment.id, y = tpm)) +
    
    # ggtitle(gene) + #, subtitle = "Log2FC and adj. p-values are indicated") +
    
    geom_col(data = longTPM.summary, width = 0.8, fill = color.scale) +
    geom_errorbar(data = longTPM.summary, aes(ymax = tpm+sd, ymin = tpm-sd), width = 0.25, linewidth = 0.4, color = "black") +
    geom_jitter(position = position_jitter(width = 0.1, height = 0), shape = 1, size = 1, color = "black") + 
    
    xlab("\nTreatment") +
    ylab("Transcripts per Million (TPM)\n")
  
  if (minimal_theme) {
    gg <- gg + theme_minimal() +
      theme(
        legend.position = "none",
        axis.line.x.bottom = element_line(linewidth = 0.5),
        axis.line.y.left = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size = 8),
        # plot.title = element_text(face = "bold", size = 8)
      )
  } else {
    gg <- gg + theme_bw() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 8),
        # plot.title = element_text(face = "bold", size = 8)
      )
  }
    
  if (length(comparisons) > 0 & show_comparisons) {
    gg <- gg + geom_signif(
      comparisons = comparisons, 
      annotation = map_pvalues_to_asterisks(LFCvalues$padj),
      tip_length = 0.05, 
      step_increase = 0.1, 
      size = 0.3,
      textsize = 3, 
      margin_top = 0.1,
      vjust = 0.4,
      # alpha = .5
    )
  }
  
  gg
}
