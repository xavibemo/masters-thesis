#!/bin/Rscript
## ****************************************************************************
## This script performs a Task Inferred from Differential Expression (TIDEs) 
## approach test to assess the significance of differential expression
## analysis (DEA) results obtained from RNA-Seq data. It takes as input the DEA results 
## and conducts permutation tests on gene expression data based on different experimental conditions.
##
## Get help: tides.R -h
## Run me: tides.R -n 1000 -o results/ -c Oxo --save_histrogram --use_all_genes
##
## AUTHOR: Xavier Benedicto Molina
## DATE: 18/12/23
##
## ****************************************************************************

## COMMAND LINE ARGUMENTS
library(optparse, quietly = TRUE)

option.list <- list(
  make_option(
    c("--n_permutations", "-n"), 
    action = "store",
    type = "numeric",
    default = 1000,
    help = "N of permutations for the permutation test (defaults to %default)"
  ),
  make_option(
    c("--condition", "-c"), 
    action = "store",
    type = "character",
    default = NA,
    help = "Experimental condition to be used must be a valid value between Oxo, PD, PI, PI_Oxo, and PI_PD"
  ),
  make_option(
    c("--out_dir", "-o"), 
    action = "store",
    type = "character",
    default = ".",
    help = "Main output directory"
  ),
  make_option(
    c("--save_histogram"), 
    action = "store_true",
    type = "logical",
    default = FALSE,
    help = "Whether to plot out the histograms for the permutation tests"
  ),
  make_option(
    c("--use_all_genes"), 
    action = "store_true",
    type = "logical",
    default = FALSE,
    help = "Whether to use all genes that comprise a metabolic taks or just the essential ones"
  )
)
# Debug
# opt <- parse_args(OptionParser(option_list = option.list), args = c("--n_permutations=10000", "--condition=PD", "--out_dir=test"))
opt <- parse_args(OptionParser(option_list = option.list))


## LIBRARIES
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
library(foreach, quietly = TRUE, warn.conflicts = FALSE)
library(doParallel, quietly = TRUE, warn.conflicts = FALSE)


## FUNCTIONS
save_histogram <- function(permutations, results, task, filename) {
  
  gg <- ggplot(permutations, aes(x = mean)) +
    ggtitle(paste0(treatment, " - ", task)) +
    ylab("") +
    xlab("Mean") +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = mean(results$log2FoldChange), color = "red3", linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(filename, plot = gg, width = 4, height = 4)
}


## MAIN
# Creating our dir
print(paste0("Setting out directory to ", opt$out_dir))
dir.create(opt$out_dir, showWarnings = FALSE)

# Reading in DEA results
print("Loading DEA results...")
all.results <- readRDS("data/all.results.shrinkedLFC.RDS")
print("Done!")

if (opt$save_histogram) {
  print(paste0("Histogram plots will be saved to ", opt$out_dir, "/histograms/"))
  dir.create(paste0(opt$out_dir, "/histograms/"), showWarnings = FALSE)
}

if (!opt$use_all_genes) {
  print("Permutation test will be performed only with essential genes!")
  essential.genes <- readRDS("data/gene.to.tasks.long.RDS")
  all.tasks <- levels(essential.genes$Task)
} else {
  print("Permutation test will be performed with all genes!")
  all.genes <- readRDS("data/pfba_genes_to_reactions_to_tasks.RDS")
  all.tasks <- levels(all.genes$task.id)
}

print(paste0("Starting permutation test with n = ", opt$n_permutations))

# Setting levels
if (opt$condition %in% levels(all.results$Treatment)) {
  treatment <- opt$condition
} else {
  stop(paste0("Treatment condition ", opt$condition, " is not valid!"))
}

print(paste0("Starting permutation test for:    ", treatment))
lfc.vector <- all.results %>% filter(Treatment == treatment) %>% pull(log2FoldChange)

# Setting up parallesisation
cl <- makeCluster(detectCores())
registerDoParallel(cl = cl)

# Main loop
final.results <- foreach(
  task = all.tasks, 
  .combine = rbind,
  .packages = c("tidyverse")
) %dopar% {
  
  if (!opt$use_all_genes) {
    results <- essential.genes %>% filter(Task == task) %>% merge(all.results %>% filter(Treatment == treatment), by = "Symbols")
  } else {
    results <- all.genes %>% filter(task.id == task) %>% merge(all.results %>% filter(Treatment == treatment), by = "Symbols")
  }
  
  lfc.mean <- mean(results$log2FoldChange)
  lfc.sd <- sd(results$log2FoldChange)
  permutations <- tibble(treatment.id = character(), task.id = character(), mean = numeric())
  
  # Inner loop using purrr::map_dbl
  permutations <- permutations %>% add_row(
    treatment.id = treatment, 
    task.id = task, 
    mean = map_dbl(1:opt$n_permutations, ~{sample(lfc.vector, size = nrow(results)) %>% mean()})
  )
  
  # Computing pvalues
  pvalue <- min(
    sum(permutations$mean >= mean(results$log2FoldChange)) / opt$n_permutations, 
    sum(permutations$mean <= mean(results$log2FoldChange)) / opt$n_permutations
  )
  
  # If selected, saving plots into new directory
  if (opt$save_histogram) {
    filename <- paste0(opt$out_dir, "/histograms/", treatment, "_", task, "_histogram.pdf")
    save_histogram(permutations, results, task, filename)
  }
  
  tmp.results <- tibble(
    treatment.id = treatment, 
    task.id = task, 
    observed.mean = lfc.mean, 
    observed.dispersion = lfc.sd,
    random.mean = mean(permutations$mean), 
    n.genes = nrow(results),
    pvalue = pvalue
  )
  return(tmp.results)
}

# Stop paralelization
stopImplicitCluster()

print("Permutation test completed!")
print("Saving results...")

filename <- paste0(opt$out_dir, "/permutation_test_results_", treatment)
write.table(final.results %>% arrange(pvalue), file = paste0(filename, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
saveRDS(final.results %>% arrange(pvalue), file = paste0(filename, ".RDS"))

print("Done!")
