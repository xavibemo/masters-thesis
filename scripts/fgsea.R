#!/bin/Rscript
## ****************************************************************************
## Script to run fast GSEA (fgsea) for all the DEA results across all conditions
## and to draw corresponding heatmaps with the results.
##
## Get help: gsea.R -h
## Run me: Rscript fgsea.R --permutations=1e5 --gmt_file=pathways.gmt --nproc=24
##
## AUHTOR: Xavier Benedicto Molina
## DATE: 3/11/23
##
## ****************************************************************************

## COMMAND LINE ARGUMENTS
# Getting parsed arguments if any
library(optparse, quietly = TRUE)

option.list <- list(
  make_option(
    c("--ranking_file", "-r"), 
    action = "store",
    type = "character",
    default = NA,
    help = ".tsv or .RDS file containing a predefined GSEA ranking (first column must be genes, second column ranking metric)"
  ),
  make_option(
    c("--permutations", "-p"), 
    action = "store",
    type = "numeric",
    default = NA,
    help = "N of permutations for fgsea"
  ),
  make_option(
    c("--gmt_file"), 
    action = "store",
    type = "character",
    default = NA,
    help = ".gmt file for fgsea"
  ),
  make_option(
    c("--nproc"), 
    action = "store",
    type = "numeric",
    default = NA,
    help = "n of tasks"
  )
)
opt <- parse_args(OptionParser(option_list = option.list))


## LIBRARIES
suppressPackageStartupMessages(library(tidyverse))
library(fgsea, quietly = TRUE)


## FUNCTIONS
# Function to convert a matrix to a list (obtained from https://biostatsquid.com/gene-set-enrichment-analysis/)
matrix_to_list <- function(pws) {
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Function to prepare a gmt_file and the backgroung genes for the hyperbolic test 
# (obtained from https://biostatsquid.com/gene-set-enrichment-analysis/)
prepareGMT <- function(gmt_file, genes_in_data) 
{
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  
  for (i in 1:dim(mat)[2]) {
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  # Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  # filter for gene sets with more than 5 genes annotated
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] 
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  return(final_list)
}


## MAIN
# Creating out dir
out_dir <- unlist(str_split(opt$gmt_file, "\\.|/"))[1:4] %>% paste(collapse = "_")
print(paste0("Creating out directory ", out_dir))
dir.create(out_dir, showWarnings = FALSE)

# Reading in RDS files
if (opt$ranking_file %>% endsWith(".RDS")) {
  ranking <- readRDS(opt$ranking_file)
} else if (opt$ranking_file %>% endsWith(".tsv")) {
  ranking <- read_tsv(opt$ranking_file, show_col_types = FALSE)
} else {
  stop("File format is not .tsv or .RDS!")
}

# Checking number of cols
if (!(ranking %>% ncol == 2)) {
  stop("Ranking file contains incorrect number of columns!")
}

print("Loading ranking file. All is well!")
print(paste0("Will use reference: ", opt$gmt_file))

# fgsea call
cols <- colnames(ranking)
genes <- ranking %>% pull(cols[1])
bg.genes <- prepareGMT(opt$gmt_file, genes)
lfc.ranking <- setNames(ranking %>% pull(cols[2]), genes)

print("Ranking and background correctly references. Starting GSEA call...")
  
GSEAres <- fgseaMultilevel(
  pathways = bg.genes,
  stats = lfc.ranking,
  scoreType = "std",
  minSize = 15,
  maxSize = 500,
  nproc = opt$nproc,
  nPermSimple = opt$permutations
)

# Formatting leading edge
GSEAres <- GSEAres %>% mutate(
  Change = ifelse(NES > 0, "Upregulated", "Downregulated"),
  leadingEdge = sapply(leadingEdge, function(i) paste(i, collapse = "; "))
)

print(paste0("Writing results table at: ", out_dir, "/fgsea_results.tsv"))

GSEAres %>% arrange(padj, NES) %>%
  write.table(paste0(out_dir, "/fgsea_results.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  
print("Creating NES density plots.")

density <- ggplot(GSEAres, aes(x = abs(NES))) + 
  geom_density(fill = "royalblue", alpha = 0.1, adjust = 1/2) + 
  theme_bw() + 
  xlim(c(0, 3)) +
  ylim(c(0, 2.5))
  
ggsave(paste0(out_dir, "/fgsea_NESdensity.pdf"), plot = density, height = 4, width = 4)

print("Done!")
