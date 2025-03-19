# Script to import Salmon quantification data from FANS pilot experiment 
# using Tximport. Followed by DE analysis using edgeR.

# Author: T. Murphy
# Date: 2025-03-04

# load required libraries. 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(tximport)
library(AnnotationHub)
library(ensembldb)
library(stringr)
library(RColorBrewer)
library(edgeR)
library(Glimma)
library(ComplexHeatmap)
library(EnhancedVolcano)


### Getting annotation data from AnnotationHub ###
# initialise annotation hub
hub <- AnnotationHub()

# generate a query for latest ensembl mouse annotations
ensdb_query <- query(hub,c("EnsDb", "musculus", "112"))
ensdb_query

# from the previous query, select most recent EnsDb object
ensdb_112 <- ensdb_query[['AH116909']]

# extract the transcript data from the EnsDb object
tx_data <- transcripts(ensdb_112, return.type = "DataFrame")
tx_data

# slice the tx_data to include only the tx_id and gene_id columns
# N.B. human readable gene name in column tx_external_name
tx_to_gene <- tx_data[, c("tx_id", "gene_id")]


###############################################################
### Importing the Salmon quantification data using tximport ###
###############################################################

# get the salmon quant files.
quants_dir <- "salmon_quant/"

quants_files <- list.files(quants_dir,
                           pattern = "_quant.sf$",
                           recursive = TRUE,
                           full.names = TRUE)

# check number of quant.sf files == number we expect
length(quants_files)

# use the file names to make sample names from which we can make metadata
sample_names <- gsub("_S._quant.sf", "", basename(quants_files))
sample_names <- gsub("_S.._quant.sf", "", sample_names) # for samples 10-12

# performing the data import using tximport
names(quants_files) <- sample_names

txi <- tximport(
  files = quants_files,
  type = "salmon",
  tx2gene = tx_to_gene,
  ignoreTxVersion = TRUE
)
#############################################
### replace ensembl_IDs with gene symbols ###
#############################################

# get the IDs present in the count matrix
ensembl_ids <- rownames(txi$counts)

# find gene symbols in our EnsDb object
gene_mapping <- genes(ensdb_112,
                      filter = GeneIdFilter(ensembl_ids),
                      columns = c("gene_id", "gene_name"))
# convert to df
gene_mapping_df <- as.data.frame(gene_mapping)

# ensure no rows are missing a value in the gene_name column 
sum(is.na(gene_mapping_df$gene_name))

# We reorder the gene_mapping_df to match the order of the ensembl_ids
# N.B ensembl ids preserves the order of our data
ids2symbols <- gene_mapping_df$gene_name[match(ensembl_ids, gene_mapping_df$gene_id)]

# replace the ensembl ids with gene symbols in the txi object
rownames(txi$counts) <- ids2symbols
rownames(txi$length) <- ids2symbols
rownames(txi$abundance) <- ids2symbols

##############################################################
### Generating sample metadata and design matrix for edgeR ###
##############################################################

# creating the sample metadata
sample_meta <- data.frame(
  sample_name = sample_names,
  animal_no = str_extract(sample_names, "^[0-9]+"),
  nuc_type = str_extract(sample_names, "(?<=_)[A-Za-z0-9]+"),
  extraction_batch = rep(c(1, 2), each = 8)
)

sample_meta

# recode variables as factors
sample_meta$animal_no <- as.factor(sample_meta$animal_no)
sample_meta$nuc_type <- as.factor(sample_meta$nuc_type)
sample_meta$extraction_batch <- as.factor(sample_meta$extraction_batch)

############################################################
### edgeR DEGList object creation with an offset matrix  ###
### to normalise for library size and transcript length. ###
############################################################

# create the raw count matrix 
fans_counts <- txi$counts
y_raw <- DGEList(fans_counts)
write.csv(fans_counts, "fans_pilot_raw_counts.csv")

# handle metadata
ymeta <- sample_meta[, colnames(sample_meta) != "sample_name"]
y_raw$samples <- cbind(y_raw$samples, ymeta)
y_raw$samples$group <- NULL

# create matrix of transcript lengths from tximport
tx_lens <- txi$length

# generate a length-normalisation matrix
norm_mat <- tx_lens / exp(rowMeans(log(tx_lens)))
norm_counts <- txi$counts / norm_mat 

# calculate effective library sizes
eff_lib <- calcNormFactors(norm_counts) * colSums(norm_counts)

# multiply length-normalised counts column-wise by effective library sizes
norm_mat <- sweep(norm_mat, 2, eff_lib, "*")
norm_mat <- log(norm_mat)

# apply norm_mat as offsets 
y <- y_raw
y <- scaleOffset(y, norm_mat) 


##################
### plotting 1 ###
##################

dir.create("./plots")

# plotting library sizes
lib_size_plt <- function(lib_sizes, y_range, title){
  
  par(mar = c(8, 6, 4, 2), mgp = c(5, 1, 0)) # set margins: bottom, left, top, right
  
  barplot(lib_sizes,
          main = title,
          las = 2,
          xlab = NULL,
          ylim = y_range,
          ylab = "Library Size")

}

# store raw library sizes
raw_libs <- colSums(y_raw$counts)

# make both plots with same y axis
y_range <- c(0, 1.4e7)

png("plots/raw_library_sizes.png", width=2000, height=3000,  res=300)
lib_size_plt(raw_libs, y_range, title = "Raw library sizes")
dev.off()

png("plots/eff_library_sizes.png", width = 2000, height = 3000, res =  300)
lib_size_plt(eff_lib, y_range, title = "Effective library sizes")
dev.off()


# correlations between samples

# heatmap of correlation coefficients
sample_heatmap <- function(y){
  corr <- cor(y$counts)
  pheatmap(corr)
}

png("raw_library_sizes.png", width = 3000, height = 2500, res=300)
sample_heatmap()
dev.off()

# MDS plot function
mds_plot <- function(y, method="logFC"){

  
  # Define colours for each type of our types
  dark_cols <- brewer.pal(4, "Dark2")
  nuc_types <- levels(as.factor(y$samples$nuc_type))
  nuc_colors <- setNames(dark_cols, nuc_types)

  
  # Plot the MDS using coloured points (pch=16 for solid circles) without text labels
  plotMDS.DGEList(y,
                  bg = nuc_colors,
                  pch = 21,
                  cex = 1.5,
                  method = method,
                  main = "MDS plot of libraries from disitnct nuclei populations", labels = NULL)

  # Add a legend in the top right
  legend("bottomright",
         legend = nuc_types,
         pt.bg = nuc_colors,
         pch = 21,
         cex = 1.5,
         title = "Nuclear Type")
}

mds_plot(y)


#########################
### EdgeR: processing ###
#########################

# low count filtering

# filter out lowly expressed genes
keep <-filterByExpr(y, group = y$samples$nuc_type)
table(keep)
y <- y[keep, ]

# Dispersion estimation & model fitting 

# specify design matrix where each type of nuclei is represented as a column
design.mat <- model.matrix(~0 + nuc_type + animal_no, data = sample_meta)
colnames(design.mat) <- gsub("nuc_type", "", colnames(design.mat))

# estimate dispersion parameters 
y <- estimateDisp(y, design.mat)

# fit gene-wise quasi negative binomial models
fit <- glmQLFit(y, design = design.mat)

#### DE testing with quasi-liklihood F-tests ###
# We compare expression of one nuc_type versus the avg of the others

nuc.contrasts <- makeContrasts(
  NeuNvsAvg = NeuN - (PU1+SOX10+SOX2)/3,
  PU1vsAvg = PU1 -(NeuN+SOX10+SOX2)/3,
  SOX10vsAvg = SOX10 -(NeuN+PU1+SOX2)/3,
  SOX2vsAvg = SOX2 - (NeuN+PU1+SOX10)/3,
  levels = design.mat
  )

qlf.NeuN <- glmQLFTest(fit, contrast = nuc.contrasts[,"NeuNvsAvg"])
topTags(qlf.NeuN)

qlf.PU1 <- glmQLFTest(fit, contrast = nuc.contrasts[, "PU1vsAvg"])
topTags(qlf.PU1)

qlf.SOX10 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX10vsAvg"] )
topTags(qlf.SOX10)

qlf.SOX2 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX2vsAvg"])
topTags(qlf.SOX2)

#################################
### Plotting 2: Volcano plots ###
#################################

plot_volcano <- function(DEgenes_df){
  # Depends on enhanced volcano package
  # configured to work with ouput from EdgeR's qlfTest function
  # Will return the plot as an object with the name taken from the first part
  # of the object supplied to the function. 
  
  plt <- EnhancedVolcano(DEgenes_df,
                  lab = rownames(DEgenes_df),
                  x = "logFC",
                  y = "PValue")
  
  return(plt)
}

NeuN_genes <- data.frame(qlf.NeuN$table)
NeuN_volcano <- plot_volcano(NeuN_genes)
