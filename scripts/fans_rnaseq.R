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
library(ComplexHeatmap)
library(EnhancedVolcano)
library(openxlsx)
library(ggplot2)
library(ggrepel)

# source functions.R
source("./scripts/functions.R")

#######################################################
### Create mapping between transcript IDs & gene IDs ##
#######################################################

### Get annotation data from AnnotationHub ###
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
rownames(fans_counts) <- make.unique(rownames(fans_counts))
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


#########################
### Exploratory plots ###
#########################

dir.create("./plots")

# store raw library sizes
raw_libs <- colSums(y_raw$counts)

# make both plots with same y axis
y_range <- c(0, 1.4e7)

# raw lib sizes
png("plots/raw_library_sizes.png", width=2000, height=3000,  res=300)
lib_size_plt(raw_libs, y_range, title = "Raw library sizes")
dev.off()

# effective library sizes
png("plots/eff_library_sizes.png", width = 2000, height = 3000, res =  300)
lib_size_plt(eff_lib, y_range, title = "Effective library sizes")
dev.off()

# heatmap
png("plots/sample_heatmap.png", width = 3000, height = 2500, res=300)
sample_heatmap(y)
dev.off()


# mds plot with customized aesthetics
png("plots/MDSplot.png", width=2000, height = 2000, res = 300)
mds_plot(y)
dev.off()


#########################
### EdgeR: processing ###
#########################

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

###############################
#### Define Contrast Matrix ###
###############################

# compare expression of one nuc_type to the avg of the others
nuc.contrasts <- makeContrasts(
  NeuNvsAvg = NeuN - (PU1+SOX10+SOX2)/3,
  PU1vsAvg = PU1 -(NeuN+SOX10+SOX2)/3,
  SOX10vsAvg = SOX10 -(NeuN+PU1+SOX2)/3,
  SOX2vsAvg = SOX2 - (NeuN+PU1+SOX10)/3,
  levels = design.mat
  )
################################################
#### DE testing with quasi-liklihood F-tests ###
################################################

# these tests use a null hypothesis that DE is absent i.e. 0 FC
qlf.NeuN <- glmQLFTest(fit, contrast = nuc.contrasts[,"NeuNvsAvg"])
qlf.PU1 <- glmQLFTest(fit, contrast = nuc.contrasts[, "PU1vsAvg"])
qlf.SOX10 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX10vsAvg"] )
qlf.SOX2 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX2vsAvg"])

###########################
### DE testing glmTreat ###
###########################

# these examine if expression is signficantly different to a threshold 

treat.NeuN <- glmTreat(fit, contrast = nuc.contrasts[,"NeuNvsAvg"], lfc = log2(1.2))
treat.PU1 <- glmTreat(fit, contrast = nuc.contrasts[,"PU1vsAvg"], lfc = log2(1.2))
treat.SOX10 <- glmTreat(fit, contrast = nuc.contrasts[,"SOX10vsAvg"], lfc = log2(1.2))
treat.SOX2 <- glmTreat(fit, contrast = nuc.contrasts[,"SOX2vsAvg"], lfc = log2(1.2))

#########################################
###  generation of volcano plots ###
#########################################

# generate dataframes from the qlfTest objects
NeuN_qlfTest <- data.frame(qlf.NeuN$table)
PU1_qlfTest <- data.frame(qlf.PU1$table)
SOX10_qlfTest <- data.frame(qlf.SOX10$table)
SOX2_qlfTest <- data.frame(qlf.SOX2$table)

# Create vector of top genes to plot as labels
topNeuN <- rownames(topTags(qlf.NeuN, n = 25, sort.by = "PValue"))
topPU1 <- rownames(topTags(qlf.PU1, n = 25, sort.by = "PValue"))
topSOX10 <- rownames(topTags(qlf.SOX10, n = 25, sort.by = "PValue"))
topSOX2 <- rownames(topTags(qlf.SOX2, n = 25, sort.by = "PValue"))

# generate dataframes for the Treat test objects
NeuN_treat <- data.frame(treat.NeuN$table)
PU1_treat <-  data.frame(treat.NeuN$table)
SOX10_treat <- data.frame(treat.SOX10$table)
SOX2_treat <- data.frame(treat.PU1$table)


# set log2FC threshold to draw on plots.
threshold= 2


# Define marker genes
# vectors of canonical marker genes for each cell type
NeuN_mg <- c("Rbfox3", "Eno2", "Nefl", "Snap25", "Meg3", "Elavl3", "Tubb3", "Map2", "Mapt")
PU1_mg <- c("Spi1", "Aif1", "C1qa", "Tmem119", "Apoe", "Dock8", "Irf8")
SOX10_mg <- c("Sox10", "Mbp", "Mog", "Gal3st1", "S100b", "Pdgfra", "Cnp")
SOX2_mg <- c("Sox2", "Sox9", "Gfap", "Atp1b2", "Slc1a3", "Slc1a2")
multiple <- c("Prnp", "Apoe", "Stx6")

# Create a named list of markers by cell type:
markers <- list(neurons = NeuN_mg, microglia = PU1_mg, oligos = SOX10_mg, astrocytes = SOX2_mg, multiple = multiple)

NeuN_volcano <- PlotVolcano(NeuN_qlfTest,selectLab = topNeuN, FCcutoff = 2, title = "NeuN vs. Avg of Others")
PU1_volcano <- PlotVolcano(PU1_qlfTest, selectLab = topPU1, FCcutoff = 2, title = "PU1 vs. Avg of Others")
SOX10_volcano <- PlotVolcano(SOX10_qlfTest, selectLab = topSOX10, FCcutoff = 2, title = "SOX10 vs. Avg of Others")
SOX2_volcano <- PlotVolcano(SOX2_qlfTest, selectLab = topSOX2, FCcutoff = 2, title = "SOX2 vs. Avg of Others")


NeuN_marker_volcano <- MarkerGeneVolcano(NeuN_qlfTest, markers, title = "NeuN vs. Avg of Others")
PU1_marker_volcano <- MarkerGeneVolcano(PU1_qlfTest, markers, title = "PU1 vs. Avg of Others")
SOX10_marker_volcano <- MarkerGeneVolcano(SOX10_qlfTest, markers, title = "SOX10 vs. Avg of Others")
SOX2_marker_volcano <- MarkerGeneVolcano(SOX2_qlfTest, markers, title = "SOX2 vs. Avg of Others")

ggsave(filename = "plots/NeuN_marker_volcano.png", plot = NeuN_marker_volcano, width = 10, height = 10, dpi = 300)
ggsave(filename = "plots/PU1_marker_volcano.png", plot = PU1_marker_volcano, width = 10, height = 10, dpi = 300)
ggsave(filename = "plots/SOX10_marker_volcano.png", plot = SOX10_marker_volcano, width = 10, height = 10, dpi = 300)
ggsave(filename = "plots/SOX2_marker_volcano.png", plot = SOX2_marker_volcano, width = 10, height = 10, dpi = 300)

################################
#### Export DE analysis data ###
################################
# raw csvs
dir.create("output")
write.csv(x = NeuN_de, file = "output/FANS_pilot_NeuN_DE_genes.csv")
write.csv(x = PU1_de, file =  "output/FANS_pilot_PU1_DE_genes.csv")
write.csv(x = SOX10_de, file = "output/FANS_pilot_SOX10_DE_genes.csv")
write.csv(x = SOX2_de, file = "output/FANS_pilot_SOX2_DE_genes.csv")

# save same data as multi-sheet excel file
df_names <- c('Sheet1' = NeuN_de,
              'Sheet2' = PU1_de,
              'Sheet3' = SOX10_de,
              'Sheet4' = SOX2_de)
write.xlsx(df_names, file = "output/FANS_pilot_DE_genes.xlsx")