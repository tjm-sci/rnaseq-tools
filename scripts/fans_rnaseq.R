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

sample_heatmap <- function(y){
  # heatmap of correlation coefficients
  corr <- cor(y$counts)
  pheatmap(corr)
}

png("plots/sample_heatmap.png", width = 3000, height = 2500, res=300)
sample_heatmap(y)
dev.off()

# MDS plot function
mds_plot <- function(y){
  
  # wrapper function for EdgeR's plotMDS with nice
  # formatting for different types of nuclei 
  
  # Define colours for each type of our types
  dark_cols <- brewer.pal(4, "Dark2")
  nuc_types <- levels(as.factor(y$samples$nuc_type))
  nuc_colors <- setNames(dark_cols, nuc_types)

  
  # Plot the MDS using coloured points (pch=16 for solid circles) 
  plotMDS.DGEList(y,
                  bg = nuc_colors,
                  pch = 21,
                  cex = 1.5,
                  method = "logFC",
                  main = "MDS plot of libraries from disitnct nuclei populations", labels = NULL)

  # Add a legend in the bottom right
  legend("bottomright",
         legend = nuc_types,
         pt.bg = nuc_colors,
         pch = 21,
         cex = 1.5,
         title = "Nuclear Type")
}

mds <- mds_plot(y)
png("plots/MDSplot.png", width=2000, height = 2000, res = 300)
mds_plot(y)
dev.off()


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
y$
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
topTags(qlf.NeuN)
qlf.PU1 <- glmQLFTest(fit, contrast = nuc.contrasts[, "PU1vsAvg"])
topTags(qlf.PU1)
qlf.SOX10 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX10vsAvg"] )
topTags(qlf.SOX10)
qlf.SOX2 <- glmQLFTest(fit, contrast = nuc.contrasts[, "SOX2vsAvg"])
topTags(qlf.SOX2)

###########################
### DE testing glmTreat ###
###########################

# these examine if expresion is signficantly different to a threshold ###

treat.NeuN <- glmTreat(fit, contrast = nuc.contrasts[,"NeuNvsAvg"], lfc = log2(1.2))
topTags(treat.NeuN)
treat.PU1 <- glmTreat(fit, contrast = nuc.contrasts[,"PU1vsAvg"], lfc = log2(1.2))
topTags(treat.PU1)
treat.SOX10 <- glmTreat(fit, contrast = nuc.contrasts[,"SOX10vsAvg"], lfc = log2(1.2))
topTags(treat.SOX10)
treat.SOX2 <- glmTreat(fit, contrast = nuc.contrasts[,"SOX2vsAvg"], lfc = log2(1.2))
topTags(treat.SOX2)

#########################################
### Enhanced Volcano Plot ###
#########################################


PlotVolcano <- function(DEgenes_df, title = ""){
  # DEgenes_df: data frame of differential expression results (e.g., from edgeR),
  # Create a vector of marker genes (the union of all marker lists) to emphasize with labels.
  select_labels <- unique(unlist(markers))
  
  # Create the base volcano plot using EnhancedVolcano (with its default point colors).
  plt <- EnhancedVolcano(
    DEgenes_df,
    lab = DEgenes_df$gene,
    x = "logFC",
    y = "PValue",
    title = title,
    pointSize = 0.75,
    subtitle = NULL,
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    legendPosition = "right",
    selectLab = select_labels
  )
  
  print(plt)
  return(plt)
}

#################################################################
### Custom Function for plotting marker genes on Volcano Plot ###
#################################################################

MarkerGeneVolcano <- function(DEgenes_df, markers, title = ""){
  # DEgenes_df: data frame of DE results (e.g., from edgeR)
  #            Must contain columns "logFC" and "PValue" and have rownames as gene symbols.
  # markers: a named list of marker gene vectors, one element per cell type,
  #          e.g. list(neurons = c("Rbfox3", "Eno2", "Nefl", "Snap25")...)
  # title: Plot title
  
  # Store gene names in a column
  DEgenes_df$gene <- rownames(DEgenes_df)
  
  # Create a new column for cell type annotation, defaulting to NA.
  DEgenes_df$cellType <- NA
  
  # Loop over the markers list and assign cell type label for genes found in each marker vector.
  for(ct in names(markers)){
    DEgenes_df$cellType[DEgenes_df$gene %in% markers[[ct]]] <- ct
  }
  
  # Make cellType a factor and specify levels to match colouring with other plots
  DEgenes_df$cellType <- factor(DEgenes_df$cellType, levels = names(markers))
  
  # Create a vector of marker genes (the union of all marker lists) to emphasize with labels.
  select_labels <- unique(unlist(markers))
  
  # Create the base volcano plot with EnhancedVolcano (using its default label placement).
  plt <- EnhancedVolcano(
    DEgenes_df,
    lab = DEgenes_df$gene,
    x = "logFC",
    y = "PValue",
    title = title,
    pointSize = 0.75,
    subtitle = NULL,
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    col = 'grey50',
    legendPosition = 'top',
    selectLab = select_labels
  )
  
  # Remove the default ggrepel label layer added by EnhancedVolcano.
  # (We assume it is the first layer that is a GeomTextRepel or GeomLabelRepel.)
  plt$layers <- plt$layers[!sapply(plt$layers, function(layer) {
    inherits(layer$geom, "GeomTextRepel") || inherits(layer$geom, "GeomLabelRepel")
  })]
  
 
  # Overlay an additional point layer for marker genes.
  # Using shape 21 (a filled circle with a border) so that the border is black
  # and the fill is based on cell type.
  marker_data <- subset(DEgenes_df, !is.na(cellType))
  plt <- plt +
    geom_point(data = marker_data,
               aes(x = logFC, y = -log10(PValue), fill = cellType),
               size = 4, shape = 21, color = "black", stroke = 1.5, show.legend = TRUE)
  
  # Assign one color per cell type using a palette (using Dark2).
  cell_types <- names(markers)
  if(length(cell_types) <= 8){
    ctcols <- setNames(brewer.pal(n = length(cell_types), name = "Dark2"), cell_types)
  } else {
    ctcols <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(cell_types)), cell_types)
  }
  
  # Add a manual fill scale for the marker overlay that creates a legend for cell types.
  plt <- plt +
    scale_fill_manual(
      name = "Cell Type",
      values = ctcols,
      guide = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 1.5))
    ) +
    theme(legend.position = "right")
  
  # remove the deault legend created by Enhanced Volcano plot
  plt <- plt + guides(color = "none", shape = "none", size = "none")
  
  # Add a new label layer.
  plt <- plt +
    geom_label_repel(
      data = subset(DEgenes_df, gene %in% select_labels),
      aes(x = logFC, y = -log10(PValue), label = gene),
      #nudge_x = 1,         # increase horizontal distance
      #nudge_y = 1,       # increase vertical distance 
      direction = 'both',
      box.padding = 0.75,
      point.padding = 1,
      segment.color = "black",
      force = 2,
      max.overlaps = Inf
    )
  
  print(plt)
  return(plt)
}

##########################################
### Define marker genes and make plots ###
##########################################

# vectors of canonical marker genes for each cell type
NeuN_mg <- c("Rbfox3", "Eno2", "Nefl", "Snap25", "Meg3")
PU1_mg <- c("Spi1", "Aif1", "C1qa", "Tmem119", "Apoe", "Dock8", "Irf8")
SOX10_mg <- c("Sox10", "Mbp", "Mog", "Gal3st1", "S100b", "Pdgfra", "Cnp", "Stx6")
SOX2_mg <- c("Sox2", "Sox9", "Gfap", "Atp1b2", "Slc1a3", "Slc1a2", "Apoe", "Prnp" )

# Create a named list of markers by cell type:
markers <- list(neurons = NeuN_mg, microglia = PU1_mg, oligos = SOX10_mg, astrocytes = SOX2_mg)

# generate dataframes from the qlfTest objects
NeuN_de <- data.frame(qlf.NeuN$table)
PU1_de <- data.frame(qlf.PU1$table)
SOX10_de <- data.frame(qlf.SOX10$table)
SOX2_de <- data.frame(qlf.SOX2$table)

# make the plots
NeuN_volcano <- PlotVolcanoMarkers(NeuN_de, markers, title = "NeuN+ vs. non-NeuN")
PU1_volcano <- PlotVolcanoMarkers(PU1_de, markers, title = "PU1+ vs. non-PU1")
SOX10_genes <-  PlotVolcanoMarkers(SOX10_de, markers,  title = "SOX10+ vs. non-NeuN averages")
SOX2_genes <- PlotVolcanoMarkers(SOX2_de, markers, title = "SOX2+ vs. non-SOX2 averages")

################################
#### Export DE analysis data ###
################################

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