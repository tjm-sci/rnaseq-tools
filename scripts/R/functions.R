# Author: T. Murphy
# Date: 2025-03-04

# Functions for processing and visualisation of FANS RNA-seq data
# These functions are sourced in the analysis script "fans_rnaseq.R"

#########################################
# plotting functions for edgeR pipeline #
#########################################


# Bar charts for library size
#############################
lib_size_plt <- function(lib_sizes, y_range, title){
  # barplot with sensible aesthetics
  par(mar = c(8, 6, 4, 2), mgp = c(5, 1, 0)) # set margins: bottom, left, top, right
  
  barplot(lib_sizes,
          main = title,
          las = 2,
          xlab = NULL,
          ylim = y_range,
          ylab = "Library Size")
  
}

# pheatmap of sample:sample correlations
########################################
sample_heatmap <- function(y){
  # heatmap of correlation coefficients
  corr <- cor(y$counts)
  pheatmap(corr)
}


# MDS plot function
###################
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

# Enhanced Volcano plot
#######################
PlotVolcano <- function(DEgenes_df, selectLab, FCcutoff = 2, title){
  
  # If there is no "gene" column, create it from rownames
  if(!"gene" %in% colnames(DEgenes_df)){
    DEgenes_df$gene <- rownames(DEgenes_df)
  }
  
  # Create the base volcano plot using EnhancedVolcano
  plt <- EnhancedVolcano(
    DEgenes_df,
    lab = DEgenes_df$gene,
    selectLab = selectLab,
    x = "logFC",
    y = "PValue",
    FCcutoff = FCcutoff,
    title = title,
    pointSize = 0.75,
    subtitle = NULL,
    boxedLabels = FALSE,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    legendPosition = "right"
  )
  
  print(plt)
  return(plt)
}

# Enhanced Volcano with colour coded, manually-specified marker genes labelled.
##############################################################################
MarkerGeneVolcano <- function(DEgenes_df, markers, FCcutoff=2, title = ""){
  
  # DEgenes_df: data frame of DE results (e.g., from edgeR)
  # Must contain columns "logFC" and "PValue" and have rownames as gene symbols.
  # markers: a named list of marker gene vectors, one element per cell type,
  #          e.g. list(neurons = c("Rbfox3", "Eno2", "Nefl", "Snap25")...)
  # title: Plot title
  # FCcutoff: Fold-change cuttoff drawn as a vertical line.
  
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
    FCcutoff = FCcutoff,
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