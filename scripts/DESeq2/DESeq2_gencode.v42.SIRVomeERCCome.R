#!/usr/bin/env Rscript

# Title: DESeq2_gencode.v42.SIRVomeERCCome
# Objective: Standard DESeq2 pipeline giving differentially expressed genes (DGE) with gencode annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(tidyverse)
  library(DESeq2)
  library(RColorBrewer)
  library(pheatmap)
  library(Glimma)
  library(ggpmisc)
  library(cowplot)
  library(ggpointdensity)
  library(rcartocolor)
  library(ggnewscale)
  library(ggrastr)
  library(SummarizedExperiment)
  library(ggrepel)
  library(ggridges)
  library(cowplot)
  library(limma)
  library(scales)
  library(extrafont)
})

# Define Plot Theme ------------------------------------------------------------

theme_VB <- function(){ 
  
  theme_classic() %+replace%    #replace elements we want to change
    
    theme(legend.position="top", 
          strip.background = element_blank(), 
          plot.background = element_blank(), 
          panel.background = element_blank(), 
          legend.background = element_rect(linewidth=0.2,
                                           linetype="solid",
                                           color ="darkgray"),
          axis.line=element_line(linewidth=0.1),
          axis.ticks=element_line(linewidth=0.1),
          axis.text=element_text(size=6), 
          axis.title=element_text(size=6), 
          legend.title = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.justification="left",
          legend.margin=margin(0,5,0,0),
          legend.box.margin=margin(1,1,-5,1),
          legend.spacing.x = unit(-0.125, "cm"),
          # legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(size = 6, hjust = 0, face = "bold"), 
          plot.subtitle = element_text(size = 5, hjust = 0), 
          plot.caption = element_text(size = 5, hjust = 0),
          text=element_text(family="Arial"),
          # plot.title.position = "plot",
          # plot.caption.position =  "plot"
    )
}

##
# Data preparation ----------------------------------------------------------
##

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 1)

# mydir is the first argument = working directory
mydir = arguments$args[1]

# Survey name
myname <- tail(unlist(strsplit(mydir,"/")), n=1)

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

# Create folders if not exists
dir.create(file.path(mydir, "DESeq2", "DGE"), showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(mydir, "DESeq2", "DTE"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(mydir, "Plots", "DESeq2", "DGE"), showWarnings = FALSE, recursive = TRUE)
# dir.create(file.path(mydir, "Plots", "DESeq2", "DTE"), showWarnings = FALSE, recursive = TRUE)

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample, "quant.sf")

# Supplement with sample IDs as "names"
names(files) <- samples$sample

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

## Gencode annotation loading ----------------------------------------------

# Prepared only once:

# gtf_gencode_df_alt <- rtracklayer::readGFF("/home/volker/reference/Gencode/gencode.v42.annotation.gff3.gz",
#                                            filter = list(type=c("gene",
#                                                                 "transcript")))
# gtf_gencode_tbl <- as_tibble(gtf_gencode_df_alt)
# 
# gtf_gencode_tbl_tags <- gtf_gencode_tbl %>% 
#   mutate(ensembl_canonical = case_when(str_detect(tag, "Ensembl_canonical") ~ TRUE,
#                                                   TRUE ~ FALSE)) %>% 
#   mutate(ensembl_basic = case_when(str_detect(tag, "basic") ~ TRUE,
#                                        TRUE ~ FALSE)) %>% 
#   mutate(appris_principal = case_when(str_detect(tag, "appris_principal") ~ TRUE,
#                                    TRUE ~ FALSE))
# 
# gtf_gencode_df_short <- gtf_gencode_tbl_tags %>% 
#   dplyr::select("type",
#                 "gene_id",
#                 "gene_name",
#                 "gene_type", 
#                 "transcript_id",
#                 "transcript_type",
#                 "transcript_name",
#                 "transcript_support_level",
#                 "ensembl_canonical",
#                 "ensembl_basic",
#                 "appris_principal") %>%
#   distinct()
# 
# gtf_gencode_df_short %>% write_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

gtf_gencode_df_short <- read_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

# Import tx2gene file which references each transcript to the corresponding gene ID
tx2gene <- read_tsv(file.path(ref_dir, "tx2gene.gencode.v42.SIRV.ERCC.tsv"))
txi <- tximport(files, 
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion = FALSE)

# Generate DESeqDataSet using samples and condition as parameters
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

# Set control condition as reference
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "control")

# Pre-filter more stringently: require at least half of the samples to have equal or more than 10 counts per gene
# Also remove ERCC and SIRV spike-in genes
keep_string <- rowSums(counts(ddsTxi) >= 10 ) >= nrow(samples)/2 & !str_detect(rownames(counts(ddsTxi)), 'ERCC') & !str_detect(rownames(counts(ddsTxi)), 'SIRV')

ddsTxi_string <- ddsTxi[keep_string,]

# Perform the DESeq analysis - on stringently pre-filtered ddsTxi
dds_string <- DESeq(ddsTxi_string)

##
# Sample-to-sample distance ----------------------------------------------------------
##

## Stringent -----------------------------------------------------------------

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Generate heatmap of the sample-to-sample distances
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, filename=(file.path(mydir, "Plots", "DESeq2", "DGE", "DGE_sampleDistanceHeatmap_stringent.pdf")))

# Principal component plot of the samples
plt <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(plt, "percentVar"))

plt_plot <- ggplot(plt, aes(PC1, PC2, color=condition)) +
  theme_classic() +	
  geom_point(size=3) +		
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.title = element_text(size = 10), legend.text = element_text(size = 12), plot.title = element_text(size = 14), plot.subtitle = element_text(size = 12), plot.caption = element_text(size = 10)) +
  labs(title = paste(myname), subtitle = "DESeq2 PCA analysis", caption = paste0("DESeq2 (", packageVersion("DESeq2"), ")"), x = paste0("PC1: ",percentVar[1],"% variance"), y = paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggsave(file.path(mydir, "Plots", "DESeq2", "DGE", "DGE_PCA_plot_stringent.pdf"), plt_plot, width = 6, height=4, device=cairo_pdf)

##
# Individual analyses ----------------------------------------------------------
##


# Build for-loop to go through all conditions vs control
mycount=2
for (i in 2:length(cond)){

  # String -----------------------------------------------------------------
  
  dds = dds_string
  FiltType = "String"
  
  # Perform DESeq2
  
  setwd((file.path(mydir, "DESeq2")))
  
  condition <- paste(cond[mycount])
  
  # Create folders

  # Create "DESeq2 cond specific" folder for PLOTS if it does not exists
  dir.create(file.path(mydir,
                       "Plots",
                       "DESeq2",
                       "DGE",
                       FiltType,
                       paste(cond[mycount])),
             showWarnings = FALSE,
             recursive = TRUE)
  
  # Create condition-specific subfolder
  dir.create(file.path(mydir,
                       "DESeq2",
                       "DGE",
                       FiltType,
                       paste(cond[mycount])),
             showWarnings = FALSE,
             recursive = TRUE)
  
  ## Normal ------------------------------------------------------------------
  
  Lfc_method = "Normal"
  
  res <- results(dds,
                 alpha=0.0001,
                 contrast=c("condition",paste(cond[mycount]),"control"))	
  
  # Convert results file to tibble
  resDF <- as_tibble(res,
                     rownames = "gene_id") %>% 
    mutate(condition_2 = condition)
  
  # Get additional IDs from Gencode annotation
  resDF_final <- left_join(resDF,
                           gtf_gencode_df_short %>% filter(type=="gene") %>% select(gene_id, gene_name, gene_type),
                           by = c("gene_id" = "gene_id"))
  
  # Generate output csv file
  message("Generate csv output file for ", paste(cond[mycount]))
  write.csv(resDF_final, 
            file=(file.path(mydir,
                            "DESeq2",
                            "DGE",
                            FiltType,
                            cond[mycount], 
                            paste0(paste(cond[mycount]),
                                   "_vs_control_DESeq2_results_",
                                   Lfc_method,
                                   ".csv"))))
  
  # Generate Glimma MD plot - normal
  message("Generate Glimma MD plots for ", paste(cond[mycount]))
  status <- as.numeric(res$padj < .0001)
  anno <- data.frame(GeneID = resDF_final$gene_id,
                     GeneName = resDF_final$gene_name,
                     GeneType = resDF_final$gene_type)
  glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE),
           groups=dds$condition, transform=TRUE,
           samples=colnames(dds), anno=anno,
           path=(file.path(mydir,
                           "DESeq2",
                           "DGE",
                           FiltType,
                           paste(cond[mycount]))),
           folder=paste0("glimma-MD_", Lfc_method), launch=FALSE)
  
  # Prepare for plotting
  resDF_final_plot <- resDF_final %>% 
    mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange > 10, 10)) %>% 
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange < -10, -10)) %>% 
    dplyr::filter(!is.na(padj))
  
  resDF_final_plot_ECDF <- resDF_final_plot %>% 
    mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                            gene_type == "protein_coding" ~ "coding",
                            gene_type == "lncRNA" ~ "lncRNA"))
  
  resDF_final_plot_sig <- subset(resDF_final_plot, padj < 0.0001 & abs(log2FoldChange) > 1) %>% 
    mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                            gene_type == "protein_coding" ~ "coding",
                            gene_type == "lncRNA" ~ "lncRNA")) %>% 
    mutate(UpDown = case_when(log2FoldChange > 1 ~ "up",
                              log2FoldChange < -1 ~ "down"))
  
  # Create unfiltered volcano plot
  myplot <- ggplot(data = resDF_final_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
    theme_VB() +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 330), clip = 'on') +
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, !gene_type %in% c("protein_coding", "lncRNA")), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="darkgray")), dpi=600) + 
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, gene_type == "protein_coding"), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="#008080")), dpi=600) + 
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, gene_type == "lncRNA"), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="#CA562C")), dpi=600) + 
    geom_hline(yintercept=-log10(0.0001), linetype="dashed", linewidth=0.2) +
    geom_hline(aes(yintercept=320), linetype="dashed", color="gray", linewidth=0.2) +
    geom_vline(xintercept = 1, linetype="dashed", linewidth=0.2) +
    geom_vline(xintercept = -1, linetype="dashed", linewidth=0.2) +	
    labs(title = paste0("Control_vs_",
                        paste0(cond[mycount]),
                        "\n"), 
         subtitle=paste0("DESeq2 (", packageVersion("DESeq2"), ") ",
                         "DGE analysis",
                         " | pre-filter: ",
                         FiltType,
                         " | lfc: ",
                         Lfc_method,
                         "\n"),
         x = "log2 FoldChange",
         y = "-log10 p.adjust") +
    scale_color_identity(name = "Gene type",
                         breaks = c("darkgray", "#008080", "#CA562C"),
                         labels = c("other", "coding", "lncRNA"),
                         guide = "legend") 
  
  myplot <- myplot + annotate("text", x=-8.5, y=335, label="Max p-value", size = 6/.pt, color = "gray")
  
  resDF_final_plot_summary <- resDF_final_plot_sig %>% 
    group_by(type, UpDown) %>% 
    summarize(n = n()) %>% 
    pivot_wider(names_from = UpDown, values_from = n) %>% 
    arrange(factor(type, levels = c('coding', 'lncRNA', 'other')))
  
  # Create filtered density plot
  myplot2 <- ggplot(data = resDF_final_plot_sig, aes(log2FoldChange)) + 
    theme_VB() + 	
    theme(legend.position="none") +
    xlim(-10, 10) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, !gene_type %in% c("protein_coding", "lncRNA")), color="darkgray", fill="darkgray", alpha = 0.2) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, gene_type == "protein_coding"), color="#008080", fill="#008080", alpha = 0.2) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, gene_type == "lncRNA"), color="#CA562C", fill="#CA562C", alpha = 0.2) +
    annotate(geom = "table", x = -Inf, y = Inf, label = (resDF_final_plot_summary), 
             table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                       core=list(fg_params=list(col = c("#008080", "#CA562C", "darkgray"),fontface=1)),
                                       colhead=list(fg_params=list( fontface=1))),
             vjust = 1, hjust = -0.1) +
    labs(x = "log2 FoldChange", y = "filtered density", caption = paste0("\nFilter: |log2FoldChange| > 1 & p.adjust < 0.0001")) 
  
  finalplot <- plot_grid(myplot, myplot2, ncol = 1, align = "v", rel_heights = c(2, 1))
  
  # Save plot
  ggsave(file.path(mydir,
                   "Plots",
                   "DESeq2",
                   "DGE",
                   FiltType,
                   paste0(cond[mycount]),
                   paste0(paste(cond[mycount]),
                          "_Volcano_Density_",
                          Lfc_method,
                          ".pdf")),
         finalplot,
         width = 7,
         height = 8,
         units = "cm",
         device=cairo_pdf,
         bg = "transparent")
  
  ###
  # ECDF plots
  ###
  
  # Generate named vector for color coding
  my_colors <- c("other" = "darkgray", "coding" = "#008080", "lncRNA" = "#CA562C")
  
  # Summarize n for type
  resDF_final_plot_all_summary <- resDF_final_plot_ECDF %>% 
    group_by(type) %>% 
    summarize(n = n()) %>% 
    arrange(factor(type, levels = c('coding', 'lncRNA', 'other')))
  
  # Kolmogorov-Smirnov Test
  ks_resDF_coding_lncRNA <- ks.test(resDF_final_plot_ECDF %>% 
                                   dplyr::filter(type == "coding") %>% 
                                   dplyr::pull(log2FoldChange), 
                                 resDF_final_plot_ECDF %>% 
                                   dplyr::filter(type == "lncRNA") %>% 
                                   dplyr::pull(log2FoldChange))
  
  ks_resDF_coding_lncRNA$p.value <- ifelse(ks_resDF_coding_lncRNA$p.value == 0, "< 2.2e-16", signif(ks_resDF_coding_lncRNA$p.value,3))
  
  # Plot
  ECDF_plot <- ggplot(resDF_final_plot_ECDF, aes(x=log2FoldChange, col=type)) + 
    theme_VB() + 	
    theme(legend.spacing.x = unit(0.1, "cm")) +
    geom_vline(xintercept = 0, linetype="dashed") +	
    stat_ecdf(geom = "step") +
    xlim(-5, 5) +
    annotate(geom = "table", x = -5, y = 1, label = (resDF_final_plot_all_summary), 
             table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                       core=list(fg_params=list(col = c("#008080", "#CA562C", "darkgray"),fontface=1)),
                                       colhead=list(fg_params=list( fontface=1))),
             vjust = 1, hjust = -0.1) +
    annotate("text", x = Inf, y = -Inf, size = 6/.pt, label = paste0("KS test \n (coding vs lncRNA) \n p-value: ",
                                                                     paste(ks_resDF_coding_lncRNA$p.value),
                                                                     "\n D: ",
                                                                     paste(signif(ks_resDF_coding_lncRNA$statistic,3))),
             vjust = -0.25, hjust = 1) +
    scale_color_manual(values=my_colors) +
    labs(title = paste0("Control_vs_",
                        paste0(cond[mycount]),
                        "\n"),
         subtitle=paste0("DESeq2 (", packageVersion("DESeq2"), ") ",
                         "DGE analysis",
                         " | pre-filter: ",
                         FiltType,
                         " | lfc: ",
                         Lfc_method,
                         "\n"),
         x = "log2 FoldChange",
         y = "ECDF",
         color = "Gene type") 
  
  
  
  # Save plot
  ggsave(file.path(mydir,
                   "Plots",
                   "DESeq2",
                   "DGE",
                   FiltType,
                   paste0(cond[mycount]),
                   paste0(paste(cond[mycount]),
                          "_ECDF_",
                          Lfc_method,
                          ".pdf")),
         ECDF_plot,
         width = 7,
         height = 6,
         units = "cm",
         device=cairo_pdf, bg = "transparent")
  
  ## lfcShrunk ---------------------------------------------------------------
  
  Lfc_method = "Shrunk"
  
  res <- lfcShrink(dds,
                   parallel = TRUE,
                   coef=paste0("condition_",paste(cond[mycount]),"_vs_control"))	
  
  # Convert results file to tibble
  resDF <- as_tibble(res,
                     rownames = "gene_id") %>% 
    mutate(condition_2 = condition)
  
  # Get additional IDs from Gencode annotation
  resDF_final <- left_join(resDF,
                           gtf_gencode_df_short %>% filter(type=="gene") %>% select(gene_id, gene_name, gene_type),
                           by = c("gene_id" = "gene_id"))
  
  # Generate output csv file
  message("Generate csv output file for ", paste(cond[mycount]))
  write.csv(resDF_final, 
            file=(file.path(mydir,
                            "DESeq2",
                            "DGE",
                            FiltType,
                            cond[mycount], 
                            paste0(paste(cond[mycount]),
                                   "_vs_control_DESeq2_results_",
                                   Lfc_method,
                                   ".csv"))))
  
  # Generate Glimma MD plot - normal
  message("Generate Glimma MD plots for ", paste(cond[mycount]))
  status <- as.numeric(res$padj < .0001)
  anno <- data.frame(GeneID = resDF_final$gene_id,
                     GeneName = resDF_final$gene_name,
                     GeneType = resDF_final$gene_type)
  glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE),
           groups=dds$condition, transform=TRUE,
           samples=colnames(dds), anno=anno,
           path=(file.path(mydir,
                           "DESeq2",
                           "DGE",
                           FiltType,
                           paste(cond[mycount]))),
           folder=paste0("glimma-MD_", Lfc_method), launch=FALSE)
  
  # Prepare for plotting
  resDF_final_plot <- resDF_final %>% 
    mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange > 10, 10)) %>% 
    mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange < -10, -10)) %>% 
    dplyr::filter(!is.na(padj))
  
  resDF_final_plot_ECDF <- resDF_final_plot %>% 
    mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                            gene_type == "protein_coding" ~ "coding",
                            gene_type == "lncRNA" ~ "lncRNA"))
  
  resDF_final_plot_sig <- subset(resDF_final_plot, padj < 0.0001 & abs(log2FoldChange) > 1) %>% 
    mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                            gene_type == "protein_coding" ~ "coding",
                            gene_type == "lncRNA" ~ "lncRNA")) %>% 
    mutate(UpDown = case_when(log2FoldChange > 1 ~ "up",
                              log2FoldChange < -1 ~ "down"))
  
  # Create unfiltered volcano plot
  myplot <- ggplot(data = resDF_final_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
    theme_VB() +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 330), clip = 'on') +
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, !gene_type %in% c("protein_coding", "lncRNA")), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="darkgray")), dpi=600) + 
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, gene_type == "protein_coding"), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="#008080")), dpi=600) + 
    ggrastr::rasterise(geom_point(data=dplyr::filter(resDF_final_plot, gene_type == "lncRNA"), 
                                  alpha=0.5,
                                  size=0.75,
                                  aes(color="#CA562C")), dpi=600) + 
    geom_hline(yintercept=-log10(0.0001), linetype="dashed", linewidth=0.2) +
    geom_hline(aes(yintercept=320), linetype="dashed", color="gray", linewidth=0.2) +
    geom_vline(xintercept = 1, linetype="dashed", linewidth=0.2) +
    geom_vline(xintercept = -1, linetype="dashed", linewidth=0.2) +	
    labs(title = paste0("Control_vs_",
                        paste0(cond[mycount]),
                        "\n"), 
         subtitle=paste0("DESeq2 (", packageVersion("DESeq2"), ") ",
                         "DGE analysis",
                         " | pre-filter: ",
                         FiltType,
                         " | lfc: ",
                         Lfc_method,
                         "\n"),
         x = "log2 FoldChange",
         y = "-log10 p.adjust") +
    scale_color_identity(name = "Gene type",
                         breaks = c("darkgray", "#008080", "#CA562C"),
                         labels = c("other", "coding", "lncRNA"),
                         guide = "legend") 
  
  myplot <- myplot + annotate("text", x=-8.5, y=335, label="Max p-value", size = 6/.pt, color = "gray")
  
  resDF_final_plot_summary <- resDF_final_plot_sig %>% 
    group_by(type, UpDown) %>% 
    summarize(n = n()) %>% 
    pivot_wider(names_from = UpDown, values_from = n) %>% 
    arrange(factor(type, levels = c('coding', 'lncRNA', 'other')))
  
  # Create filtered density plot
  myplot2 <- ggplot(data = resDF_final_plot_sig, aes(log2FoldChange)) + 
    theme_VB() + 	
    theme(legend.position="none") +
    xlim(-10, 10) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, !gene_type %in% c("protein_coding", "lncRNA")), color="darkgray", fill="darkgray", alpha = 0.2) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, gene_type == "protein_coding"), color="#008080", fill="#008080", alpha = 0.2) +
    geom_density(data=dplyr::filter(resDF_final_plot_sig, gene_type == "lncRNA"), color="#CA562C", fill="#CA562C", alpha = 0.2) +
    annotate(geom = "table", x = -Inf, y = Inf, label = (resDF_final_plot_summary), 
             table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                       core=list(fg_params=list(col = c("#008080", "#CA562C", "darkgray"),fontface=1)),
                                       colhead=list(fg_params=list( fontface=1))),
             vjust = 1, hjust = -0.1) +
    labs(x = "log2 FoldChange", y = "filtered density", caption = paste0("\nFilter: |log2FoldChange| > 1 & p.adjust < 0.0001")) 
  
  finalplot <- plot_grid(myplot, myplot2, ncol = 1, align = "v", rel_heights = c(2, 1))
  
  # Save plot
  ggsave(file.path(mydir,
                   "Plots",
                   "DESeq2",
                   "DGE",
                   FiltType,
                   paste0(cond[mycount]),
                   paste0(paste(cond[mycount]),
                          "_Volcano_Density_",
                          Lfc_method,
                          ".pdf")),
         finalplot,
         width = 7,
         height = 8,
         units = "cm",
         device=cairo_pdf,
         bg = "transparent")
  
  ###
  # ECDF plots
  ###
  
  # Generate named vector for color coding
  my_colors <- c("other" = "darkgray", "coding" = "#008080", "lncRNA" = "#CA562C")
  
  # Summarize n for type
  resDF_final_plot_all_summary <- resDF_final_plot_ECDF %>% 
    group_by(type) %>% 
    summarize(n = n()) %>% 
    arrange(factor(type, levels = c('coding', 'lncRNA', 'other')))
  
  # Kolmogorov-Smirnov Test
  ks_resDF_coding_lncRNA <- ks.test(resDF_final_plot_ECDF %>% 
                                   dplyr::filter(type == "coding") %>% 
                                   dplyr::pull(log2FoldChange), 
                                 resDF_final_plot_ECDF %>% 
                                   dplyr::filter(type == "lncRNA") %>% 
                                   dplyr::pull(log2FoldChange))
  
  ks_resDF_coding_lncRNA$p.value <- ifelse(ks_resDF_coding_lncRNA$p.value == 0, "< 2.2e-16", signif(ks_resDF_coding_lncRNA$p.value,3))
  
  # Plot
  ECDF_plot <- ggplot(resDF_final_plot_ECDF, aes(x=log2FoldChange, col=type)) + 
    theme_VB() + 	
    theme(legend.spacing.x = unit(0.1, "cm")) +
    geom_vline(xintercept = 0, linetype="dashed") +	
    stat_ecdf(geom = "step") +
    xlim(-5, 5) +
    annotate(geom = "table", x = -5, y = 1, label = (resDF_final_plot_all_summary), 
             table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                       core=list(fg_params=list(col = c("#008080", "#CA562C", "darkgray"),fontface=1)),
                                       colhead=list(fg_params=list( fontface=1))),
             vjust = 1, hjust = -0.1) +
    annotate("text", x = Inf, y = -Inf, size = 6/.pt, label = paste0("KS test \n (coding vs lncRNA) \n p-value: ",
                                                                     paste(ks_resDF_coding_lncRNA$p.value),
                                                                     "\n D: ",
                                                                     paste(signif(ks_resDF_coding_lncRNA$statistic,3))),
             vjust = -0.25, hjust = 1) +
    scale_color_manual(values=my_colors) +
    labs(title = paste0("Control_vs_",
                        paste0(cond[mycount]),
                        "\n"),
         subtitle=paste0("DESeq2 (", packageVersion("DESeq2"), ") ",
                         "DGE analysis",
                         " | pre-filter: ",
                         FiltType,
                         " | lfc: ",
                         Lfc_method,
                         "\n"),
         x = "log2 FoldChange",
         y = "ECDF",
         color = "Gene type") 
  
  
  
  # Save plot
  ggsave(file.path(mydir,
                   "Plots",
                   "DESeq2",
                   "DGE",
                   FiltType,
                   paste0(cond[mycount]),
                   paste0(paste(cond[mycount]),
                          "_ECDF_",
                          Lfc_method,
                          ".pdf")),
         ECDF_plot,
         width = 7,
         height = 6,
         units = "cm",
         device=cairo_pdf, bg = "transparent")
  
  if (mycount == 2) {
    DGE_cat <- resDF_final
  } else {
    DGE_cat <- bind_rows(DGE_cat,
                         resDF_final)
  }
  
  # Counter increment
  mycount = mycount + 1
  
}

# Combined analysis  -------------------------------------------------------

write.csv(DGE_cat, file=file.path(mydir, "DESeq2", "DGE", paste0("DESeq2_DGE_cat",
                                                                ".csv")))

## DGE Ridgeplot -------------------------------------------------------

# Ridgeplot
DGE_ridgeplot <- DGE_cat %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  ggplot(aes(x=log2FoldChange,
             y=condition_2
  )) + 
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_rect(linewidth=0.1, 
                                         linetype="solid", 
                                         colour ="black"), 
        panel.grid.major.y = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_rect(aes(xmin=-1, xmax=1, ymin=-Inf, ymax=Inf),
            alpha=0.1,
            color=NA,
            fill="#f0f0f0") +
  geom_vline(xintercept = -1, linetype="dashed", color="darkgray", linewidth = 0.5) +
  geom_vline(xintercept = 1, linetype="dashed", color="darkgray", linewidth = 0.5) +
  stat_density_ridges(aes(height = after_stat(ndensity),
                          fill = factor(after_stat(quantile))),
                      scale = 0.9,
                      quantiles = c(0.25,0.5,0.75),
                      calc_ecdf = TRUE,
                      geom = "density_ridges_gradient",
                      rel_min_height = 0.01) +
  geom_vline(xintercept = 0, linetype="solid", color="black", linewidth = 0.5) +
  scale_fill_manual(values = c("#bfd3e6", "#8c96c6", "#88419d", "#4d004b"),
                    name = "Quantile") +
  scale_x_continuous(breaks=c(-7.5, -5, -2.5, -1, 0, 1, 2.5, 5, 7.5),
                     limits=c(-7.5, 7.5)) +
  labs(x = "log2 FoldChange",
       y = "Condition") 

## DGE NMD Counts ----------------------------------------------------------

DGE_cat_counts_absolute <- DGE_cat %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  mutate(UpDown = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up",
                            log2FoldChange < -1 & padj < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  dplyr::count(condition_2, UpDown) %>% 
  mutate(UpDown = fct_rev(fct_inorder(as_factor(UpDown)))) %>% 
  group_by(condition_2) %>% 
  mutate(n_per = round(n / sum(n), 2)) %>% 
  ungroup() %>% 
  # mutate(condition_2 = fct_expand(condition_2, "Placeholder", after = Inf)) %>% 
  ggplot(aes(x=n,
             y=condition_2, 
             fill = UpDown)) + 
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_rect(linewidth=0.1, 
                                         linetype="solid", 
                                         colour ="black"), 
        panel.grid.major.y = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_bar(stat="identity",
           color="black",
           linewidth = 0.2) +
  geom_text(aes(label = case_when(n_per > 0.2 & UpDown != "ns" ~ paste0(n,
                                                                        "\n",
                                                                        n_per*100,
                                                                        "%"),
                                  TRUE ~ "")),
            color="white",
            size = 6*0.36,
            position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = case_when(n_per > 0.2 & UpDown == "ns" ~ paste0(n,
                                                                        "\n",
                                                                        n_per*100,
                                                                        "%"),
                                  TRUE ~ "")),
            color="black",
            size = 6*0.36,
            position = position_stack(vjust = 0.5)) +
  scale_y_discrete(drop=F) +
  scale_fill_manual(values = c("up" = "#4d004b",
                               "ns" = "gray",
                               "down" = "#0570b0")) +
  labs(x = "Number of identified genes (n)",
       y = "",
       color="Significantly Up-/Down-regulated",
       caption=paste0("\nDESeq2-DGE - Filter: |log2FoldChange| > 1 & padj < 0.0001"))


DGE_ridgeplot_counts <- plot_grid(DGE_ridgeplot + 
                                    labs(title="DESeq2 DGE - gencode.v42.SIRVomeERCCome") +
                                    guides(fill = guide_legend(title.position = "top",
                                                               title.hjust = 0.5,
                                                               label.position = "bottom",
                                                               label.hjust = 0.5)),
                                  DGE_cat_counts_absolute + 
                                    theme(axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank()) +
                                    guides(fill = guide_legend(title = "Sig. regulated",
                                                               title.position = "top",
                                                               title.hjust = 0.5,
                                                               label.position = "bottom",
                                                               label.hjust = 0.5,
                                                               reverse = TRUE)),
                                  nrow=1,
                                  align = "h", 
                                  axis = "bt",
                                  rel_widths = c(2,1))

ggsave(file.path(mydir, "Plots", "DESeq2", "DGE", "DESeq2_DGE_summary.pdf"),
       width = 12,
       height=12,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

writeLines(capture.output(sessionInfo()), paste0(mydir, "/DESeq2/DGE/DESeq2_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
