#!/usr/bin/env Rscript

# Title: edgeR_DTE_gencode.v42.SIRVomeERCCome
# Objective: Standard edgeR pipeline on the transcript level giving differentially expressed transcripts (DTE) with gencode.v42.SIRVomeERCCome annotations and first output analyses
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(tidyverse)
  # library(RColorBrewer)
  # library(pheatmap)
  library(Glimma)
  # library(ggpmisc)
  # library(cowplot)
  # library(ggpointdensity)
  # library(rcartocolor)
  # library(ggnewscale)
  library(ggrastr)
  # library(SummarizedExperiment)
  # library(ggrepel)
  # library(ggridges)
  # library(cowplot)
  # library(limma)
  # library(scales)
  # library(extrafont)
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

# Load gencode short annotation
gtf_gencode_df_short <- read_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

# Create folders if not exists
dir.create(file.path(mydir, "edgeR", "DTE"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(mydir, "Plots", "edgeR", "DTE"), showWarnings = FALSE, recursive = TRUE)

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE) %>% 
  mutate(condition = fct_inorder(condition)) %>% 
  mutate(condition = relevel(condition, ref="control"))

# Get unique conditions
cond <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate files object used for tximport
files <- file.path(mydir, "Salmon", samples$sample)

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(files)))

# Import transcript-level counts 
catch <- catchSalmon(paths = files)

# Account for the mapping ambiguity
scaled.counts <- catch$counts/catch$annotation$Overdispersion

# Create DGEList objec
DGEList <- DGEList(counts = scaled.counts,
                   samples = samples,
                   group = samples$condition,
                   genes = catch$annotation)

# How many transcripts and how many samples?
as_tibble(dim(DGEList)) %>% 
  add_column(name = c("transcripts", "samples"), .before = "value")

# Perform filtering of lowly expressed transcripts
keep <- filterByExpr(DGEList)
# Filter out spike-in transcripts
keep_woSpike <- keep == TRUE & !str_detect(names(keep), 'ERCC') & !str_detect(names(keep), 'SIRV')

print("Sufficiently expressed transcripts")
table(keep_woSpike)

DGEList_filt <- DGEList[keep_woSpike, , keep.lib.sizes=FALSE]

# Perform normalization
DGEList_filt_norm <- normLibSizes(DGEList_filt)

# Define design matrix
design <- model.matrix(~ 0 + group,data = DGEList_filt_norm$samples)
colnames(design) <- gsub("group", "", colnames(design))

# Dispersion estimation
DGEList_filt_norm_NB <- estimateDisp(DGEList_filt_norm, design, robust=TRUE)
print("Common dispersion")
DGEList_filt_norm_NB$common.dispersion

# Plot BCV
pdf(file = file.path(mydir, "Plots", "edgeR", "DTE", "NB_dispersion_BCV_plot.pdf"))
plotBCV(DGEList_filt_norm_NB)
dev.off()


# Quasi-likelihood (QL) pipeline
fit <- glmQLFit(DGEList_filt_norm_NB, design, robust=TRUE)

# Plot QLDisp
pdf(file = file.path(mydir, "Plots", "edgeR", "DTE", "QL_dispersions_plot.pdf"))
plotQLDisp(fit)
dev.off()

# MDS plot
htmlwidgets::saveWidget(glimmaMDS(DGEList_filt_norm_NB), 
                        file.path(mydir,
                                  "edgeR",
                                  "DTE",
                                  "glimma-MDS.html"))

##
# Individual analyses ----------------------------------------------------------
##

# Build for-loop to go through all conditions vs control
mycount=2
for (i in 2:length(cond)){
  
  setwd((file.path(mydir, "edgeR")))
  
  condition <- paste(cond[mycount])
  
  # Create folders
  
  # Create "edgeR cond specific" folder for PLOTS if it does not exists
  dir.create(file.path(mydir,
                       "Plots",
                       "edgeR",
                       "DTE",
                       paste(cond[mycount])),
             showWarnings = FALSE,
             recursive = TRUE)
  
  # Create condition-specific subfolder
  dir.create(file.path(mydir,
                       "edgeR",
                       "DTE",
                       paste(cond[mycount])),
             showWarnings = FALSE,
             recursive = TRUE)
  
  # Make contrasts
  contrast <- makeContrasts(paste(cond[mycount], "-", "control"),levels=design)
  
  # Test using QL F-test
  qlf <- glmQLFTest(fit, contrast=contrast)
  
  # Obtain summary
  is.de <- decideTests(qlf, p.value=0.0001)
  summary(is.de)
  
  # Get data frame and join with annotation
  qlf_df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none")) %>% 
    rownames_to_column(var="transcript_id") %>% 
    left_join(gtf_gencode_df_short) %>% 
    relocate(gene_id, gene_name, transcript_id, transcript_name, transcript_type) %>% 
    mutate(condition_2 = condition)
  
  # Generate output csv file
  message("Generate csv output file for ", paste(cond[mycount]))
  write.csv(qlf_df, 
            file=(file.path(mydir,
                            "edgeR",
                            "DTE",
                            cond[mycount], 
                            paste0(paste(cond[mycount]),
                                   "_vs_control_edgeR_DTE_results",
                                   ".csv"))))
  
  # Generate Glimma MD plot - normal
  message("Generate Glimma MD plots for ", paste(cond[mycount]))
  
  anno <- data.frame(GeneID = qlf_df$gene_id,
                     GeneName = qlf_df$gene_name,
                     TranscriptID = qlf_df$transcript_id,
                     TranscriptName = qlf_df$transcript_name,
                     TranscriptType = qlf_df$transcript_type,
                     FDR = num(qlf_df$FDR,
                               digits = 2, 
                               notation = "sci"))
  
  htmlwidgets::saveWidget(glimmaMA(qlf,
                                   dge=DGEList_filt_norm_NB, 
                                   anno=anno,
                                   main=paste(condition, " -vs- control"),
                                   display.columns=c("GeneName", "TranscriptID", "TranscriptName", "TranscriptType", "FDR"),
                                   width = 1200,
                                   height = 1200), 
                          file.path(mydir,
                                    "edgeR",
                                    "DTE",
                                    paste(cond[mycount]),
                                    "glimma-MD.html"))
  
  # 
  # 
  # 
  # glMDPlot(qlf, 
  #          status=is.de,
  #          anno=anno,
  #          counts=DGEList_filt_norm_NB$counts,
  #          groups=samples$condition,
  #          # transform=TRUE,
  #          samples=colnames(DGEList_filt_norm_NB),
  #          main=paste(condition, " -vs- control"),
  #          display.columns=c("GeneName", "TranscriptID", "TranscriptName", "TranscriptType"),
  #          sample.cols=
  #          path=(file.path(mydir,
  #                          "edgeR",
  #                          "DTE_NMDRegHumanTxome",
  #                          paste(cond[mycount]))),
  #          folder=paste0("glimma-MD"), launch=FALSE)
  
  # Prepare for plotting
  # qlf_df_plot <- qlf_df %>% 
  #   mutate(FDR = replace(FDR, FDR == 0, 1e-320)) %>% 
  #   mutate(logFC = replace(logFC, logFC > 10, 10)) %>% 
  #   mutate(logFC = replace(logFC, logFC < -10, -10)) %>% 
  #   dplyr::filter(!is.na(FDR))
  # 
  # qlf_df_plot_ECDF <- qlf_df_plot %>% 
  #   mutate(type = case_when(!transcript_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other", 
  #                           transcript_biotype == "protein_coding" ~ "coding",
  #                           transcript_biotype == "nonsense_mediated_decay" ~ "NMD",
  #                           transcript_biotype == "lncRNA" ~ "lncRNA"
  #                           ))
  # 
  # qlf_df_plot_sig <- subset(qlf_df_plot, FDR < 0.0001 & abs(logFC) > 1) %>% 
  #   mutate(type = case_when(!transcript_biotype %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other", 
  #                           transcript_biotype == "protein_coding" ~ "coding",
  #                           transcript_biotype == "nonsense_mediated_decay" ~ "NMD",
  #                           transcript_biotype == "lncRNA" ~ "lncRNA"))%>% 
  #   mutate(UpDown = case_when(logFC > 1 ~ "up",
  #                             logFC < -1 ~ "down"))
  # 
  # # Create unfiltered volcano plot
  # myplot <- ggplot(data = qlf_df_plot, aes(x=logFC, y=-log10(FDR))) + 
  #   theme_VB() +
  #   coord_cartesian(xlim = c(-10, 10), ylim = c(0, 330), clip = 'on') +
  #   ggrastr::rasterise(geom_point(data=dplyr::filter(qlf_df_plot, !transcript_biotype %in% c("protein_coding", "nonsense_mediated_decay")), 
  #                                 alpha=0.5,
  #                                 size=0.75,
  #                                 aes(color="darkgray")), dpi=600) + 
  #   ggrastr::rasterise(geom_point(data=dplyr::filter(qlf_df_plot, transcript_biotype == "protein_coding"), 
  #                                 alpha=0.5,
  #                                 size=0.75,
  #                                 aes(color="#008080")), dpi=600) + 
  #   ggrastr::rasterise(geom_point(data=dplyr::filter(qlf_df_plot, transcript_biotype == "nonsense_mediated_decay"), 
  #                                 alpha=0.5,
  #                                 size=0.75,
  #                                 aes(color="#b2182b")), dpi=600) + 
  #   geom_hline(yintercept=-log10(0.0001), linetype="dashed", linewidth=0.2) +
  #   geom_hline(aes(yintercept=320), linetype="dashed", color="gray", linewidth=0.2) +
  #   geom_vline(xintercept = 1, linetype="dashed", linewidth=0.2) +
  #   geom_vline(xintercept = -1, linetype="dashed", linewidth=0.2) +	
  #   labs(title = paste0("Control_vs_",
  #                       paste0(cond[mycount]),
  #                       "\n"), 
  #        subtitle=paste0("edgeR (", packageVersion("edgeR"), ") ",
  #                        "DTE analysis",
  #                        "\n"),
  #        x = "log2 FoldChange",
  #        y = "-log10 FDR") +
  #   scale_color_identity(name = "Transcript type",
  #                        breaks = c("darkgray", "#008080", "#b2182b"),
  #                        labels = c("other", "coding", "NMD"),
  #                        guide = "legend") 
  # 
  # myplot <- myplot + annotate("text", x=-8.5, y=335, label="Max p-value", size = 6/.pt, color = "gray")
  # 
  # resDF_final_plot_summary <- resDF_final_plot_sig %>% 
  #   group_by(type, UpDown) %>% 
  #   summarize(n = n()) %>% 
  #   pivot_wider(names_from = UpDown, values_from = n) %>% 
  #   arrange(factor(type, levels = c('coding', 'NMD', 'other')))
  # 
  # # Create filtered density plot
  # myplot2 <- ggplot(data = resDF_final_plot_sig, aes(log2FoldChange)) + 
  #   theme_VB() + 	
  #   theme(legend.position="none") +
  #   xlim(-10, 10) +
  #   geom_density(data=dplyr::filter(resDF_final_plot_sig, !transcript_type %in% c("protein_coding", "nonsense_mediated_decay")), color="darkgray", fill="darkgray", alpha = 0.2) +
  #   geom_density(data=dplyr::filter(resDF_final_plot_sig, transcript_type == "protein_coding"), color="#008080", fill="#008080", alpha = 0.2) +
  #   geom_density(data=dplyr::filter(resDF_final_plot_sig, transcript_type == "nonsense_mediated_decay"), color="#b2182b", fill="#b2182b", alpha = 0.2) +
  #   annotate(geom = "table", x = -Inf, y = Inf, label = (resDF_final_plot_summary), 
  #            table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
  #                                      core=list(fg_params=list(col = c("#008080", "#b2182b", "darkgray"),fontface=1)),
  #                                      colhead=list(fg_params=list( fontface=1))),
  #            vjust = 1, hjust = -0.1) +
  #   labs(x = "log2 FoldChange", y = "filtered density", caption = paste0("\nFilter: |log2FoldChange| > 1 & p.adjust < 0.0001")) 
  # 
  # finalplot <- plot_grid(myplot, myplot2, ncol = 1, align = "v", rel_heights = c(2, 1))
  # 
  # # Save plot
  # ggsave(file.path(mydir,
  #                  "Plots",
  #                  "DESeq2",
  #                  "DTE_BatchCorr",
  #                  FiltType,
  #                  paste0(cond[mycount]),
  #                  paste0(paste(cond[mycount]),
  #                         "_Volcano_Density_",
  #                         Lfc_method,
  #                         ".pdf")),
  #        finalplot,
  #        width = 7,
  #        height = 8,
  #        units = "cm",
  #        device=cairo_pdf,
  #        bg = "transparent")
  
  if (mycount == 2) {
    DTE_cat <- qlf_df
  } else {
    DTE_cat <- bind_rows(DTE_cat,
                         qlf_df)
  }
  
  # Counter increment
  mycount = mycount + 1
  
}

# Combined analysis  -------------------------------------------------------

write.csv(DTE_cat, file=file.path(mydir, "edgeR", "DTE", paste0("edgeR_DTE_cat",
                                                                           ".csv")))
