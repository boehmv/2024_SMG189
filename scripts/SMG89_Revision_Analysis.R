#!/usr/bin/env Rscript

# Title: SMG89_Revision_Analysis
# Objective: Top-level analysis script for "SMG1:SMG8:SMG9-complex integrity maintains robustness of nonsense-mediated mRNA decay"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(tidyverse)        # Generally required
library(ggh4x)            # Allows fixing row/column height!
library(TidyMultiqc)      # Required to load in multiQC data
library(ggdist)           # Required for certain geoms
library(scales)           # Transformation (log-scale)
library(viridis)          # Viridis color palette
library(patchwork)        # Allows easy plot stitching
library(ComplexHeatmap)   # For indeed more complex heatmaps and UpSet plots
library(ggridges)         # For Ridge plots!
library(ggcoverage)       # Produce fake IGV snapshots
library(ComplexHeatmap)   # Produce Heatmaps
library(circlize)         # Color for Heatmaps
library(gprofiler2)       # For GO analysis
library(tximeta)          # To import salmon counts for Heatmap-Clustering
library(fishpond)         # To scale salmon counts for Heatmap-Clustering
library(ggpmisc)          # To add linear regression to scatter plots
library(poolr)            # Combine p-values
library(nVennR)           # For overlaps with nVenns
library(GGally)           # For correlation plots
library(rstatix)          # Correlation testing - pipe-friendly
library(ggpubr)           # Add significance to boxplots
library(ggrepel)          # For nice labels in plots
library(tximport)
library(DESeq2)

# Set column width for plotting
cw1 = 5.5
cw2 = 12
cw3 = 18

# Define Arial as standard for ComplexHeatmap
pushViewport(viewport(gp = gpar(fontfamily = "Arial")))

# Define function for plotting "n"
n_fun <- function(x){
  return(data.frame(y = my_n,
                    label = paste0("n = ",
                                   length(x)
                    )))
}

# Load additional data sources --------------------------------------------

# Load experimentSet overview 

SMG89_datasets <- read_csv("~/2023_SMG89_project/SMG89_datasets.csv", 
                            trim_ws = TRUE) %>% 
  mutate(experiment_path = paste0(path,"/experiment.txt")) %>% 
  mutate(DESeq2_DGE_path = paste0(path,"/DESeq2/DGE/DESeq2_DGE_cat.csv")) %>% 
  mutate(Swish_DGE_path = paste0(path,"/Swish/DGE/Swish_DGE_cat.csv")) %>% 
  mutate(Swish_DTE_path = paste0(path,"/Swish/DTE/Swish_DTE_cat.csv")) %>% 
  mutate(edgeR_DTE_path = paste0(path,"/edgeR/DTE/edgeR_DTE_cat.csv")) %>% 
  # mutate(Swish_DTU_path = paste0(path,"/Swish/DTU/Swish_DTU_cat.csv")) %>% 
  mutate(ISAR_path = paste0(path,"/ISAR/SwitchList_filt_Analyzed.csv")) %>% 
  mutate(LeafCutter_path = paste0(path,"/leafcutter/leafcutter_AS_cat.csv")) %>% 
  mutate(experimentSet = fct_inorder(experimentSet)) %>% 
  mutate(publicationName = fct_inorder(publicationName))

# Load sample and condition mapping
SMG89_datasets_experiments <- readr::read_tsv(SMG89_datasets$experiment_path, 
                                             id = "experiment_path",
                                             col_names = c("sample", "condition_2")) %>% 
  left_join(SMG89_datasets %>% dplyr::select(-c(DESeq2_DGE_path, Swish_DGE_path, Swish_DTE_path, edgeR_DTE_path, ISAR_path))) %>% 
  arrange((experimentSet)) %>% 
  mutate(condition_2 = fct_inorder(condition_2))

# Load simplified gencode annotation
gtf_gencode_df_short <- read_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

# Load DESeq2 DGE data -----------------------------------------------------------

##
##  DESeq2 DGE  files   -----------------------------------------------------------
##

SMG89_DESeq2_DGE_combined <- readr::read_csv(SMG89_datasets %>% 
                                         pull(DESeq2_DGE_path),
                                       id = "DESeq2_DGE_path",
                                       col_select = (-1)) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other",
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(SMG89_datasets %>% dplyr::select(DESeq2_DGE_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-DESeq2_DGE_path)

# Load edgeR DTE data -----------------------------------------------------------

##
## edgeR DTE files   -----------------------------------------------------------
##

SMG89_edgeR_DTE_combined <- SMG89_datasets %>% 
  pull(edgeR_DTE_path) %>% 
  map_dfr(read_csv, col_select = (-1), id = 'edgeR_DTE_path') %>% 
  # filter(keep == TRUE) %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other",
                          transcript_type == "protein_coding" ~ "coding",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "lncRNA" ~ "lncRNA")) %>% 
  left_join(SMG89_datasets %>% dplyr::select(edgeR_DTE_path, experimentSet, publicationName)) %>% 
  mutate(condition_2 = fct_inorder(condition_2),
         experimentSet = fct_inorder(experimentSet),
         publicationName = fct_inorder(publicationName)) %>%
  dplyr::select(-edgeR_DTE_path)


# Save pre-parsed datasources - CSV ---------------------------------------------

SMG89_DESeq2_DGE_combined  %>%  write_csv(file.path("/home/volker/2023_SMG89_project", "SMG89_DESeq2_DGE_combined.csv"))
SMG89_edgeR_DTE_combined  %>%  write_csv(file.path("/home/volker/2023_SMG89_project", "SMG89_edgeR_DTE_combined.csv"))

# Save pre-parsed datasources - RDS ---------------------------------------------

## Save cluster information -------------------------------------------
save(SMG89_DESeq2_DGE_combined, 
     SMG89_edgeR_DTE_combined,
     file = paste0("/home/volker/2023_SMG89_project/",Sys.Date(),"_SMG189_datasources.rds"))

load("/home/volker/2023_SMG89_project/2024-10-28_SMG189_datasources.rds")

##
##
# Revision #1 -------------------------------------------------------------
##
##

##
## Figure 1 - KID ----------------------------------------------------------------
##

### DGE - Volcano/Density ---------------------------------------------------------------------

# Prepare for plotting
SMG89_delKID_DGE_for_plot <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 == "SMG8_delKID_0uM") %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange > 10, 10)) %>% 
  mutate(log2FoldChange = replace(log2FoldChange, log2FoldChange < -10, -10)) %>% 
  dplyr::filter(!is.na(padj))

SMG89_delKID_DGE_for_plot_ECDF <- SMG89_delKID_DGE_for_plot %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA"))

SMG89_delKID_DGE_for_plot_sig <- subset(SMG89_delKID_DGE_for_plot, padj < 0.0001 & abs(log2FoldChange) > 1) %>% 
  mutate(type = case_when(!gene_type %in% c("protein_coding", "lncRNA") ~ "other", 
                          gene_type == "protein_coding" ~ "coding",
                          gene_type == "lncRNA" ~ "lncRNA")) %>% 
  mutate(UpDown = case_when(log2FoldChange > 1 ~ "up",
                            log2FoldChange < -1 ~ "down"))

# Create unfiltered volcano plot
SMG89_delKID_DGE_plot <- ggplot(data = SMG89_delKID_DGE_for_plot, aes(x=log2FoldChange, y=-log10(padj))) + 
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 330), clip = 'on') +
  ggrastr::rasterise(geom_point(data=dplyr::filter(SMG89_delKID_DGE_for_plot, !gene_type %in% c("protein_coding", "lncRNA")), 
                                alpha=0.5,
                                size=0.75,
                                aes(color="#49C1AD")), dpi=600) + 
  ggrastr::rasterise(geom_point(data=dplyr::filter(SMG89_delKID_DGE_for_plot, gene_type == "protein_coding"), 
                                alpha=0.5,
                                size=0.75,
                                aes(color="#3E356B")), dpi=600) + 
  ggrastr::rasterise(geom_point(data=dplyr::filter(SMG89_delKID_DGE_for_plot, gene_type == "lncRNA"), 
                                alpha=0.5,
                                size=0.75,
                                aes(color="#357BA2")), dpi=600) + 
  geom_hline(yintercept=-log10(0.0001), linetype="dashed", linewidth=0.2) +
  geom_hline(aes(yintercept=320), linetype="dashed", color="gray", linewidth=0.2) +
  geom_vline(xintercept = 1, linetype="dashed", linewidth=0.2) +
  geom_vline(xintercept = -1, linetype="dashed", linewidth=0.2) +	
  geom_label_repel(aes(label=case_when(gene_name %in% c("GADD45A", "GADD45B", "GABARAPL1", "GAS5", "SNHG12", "ZFAS1") ~ gene_name,
                                       TRUE ~"")),
                   size = 4*0.36,
                   fill="lightgray",
                   max.overlaps = Inf,
                   label.padding = unit(0.1, "lines"),
                   label.size = 0,
                   parse=F) +
  labs(x = "log2 FoldChange",
       y = "-log10 padj") +
  scale_color_identity(name = "Gene type",
                       breaks = c("#3E356B", "#357BA2", "#49C1AD"),
                       labels = c("coding", "lncRNA", "other"),
                       guide = "legend") 

SMG89_delKID_DGE_plot <- SMG89_delKID_DGE_plot + annotate("text", x=-7.5, y=335, label="Max p-value", size = 6/.pt, color = "gray") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols =  unit(45, "mm"))

SMG89_delKID_DGE_for_plot_sig_summary <- SMG89_delKID_DGE_for_plot_sig %>% 
  group_by(type, UpDown) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = UpDown, values_from = n) %>% 
  arrange(factor(type, levels = c('coding', 'lncRNA', 'other')))

# Create filtered density plot
SMG89_delKID_DGE_plot2 <- ggplot(data = SMG89_delKID_DGE_for_plot_sig, aes(log2FoldChange)) + 
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +	
  theme(legend.position="none") +
  xlim(-10, 10) +
  geom_density(data=dplyr::filter(SMG89_delKID_DGE_for_plot_sig, !gene_type %in% c("protein_coding", "lncRNA")), color="#49C1AD", fill="#49C1AD", alpha = 0.2) +
  geom_density(data=dplyr::filter(SMG89_delKID_DGE_for_plot_sig, gene_type == "protein_coding"), color="#3E356B", fill="#3E356B", alpha = 0.2) +
  geom_density(data=dplyr::filter(SMG89_delKID_DGE_for_plot_sig, gene_type == "lncRNA"), color="#357BA2", fill="#357BA2", alpha = 0.2) +
  annotate(geom = "table", x = -Inf, y = Inf, label = (SMG89_delKID_DGE_for_plot_sig_summary), 
           table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                     core=list(fg_params=list(col = c("#3E356B", "#357BA2", "#49C1AD"),fontface=1)),
                                     colhead=list(fg_params=list( fontface=1))),
           vjust = 1, hjust = -0.1) +
  labs(x = "log2 FoldChange", y = "filtered density", caption = paste0("\nFilter: |log2FoldChange| > 1 & p.adjust < 0.0001")) +
  force_panelsizes(rows = unit(20, "mm"),
                   cols =  unit(45, "mm"))


# Save plot
ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_delKID_DGE_finalplot_Volcano_", Sys.Date(), ".pdf")),
       SMG89_delKID_DGE_plot,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Save plot
ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_delKID_DGE_finalplot_Density_", Sys.Date(), ".pdf")),
       SMG89_delKID_DGE_plot2,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### GO DGE  -----------------------------------------------------------------

# Background
GO_DGE_SMG89_delKID <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 == "SMG8_delKID_0uM") %>% 
  distinct(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Define up/down regulated classes
SMG89_delKID_DGE_for_plot_sig <- SMG89_delKID_DGE_for_plot %>% 
  mutate(DGE_class = case_when(padj < 0.0001 & log2FoldChange > 1 ~ "up",
                               padj < 0.0001 & log2FoldChange < -1 ~ "down",
                               TRUE ~ "n.s."))

# Up-regulated - Ordered by padj!
SMG89_delKID_DGE_for_plot_sig_up <- SMG89_delKID_DGE_for_plot_sig %>% 
  filter(DGE_class == "up") %>% 
  arrange((padj)) %>% 
  dplyr::select(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Up-regulated - Ordered by padj!
SMG89_delKID_DGE_for_plot_sig_down <- SMG89_delKID_DGE_for_plot_sig %>% 
  filter(DGE_class == "down") %>% 
  arrange((padj)) %>% 
  dplyr::select(gene_id) %>% 
  separate(gene_id, c("gene_id", "Version")) %>% 
  pull(gene_id)

# Up GO
gostres_SMG89_delKID_up_GO_list <- gost(query = SMG89_delKID_DGE_for_plot_sig_up,
                                        custom_bg = GO_DGE_SMG89_delKID,
                                        sources = "GO:BP",
                                        domain_scope = "custom_annotated",
                                        organism = "hsapiens",
                                        ordered_query = TRUE,
                                        correction_method = c("fdr"),
                                        evcodes = TRUE,
                                        # as_short_link = TRUE,
                                        significant = FALSE)

gostres_SMG89_delKID_up_GO_list_result <- gostres_SMG89_delKID_up_GO_list$result

gostres_SMG89_delKID_up_GO_list_result %>% 
  dplyr::count()

gostres_SMG89_delKID_up_GO_list_result %>% 
  arrange(p_value) %>% 
  mutate(term_id = fct_inorder(term_id)) %>% 
  ggplot(aes(y=-log10(p_value),
             x=term_id)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  coord_cartesian(ylim = c(0, 2), clip = 'off') +
  scale_x_discrete(expand=c(0.1,0.1)) +
  ggrastr::rasterise(geom_point(size=0.25, alpha=0.5), dpi=600, dev = "cairo") +
  geom_hline(yintercept = -log10(0.05),
             color="darkred") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="GO:BP terms",
       y="-log10(p-value)") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols =  unit(20, "mm"))

# Save plot
ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_delKID_DGE_up_GO", Sys.Date(), ".pdf")),
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Down GO
gostres_SMG89_delKID_down_GO_list <- gost(query = SMG89_delKID_DGE_for_plot_sig_down,
                                          custom_bg = GO_DGE_SMG89_delKID,
                                          sources = "GO:BP",
                                          domain_scope = "custom_annotated",
                                          organism = "hsapiens",
                                          ordered_query = TRUE,
                                          correction_method = c("fdr"),
                                          evcodes = TRUE,
                                          # as_short_link = TRUE,
                                          significant = FALSE)

gostres_SMG89_delKID_down_GO_list_result <- gostres_SMG89_delKID_down_GO_list$result

gostres_SMG89_delKID_down_GO_list_result %>% 
  dplyr::count()

gostres_SMG89_delKID_down_GO_list_result %>% 
  arrange(p_value) %>% 
  mutate(term_id = fct_inorder(term_id)) %>% 
  ggplot(aes(y=-log10(p_value),
             x=term_id)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=6), 
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  coord_cartesian(ylim = c(0, 2), clip = 'off') +
  scale_x_discrete(expand=c(0.1,0.1)) +
  ggrastr::rasterise(geom_point(size=0.25, alpha=0.5), dpi=600, dev = "cairo") +
  geom_hline(yintercept = -log10(0.05),
             color="darkred") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x="GO:BP terms",
       y="-log10(p-value)") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols =  unit(20, "mm"))

# Save plot
ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_delKID_DGE_down_GO", Sys.Date(), ".pdf")),
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### DTE ECDF ----------------------------------------------------------------

SMG89_edgeR_DTE_deltaKID_sig <- SMG89_edgeR_DTE_combined %>%
  filter(condition_2 == "SMG8_delKID_0uM") %>% 
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  filter(abs(log2FC) > 1 & qvalue < 0.0001) %>% 
  mutate(UpDown = case_when(log2FC > 1 ~ "up",
                            log2FC < -1 ~ "down")) %>% 
  dplyr::count(type,UpDown) %>% 
  pivot_wider(names_from = UpDown, values_from = n) %>% 
  arrange(factor(type, levels = c('coding', 'lncRNA', "NMD", 'other')))

# Kolmogorov-Smirnov Test
ks_SMG89_edgeR_DTE_deltaKID_sig <- ks.test(SMG89_edgeR_DTE_combined %>%
                                             filter(condition_2 == "SMG8_delKID_0uM") %>% 
                                             dplyr::filter(type == "coding") %>% 
                                             dplyr::pull(logFC), 
                                           SMG89_edgeR_DTE_combined %>%
                                             filter(condition_2 == "SMG8_delKID_0uM") %>%
                                             dplyr::filter(type == "NMD") %>% 
                                             dplyr::pull(logFC))

ks_SMG89_edgeR_DTE_deltaKID_sig$p.value <- ifelse(ks_SMG89_edgeR_DTE_deltaKID_sig$p.value == 0, "< 2.2e-16", signif(ks_SMG89_edgeR_DTE_deltaKID_sig$p.value,3))

SMG89_edgeR_DTE_combined %>%
  filter(condition_2 == "SMG8_delKID_0uM") %>% 
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  ggplot(aes(x=log2FC, color=type)) + 
  theme_classic() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.1,0.1,0.1,0.1)),
        legend.position = "right",
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_vline(xintercept = 0, linetype="solid", color="gray") +	
  geom_vline(xintercept = 1, linetype="dashed", color="darkgray") +	
  geom_vline(xintercept = -1, linetype="dashed", color="darkgray") +	
  ggrastr::rasterise(stat_ecdf(linewidth=0.75,
                               alpha=0.75,
                               geom = "line"), dpi=600, dev = "cairo") +
  scale_color_manual(values=c("coding" = "#3E356B",
                              "NMD" = "darkred",
                              "lncRNA" = "#357BA2",
                              "other" = "darkgray")) +
  scale_x_continuous(limits = c(-3.5, 3.5),
                     breaks = c(-2.5, -1, 0, 1, 2.5)) +
  annotate(geom = "table", x = Inf, y = -Inf, label = (SMG89_edgeR_DTE_deltaKID_sig), 
           table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                     core=list(fg_params=list(col = c("#3E356B", "#357BA2", "darkred", "#49C1AD"),fontface=1)),
                                     colhead=list(fg_params=list( fontface=1))),
           vjust = 1, hjust = -0.1) +
  annotate(geom ="text", x = -1.25, y = 0.6, size = 5*0.36, label = paste0("KS test \n (coding vs NMD) \n p-value: ",
                                                                           paste(ks_SMG89_edgeR_DTE_deltaKID_sig$p.value),
                                                                           "\n D: ",
                                                                           paste(signif(ks_SMG89_edgeR_DTE_deltaKID_sig$statistic,3))),
           vjust = -0.25, hjust = 1) +
  labs(x = "log2 FoldChange",
       y = "ECDF",
       color="Condition") +
  force_panelsizes(rows = unit(30, "mm"),
                   cols = unit(45, "mm")) +
  guides(color=guide_legend(nrow=4,byrow=TRUE))

# Save plot
ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_delKID_DTE_ECDF_", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### Transcript properties ---------------------------------------------------

gtf_gencode_v42_cluster_length_GC_mfe <- read_csv(file.path("/home/volker/2023_UPF1project", "gtf_gencode_v42_cluster_length_GC_mfe.csv"))

# Define up/down regulated classes
SMG89_edgeR_DTE_deltaKID_sigClass <- SMG89_edgeR_DTE_combined %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  filter(condition_2 == "SMG8_delKID_0uM") %>% 
  mutate(DTE_class = case_when(qvalue < 0.0001 & log2FC > 1 ~ "up",
                               qvalue < 0.0001 & log2FC < -1 ~ "down",
                               TRUE ~ "n.s."))

SMG89_edgeR_DTE_deltaKID_sigClass_prop <- SMG89_edgeR_DTE_deltaKID_sigClass %>% 
  left_join(gtf_gencode_v42_cluster_length_GC_mfe) %>% 
  dplyr::select(transcript_id,
                DTE_class,
                nexon,
                tx_len,
                cds_len,
                utr5_len,
                utr3_len,
                GC_tx,
                GC_utr5,
                GC_cds,
                GC_utr3,
                utr3_mfe_nt) %>% 
  mutate(DTE_class = fct_rev(fct_relevel(DTE_class, 
                                         "up",
                                         "n.s.",
                                         "down"))) %>% 
  pivot_longer(cols = -c(transcript_id,
                         DTE_class),
               names_to = "variable",
               values_to = "values") %>% 
  mutate(variable = case_when(variable == "nexon" ~ "number of exons (n)",
                              variable == "tx_len" ~ "transcript length (nt)",
                              variable == "utr5_len" ~ "5′UTR length (nt)",
                              variable == "cds_len" ~ "CDS length (nt)",
                              variable == "utr3_len" ~ "3′UTR length (nt)",
                              variable == "GC_tx" ~ "transcript GC%",
                              variable == "GC_utr5" ~ "5′UTR GC%",
                              variable == "GC_cds" ~ "CDS GC%",
                              variable == "GC_utr3" ~ "3′UTR GC%",
                              variable == "utr3_mfe_nt" ~ "3′UTR –ΔG/nt")) %>% 
  mutate(variable = fct_relevel(variable,
                                "transcript length (nt)",
                                "5′UTR length (nt)",
                                "CDS length (nt)",
                                "3′UTR length (nt)",
                                "transcript GC%",
                                "5′UTR GC%",
                                "CDS GC%",
                                "3′UTR GC%",
                                "number of exons (n)",
                                "3′UTR –ΔG/nt")) 


# Fix for boxplot scales
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

SMG89_edgeR_DTE_deltaKID_sigClass_prop 


SMG89_edgeR_DTE_deltaKID_sigClass_prop %>% 
  ggplot(aes(x=DTE_class,
             y=values,
             fill=DTE_class)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = 'darkgray', linewidth = 0.1),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        strip.text = element_text(size = 6),
        text=element_text(family="Arial")) +
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", fatten=2,
               linewidth=0.25) +
  scale_fill_manual(values=c("up" = "#b2182b",
                             "n.s." = "gray50",
                             "down" = "#2166ac")) +
  facet_wrap(~variable, ncol=3, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x="",
       y="",
       fill="Sig. regulated\ntranscript class") +
  force_panelsizes(rows = unit(3*3+0.6, "mm"),
                   cols = unit(6*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_delKID_DTE_properties_boxplot.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

# Significance testing using Dunn Test (https://rpkgs.datanovia.com/rstatix/reference/dunn_test.html)
SMG89_edgeR_DTE_deltaKID_sigClass_prop_dunn <- SMG89_edgeR_DTE_deltaKID_sigClass_prop %>% 
  group_by(variable) %>% 
  dunn_test(values ~ DTE_class,
            detailed = TRUE)

# Effect sizes using wilcox effsize (https://rpkgs.datanovia.com/rstatix/reference/wilcox_effsize.html)
SMG89_edgeR_DTE_deltaKID_sigClass_prop_effSize <- SMG89_edgeR_DTE_deltaKID_sigClass_prop %>% 
  group_by(variable) %>% 
  wilcox_effsize(values ~ DTE_class)

# Combine both significance and effect sizes
SMG89_edgeR_DTE_deltaKID_sigClass_combined <- SMG89_edgeR_DTE_deltaKID_sigClass_prop_dunn %>% 
  add_x_position(x = "DTE_class") %>% 
  left_join(SMG89_edgeR_DTE_deltaKID_sigClass_prop_effSize)

# Extract relevant information
SMG89_edgeR_DTE_deltaKID_sigClass_combined_Exp <- as_tibble(SMG89_edgeR_DTE_deltaKID_sigClass_combined) %>% 
  dplyr::filter(group1 == "n.s." | group2 == "n.s.") %>% 
  mutate(values = "A")  %>% 
  mutate(DTE_class = case_when(group1 == "n.s." ~ group2,
                               group1 != "n.s." ~ group1)) %>% 
  dplyr::select(DTE_class,
                values,
                variable,
                p.adj,
                p.adj.signif,
                effsize,
                magnitude) %>% 
  arrange(p.adj)

##
## Figure 2 - KOs ----------------------------------------------------------
##


### DGE - overview -----------------------------------------------------------


#### Prepare for Plot -----------------------------------------------------------

SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                              "UPF1_Nter_12h",
                              "UPF1_FKBP_HCT_12h",
                              "UPF1_FKBP_HEK_12h",
                              "SMG8_KO_0uM",
                              "SMG9_KO_0uM",
                              "SMG8_delKID_0uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                   "SMG7_KO_34_SMG6_KD",
                                   "UPF1_FKBP_HEK_12h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_Nter_12h",
                                   "SMG8_KO_0uM",
                                   "SMG9_KO_0uM",
                                   "SMG8_delKID_0uM"))) %>% 
  mutate(significant = case_when(padj < 0.0001 & abs(log2FoldChange) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(log2FoldChange>0 ~ "up",
                             log2FoldChange<0 ~ "down")) %>% 
  # mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_gene = length(unique(gene_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, type) %>% 
  mutate(detected_gene_type = length(unique(gene_id))) %>% 
  ungroup()

##### Rev_1_F2_Barplot -----------------------------------------------------------

SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, type, up_down, detected_gene_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(type_up_down = paste0(type, "-", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding-up",
                                     "coding-down",
                                     "lncRNA-up",
                                     "lncRNA-down",
                                     "other-up",
                                     "other-down"))) %>% 
  mutate(n = case_when(up_down == "down" ~ -n,
                       up_down == "up" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n,
             fill=type)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  geom_col(color="black",
           linewidth = 0.2) +
  scale_fill_viridis_d(option="G",
                       begin=0.25,
                       end=0.75) +
  labs(y="",
       x="Number of \nsig. regulated genes",
       fill="Gene biotype",
       size="Sig. DGE events") +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot$condition_2))*3+0.4, "mm"),
                   cols = unit(20, "mm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE,
                           title.position = "top"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure2_DGE_Barplot_perType", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### Rev_1_F2_Examples -----------------------------------------------------------

# Define representative NMD-targeted genes
standard_heatmap_genes <- c("ZFAS1",
                            "SNHG12",
                            "GAS5",
                            "GADD45B",
                            "GADD45A")

SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot %>%  
  filter(gene_name %in% standard_heatmap_genes) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, gene_name) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "ZFAS1",
                                 "SNHG12",
                                 "GAS5",
                                 "GADD45B",
                                 "GADD45A")) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=gene_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="DGE log2FC") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                # barwidth = unit(8, "mm"),
                                # barheight = unit(2, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "bottom")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot$condition_2))*3+0.6, "mm"),
                   cols = unit(5*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure2_DGE_Dotplot_Examples", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### DTE - overview -----------------------------------------------------------


#### Prepare for Plot -----------------------------------------------------------

SMG89_Rev_1_F2_DESeq2_DTE_combined_forPlot <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM"))) %>% 
  dplyr::rename("log2FoldChange" = "logFC",
                "padj" = "FDR") %>% 
  mutate(significant = case_when(padj < 0.0001 & abs(log2FoldChange) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(log2FoldChange>0 ~ "up",
                             log2FoldChange<0 ~ "down")) %>% 
  # mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_transcript = length(unique(transcript_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, type) %>% 
  mutate(detected_transcript_type = length(unique(transcript_id))) %>% 
  ungroup()


##### Rev1_F2_Lollipop --------------------------------------------------------

SMG89_KO_Lollipop_DTE_median_log2FC <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM")) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM"))) %>% 
  filter(type %in% c("NMD")) %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  # filter(!type %in% c("other","lncRNA")) %>% 
  group_by(condition_2, type) %>% 
  summarize(median_log2FC=median(log2FC, na.rm = TRUE)) 

SMG89_KO_Lollipop_DTE_sig <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM")) %>% 
  filter(type %in% c("NMD")) %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  # filter(!type %in% c("other","lncRNA")) %>%
  filter(qvalue < 0.0001 & abs(log2FC > 1)) %>% 
  mutate(up_down = case_when(log2FC > 1 ~ "up",
                             log2FC < -1 ~ "down")) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM"))) %>% 
  group_by(condition_2, type) %>% 
  summarize(n=n()) %>% 
  left_join(SMG89_KO_Lollipop_DTE_median_log2FC) %>% 
  # mutate(n = case_when(type == "coding" ~ -n,
  # type == "NMD" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  # geom_linerange(aes(xmin=log2FC_lower,
  #                    xmax=log2FC_upper),
  #                position = position_dodge(width = 0.9),
  #                color="black",
  #                linewidth=0.25) +
  geom_linerange(aes(xmin=0, xmax=n),
                 linewidth=0.25,
                 position = position_dodge(width=0.75)) +
  geom_point(aes(
    # size=abs(median_log2FC),
    # shape=case_when(median_log2FC > 0 ~ "positive",
    #                 median_log2FC < 0 ~ "negative"),
    fill=abs(median_log2FC)),
    shape=21,
    size=3,
    position = position_dodge(width=0.75),
    stroke=0.2) +
  # scale_shape_manual(values = c("positive" = 21,
  #                               "negative" = 23)) +
  # theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_fill_manual(values=c("coding" = "#6A00A8FF",
  #                            "NMD" = "#60CEACFF")) +
  # scale_x_continuous(limits=c(-4,4)) +
  # scale_size(limits = c(0,1.5),
  #            range = c(0.25,3)) +
  scale_fill_gradient(low = ("white"),
                      high = ("#B2182B"),
                      na.value = "grey80") +
  labs(y="",
       x="Number of\nsig. regulated\nNMD-annotated transcripts",
       fill="Median\nDTE\nlog2FC",
       shape="Median log2FC",
       size="Median log2FC") +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F2_DESeq2_DGE_combined_forPlot$condition_2))*3+0.4, "mm"),
                   cols = unit(20, "mm")) +
  scale_x_continuous(expand = expansion(c(0.15, 0.15)),
                     breaks = c(0, 1500, 3000)) +
  # theme(axis.text.y=element_blank(),
  #       axis.ticks.y = element_blank()) +
  guides(size=guide_legend(nrow=2,byrow=TRUE),
         fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure2_DTE_Lollipop", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### Rev_1_F2_Examples -----------------------------------------------------------

# Define representative NMD-targeted transcripts
standard_heatmap_transcripts <- c("TMEM222-202",
                                  "GABARAPL1-208",
                                  "SRSF2-208",
                                  "SRSF2-204",
                                  "RPL3-204")

SMG89_Rev_1_F2_DESeq2_DTE_combined_forPlot %>%  
  filter(transcript_name %in% standard_heatmap_transcripts) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, transcript_name) %>% 
  mutate(transcript_name = fct_relevel(transcript_name,
                                       "TMEM222-202",
                                       "SRSF2-208",
                                       "SRSF2-204",
                                       "RPL3-204",
                                       "GABARAPL1-208"
                                       )) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=transcript_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                # barwidth = unit(8, "mm"),
                                # barheight = unit(2, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "bottom")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F2_DESeq2_DTE_combined_forPlot$condition_2))*3+0.6, "mm"),
                   cols = unit(5*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure2_DTE_Dotplot_Examples", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Compare 8-KO vs 9-KO ----------------------------------------------------

#### DGE ---------------------------------------------------------------------

SMG89_DESeq2_DGE_significant_genes <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(padj < 0.0001) %>%
  pull(gene_id)
  

##### nVenn --------------------------------------------------------

SMG89_DESeq2_DGE_set1 <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 == "SMG8_KO_0uM") %>%
  filter(padj < 0.0001) %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  pull(gene_id_plain)

SMG89_DESeq2_DGE_set2 <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 == "SMG9_KO_0uM") %>%
  filter(padj < 0.0001) %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  pull(gene_id_plain)

SMG89_DGE_sig_myNV_up <- plotVenn(list(SMG89_DESeq2_DGE_set1,
                                       SMG89_DESeq2_DGE_set2), 
                                  sNames=c("SMG8", "SMG9"), showPlot = FALSE)

showSVG(SMG89_DGE_sig_myNV_up,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#09315C", "#51A8DB"))

##### Scatter plot ------------------------------------------------------------

SMG89_KO_comparison_numbers <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(padj < 0.0001) %>% 
  dplyr::select(gene_id, gene_name, type, log2FoldChange, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = log2FoldChange) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG9_KO_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG9_KO_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM < 0 ~ "both_down",
                               TRUE ~ "")) %>% 
  filter(direction != "") %>% 
  dplyr::count(direction) %>% 
  mutate(direction = (fct_relevel(direction,
                                  "both_up",
                                  "both_down",
                                  "diverging"))) %>% 
  arrange(direction)

# Plot
SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(padj < 0.0001) %>% 
  dplyr::select(gene_id, gene_name, type, log2FoldChange, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = log2FoldChange) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG9_KO_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG9_KO_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM < 0 ~ "both_down")) %>% 
  mutate(direction = fct_relevel(direction,
                                 "both_up",
                                 "both_down",
                                 "diverging")) %>% 
  ggplot(aes(x=SMG8_KO_0uM,
             y=SMG9_KO_0uM)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(-10,10)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  scale_fill_manual(values=c("both_up" = "#B2182B",
                             "both_down" = "#2166AC",
                             "diverging" = "gray70")) +
  ggrastr::rasterise(geom_point(aes(fill=direction),
                                shape=21,
                                size=1.25), dpi=600) +
  stat_poly_line(show.legend=FALSE,
                 color="gray50") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "p.value.label")),
               size = 5*0.36,
               color="gray20",
               label.y = "top", label.x = "right") +
  labs(subtitle= "Significant DGE events (padj < 0.0001)",
       x="log2FC(SMG8-KO)",
       y="log2FC(SMG9-KO)") +
  annotate(geom = "table", x = Inf, y = -5, label = (SMG89_KO_comparison_numbers), 
           table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                     core=list(fg_params=list(col = c("#B2182B", "#2166AC", "gray40"),fontface=1)),
                                     colhead=list(fg_params=list( fontface=1))),
           vjust = 1, hjust = -0.1) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_KO_DGE_comparison_scatter", Sys.Date(), ".pdf")),
       width = cw2+1,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### DTE ---------------------------------------------------------------------

#### nVenn --------------------------------------------------------

SMG89_edgeR_DTE_set1 <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 == "SMG8_KO_0uM") %>%
  filter(FDR < 0.0001) %>% 
  pull(transcript_id)

SMG89_edgeR_DTE_set2 <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 == "SMG9_KO_0uM") %>%
  filter(FDR < 0.0001) %>% 
  pull(transcript_id)

SMG89_DTE_sig_myNV_up <- plotVenn(list(SMG89_edgeR_DTE_set1,
                                       SMG89_edgeR_DTE_set2), 
                                  sNames=c("SMG8", "SMG9"), showPlot = FALSE)

showSVG(SMG89_DTE_sig_myNV_up,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#09315C", "#51A8DB"))


#### Scatter plot ------------------------------------------------------------

SMG89_KO_DTE_comparison_numbers <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(FDR < 0.0001) %>% 
  dplyr::select(transcript_id, transcript_name, gene_name, type, logFC, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = logFC) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG9_KO_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG9_KO_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM < 0 ~ "both_down",
                               TRUE ~ "")) %>% 
  filter(direction != "") %>% 
  dplyr::count(direction) %>% 
  mutate(direction = (fct_relevel(direction,
                                  "both_up",
                                  "both_down",
                                  "diverging"))) %>% 
  arrange(direction)

SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(FDR < 0.0001) %>% 
  dplyr::select(transcript_id, transcript_name, gene_name, type, logFC, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = logFC) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG9_KO_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG9_KO_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG9_KO_0uM < 0 ~ "both_down")) %>% 
  mutate(direction = fct_relevel(direction,
                                 "both_up",
                                 "both_down",
                                 "diverging")) %>% 
  ggplot(aes(x=SMG8_KO_0uM,
             y=SMG9_KO_0uM)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(-10,10)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  scale_fill_manual(values=c("both_up" = "#B2182B",
                             "both_down" = "#2166AC",
                             "diverging" = "gray70")) +
  ggrastr::rasterise(geom_point(aes(fill=direction),
                                shape=21,
                                size=1.25), dpi=600) +
  stat_poly_line(show.legend=FALSE,
                 color="gray50") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "p.value.label")),
               size = 5*0.36,
               color="gray20",
               label.y = "top", label.x = "right") +
  labs(subtitle= "Significant DTE events (padj < 0.0001)",
       x="log2FC(SMG8-KO)",
       y="log2FC(SMG9-KO)") +
  annotate(geom = "table", x = Inf, y = -5, label = (SMG89_KO_DTE_comparison_numbers), 
           table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                     core=list(fg_params=list(col = c("#B2182B", "#2166AC", "gray40"),fontface=1)),
                                     colhead=list(fg_params=list( fontface=1))),
           vjust = 1, hjust = -0.1) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_KO_DTE_comparison_scatter", Sys.Date(), ".pdf")),
       width = cw2+1,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


#Test
SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG9_KO_0uM")) %>% 
  filter(FDR < 0.0001) %>% 
  arrange(desc(logFC))


## Compare 8-KO vs delKID ----------------------------------------------------

### DGE ---------------------------------------------------------------------

#### nVenn --------------------------------------------------------

SMG89_DESeq2_DGE_KID_set1 <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 == "SMG8_KO_0uM") %>%
  filter(padj < 0.0001) %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  pull(gene_id_plain)

SMG89_DESeq2_DGE_KID_set2 <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 == "SMG8_delKID_0uM") %>%
  filter(padj < 0.0001) %>% 
  separate(gene_id, c("gene_id_plain", NA)) %>% 
  pull(gene_id_plain)

SMG89_DGE_KID_sig_myNV_up <- plotVenn(list(SMG89_DESeq2_DGE_KID_set1,
                                           SMG89_DESeq2_DGE_KID_set2), 
                                      sNames=c("SMG8", "delKID"), showPlot = FALSE)

showSVG(SMG89_DGE_KID_sig_myNV_up,
        opacity=0.3,
        borderWidth = 3,
        systemShow = TRUE,
        labelRegions = F,
        fontScale = 1.5,
        setColors = c("#09315C", "#51A8DB"))

SMG89_KID_comparison_numbers <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG8_delKID_0uM")) %>% 
  filter(padj < 0.0001) %>% 
  dplyr::select(gene_id, gene_name, type, log2FoldChange, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = log2FoldChange) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG8_delKID_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG8_delKID_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG8_delKID_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG8_delKID_0uM < 0 ~ "both_down",
                               TRUE ~ "")) %>% 
  filter(direction != "") %>% 
  dplyr::count(direction) %>% 
  mutate(direction = (fct_relevel(direction,
                                  "both_up",
                                  "both_down",
                                  "diverging"))) %>% 
  arrange(direction)

SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG8_delKID_0uM")) %>% 
  filter(padj < 0.0001) %>% 
  dplyr::select(gene_id, gene_name, type, log2FoldChange, condition_2) %>% 
  pivot_wider(names_from = condition_2,
              values_from = log2FoldChange) %>% 
  mutate(direction = case_when(SMG8_KO_0uM > 0 & SMG8_delKID_0uM > 0 ~ "both_up",
                               SMG8_KO_0uM < 0 & SMG8_delKID_0uM > 0 ~ "diverging",
                               SMG8_KO_0uM > 0 & SMG8_delKID_0uM < 0 ~ "diverging",
                               SMG8_KO_0uM < 0 & SMG8_delKID_0uM < 0 ~ "both_down")) %>% 
  mutate(direction = fct_relevel(direction,
                                 "both_up",
                                 "both_down",
                                 "diverging")) %>% 
  ggplot(aes(x=SMG8_KO_0uM,
             y=SMG8_delKID_0uM)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  scale_x_continuous(limits = c(-10,10)) +
  scale_y_continuous(limits = c(-10,10)) +
  geom_hline(yintercept=0,
             linewidth=0.2) +
  geom_vline(xintercept=0,
             linewidth=0.2) +
  scale_fill_manual(values=c("both_up" = "#B2182B",
                             "both_down" = "#2166AC",
                             "diverging" = "gray70")) +
  ggrastr::rasterise(geom_point(aes(fill=direction),
                                shape=21,
                                size=1.25), dpi=600) +
  stat_poly_line(show.legend=FALSE,
                 color="gray50") +
  stat_poly_eq(mapping = use_label(c("adj.R2", "p.value.label")),
               size = 5*0.36,
               color="gray20",
               label.y = "top", label.x = "right") +
  labs(subtitle= "Significant DGE events (padj < 0.0001)",
       x="log2FC(SMG8-KO)",
       y="log2FC(SMG8-delKID)") +
  annotate(geom = "table", x = Inf, y = -5, label = (SMG89_KID_comparison_numbers), 
           table.theme = ttheme_gtbw(base_size = 5, base_colour = "black", base_family = "Arial", padding = unit(c(1, 1), "mm"), 
                                     core=list(fg_params=list(col = c("#B2182B", "#2166AC", "gray40"),fontface=1)),
                                     colhead=list(fg_params=list( fontface=1))),
           vjust = 1, hjust = -0.1) +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_8KO_delKID_DGE_comparison_scatter", Sys.Date(), ".pdf")),
       width = cw2+1,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##
## Figure 5 - KOs + SMG1i ----------------------------------------------------------
##

### PCA plot -----------------------------------------------------------

# Define directory
mydir = "/home/volker/gencode.v42.datasets/2023_11_SMG89_SK"

# Indicate reference folder in which the tx2gene file is located
ref_dir="/home/volker/reference"

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

dds = dds_string

# Extracting count data transformed values
vsd <- vst(dds, blind=FALSE)

# Get sample-to-sample distance
sampleDists <- dist(t(assay(vsd)))

# Principal component plot of the samples
plt <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(plt, "percentVar"))

plt_plot <- ggplot(plt %>% mutate(condition =(fct_inorder(as_factor(condition)))), aes(PC1, PC2, fill=condition)) +
  theme_classic() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        strip.text.x = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  geom_point(shape=21,
             alpha=0.75,
             color="black",
             size=1) +
  # ggrepel::geom_text_repel(data=plt %>% distinct(condition, .keep_all = TRUE),
  #                          aes(label=condition),
  #                          box.padding = 0.25,
  #                          min.segment.length = Inf,
  #                          size=5*0.36) +
  scale_fill_manual(values=c("control" = "#000000",
                             "control_01uM" = "#D3D3D3",
                             "control_1uM" = "#A9A9A9",
                             "SMG8_KO_0uM" = "#8c96c6",
                             "SMG8_KO_01uM" = "#8856a7",
                             "SMG8_KO_1uM" = "#810f7c",
                             "SMG9_KO_0uM" = "#7bccc4",
                             "SMG9_KO_01uM" = "#43a2ca",
                             "SMG9_KO_1uM" = "#0868ac",
                             "SMG8_delKID_0uM" = "#fc8d59",
                             "SMG8_delKID_01uM" = "#e34a33",
                             "SMG8_delKID_1uM" = "#b30000")) +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  force_panelsizes(rows = unit(40, "mm"),
                   cols =  unit(40, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", "A3_DESeq2_DGE_PCA.pdf"),
       plt_plot,
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


### DGE - overview -----------------------------------------------------------


#### Prepare for Plot -----------------------------------------------------------

SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "control_01uM",
                            "control_1uM",
                            "SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  mutate(significant = case_when(padj < 0.0001 & abs(log2FoldChange) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(log2FoldChange>0 ~ "up",
                             log2FoldChange<0 ~ "down")) %>% 
  # mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_gene = length(unique(gene_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, type) %>% 
  mutate(detected_gene_type = length(unique(gene_id))) %>% 
  ungroup()

##### Rev_1_F2_Barplot -----------------------------------------------------------

SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot %>% 
  filter(significant == TRUE) %>% 
  group_by(condition_2, type, up_down, detected_gene_type) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  mutate(type_up_down = paste0(type, "-", up_down)) %>% 
  mutate(type_up_down = (fct_relevel(type_up_down,
                                     "coding-up",
                                     "coding-down",
                                     "lncRNA-up",
                                     "lncRNA-down",
                                     "other-up",
                                     "other-down"))) %>% 
  mutate(n = case_when(up_down == "down" ~ -n,
                       up_down == "up" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n,
             fill=type)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  geom_col(color="black",
           linewidth = 0.2) +
  scale_x_continuous(breaks = c(-2000, 0, 2000)) +
  scale_fill_viridis_d(option="G",
                       begin=0.25,
                       end=0.75) +
  labs(y="",
       x="Number of \nsig. regulated genes",
       fill="Gene biotype",
       size="Sig. DGE events") +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot$condition_2))*3+0.6, "mm"),
                   cols = unit(20, "mm")) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE,
                           title.position = "top"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure5_DGE_Barplot_perType", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### Rev_1_F2_Examples -----------------------------------------------------------

# Define representative NMD-targeted genes
standard_heatmap_genes <- c("ZFAS1",
                            "SNHG12",
                            "GAS5",
                            "GADD45B",
                            "GADD45A")

SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot %>%  
  filter(gene_name %in% standard_heatmap_genes) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, gene_name) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "ZFAS1",
                                 "SNHG12",
                                 "GAS5",
                                 "GADD45B",
                                 "GADD45A")) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=gene_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="DGE log2FC") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                # barwidth = unit(8, "mm"),
                                # barheight = unit(2, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "bottom")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot$condition_2))*3+0.6, "mm"),
                   cols = unit(5*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure5_DGE_Dotplot_Examples", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### DTE - overview -----------------------------------------------------------


#### Prepare for Plot -----------------------------------------------------------

SMG89_Rev_1_F5_DESeq2_DTE_combined_forPlot <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM",
                            "control_01uM",
                            "SMG8_KO_01uM",
                            "SMG9_KO_01uM",
                            "SMG8_delKID_01uM",
                            "control_1uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_1uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_inorder(condition_2)) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  dplyr::rename("log2FoldChange" = "logFC",
                "padj" = "FDR") %>% 
  mutate(significant = case_when(padj < 0.0001 & abs(log2FoldChange) > 1 ~ TRUE,
                                 TRUE ~ FALSE)) %>% 
  mutate(up_down = case_when(log2FoldChange>0 ~ "up",
                             log2FoldChange<0 ~ "down")) %>% 
  # mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  group_by(experimentSet) %>% 
  mutate(detected_transcript = length(unique(transcript_id))) %>% 
  ungroup() %>% 
  group_by(experimentSet, type) %>% 
  mutate(detected_transcript_type = length(unique(transcript_id))) %>% 
  ungroup()


##### Rev_1_F2_Lollipop --------------------------------------------------------

SMG89_SMG1i_KO_Lollipop_DTE_median_log2FC <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM",
                            "control_01uM",
                            "SMG8_KO_01uM",
                            "SMG9_KO_01uM",
                            "SMG8_delKID_01uM",
                            "control_1uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_1uM")) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  filter(type %in% c("NMD")) %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  # filter(!type %in% c("other","lncRNA")) %>% 
  group_by(condition_2, type) %>% 
  summarize(median_log2FC=median(log2FC, na.rm = TRUE)) 

SMG89_SMG1i_KO_Lollipop_DTE_sig <- SMG89_edgeR_DTE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "SMG8_KO_0uM",
                            "SMG9_KO_0uM",
                            "SMG8_delKID_0uM",
                            "control_01uM",
                            "SMG8_KO_01uM",
                            "SMG9_KO_01uM",
                            "SMG8_delKID_01uM",
                            "control_1uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_1uM")) %>% 
  filter(type %in% c("NMD")) %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  # filter(!type %in% c("other","lncRNA")) %>%
  filter(qvalue < 0.0001 & abs(log2FC > 1)) %>% 
  mutate(up_down = case_when(log2FC > 1 ~ "up",
                             log2FC < -1 ~ "down")) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  group_by(condition_2, type) %>% 
  summarize(n=n()) %>% 
  left_join(SMG89_SMG1i_KO_Lollipop_DTE_median_log2FC) %>% 
  # mutate(n = case_when(type == "coding" ~ -n,
  # type == "NMD" ~ n)) %>% 
  ggplot(aes(y=condition_2,
             x=n)) +
  theme(legend.position="top", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  # geom_linerange(aes(xmin=log2FC_lower,
  #                    xmax=log2FC_upper),
  #                position = position_dodge(width = 0.9),
  #                color="black",
  #                linewidth=0.25) +
  geom_linerange(aes(xmin=0, xmax=n),
                 linewidth=0.25,
                 position = position_dodge(width=0.75)) +
  geom_point(aes(
    # size=abs(median_log2FC),
    # shape=case_when(median_log2FC > 0 ~ "positive",
    #                 median_log2FC < 0 ~ "negative"),
    fill=abs(median_log2FC)),
    shape=21,
    size=3,
    position = position_dodge(width=0.75),
    stroke=0.2) +
  # scale_shape_manual(values = c("positive" = 21,
  #                               "negative" = 23)) +
  # theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # scale_fill_manual(values=c("coding" = "#6A00A8FF",
  #                            "NMD" = "#60CEACFF")) +
  # scale_x_continuous(limits=c(-4,4)) +
  # scale_size(limits = c(0,1.5),
  #            range = c(0.25,3)) +
  scale_fill_gradient(low = ("white"),
                      high = ("#B2182B"),
                      na.value = "grey80") +
  labs(y="",
       x="Number of\nsig. regulated\nNMD-annotated transcripts",
       fill="Median\nDTE\nlog2FC",
       shape="Median log2FC",
       size="Median log2FC") +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F5_DESeq2_DGE_combined_forPlot$condition_2))*3+0.4, "mm"),
                   cols = unit(20, "mm")) +
  scale_x_continuous(expand = expansion(c(0.15, 0.15)),
                     breaks = c(0, 1500, 3000)) +
  # theme(axis.text.y=element_blank(),
  #       axis.ticks.y = element_blank()) +
  guides(size=guide_legend(nrow=2,byrow=TRUE),
         fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure5_DTE_Lollipop", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

##### Rev_1_F2_Examples -----------------------------------------------------------

# Define representative NMD-targeted transcripts
standard_heatmap_transcripts <- c("TMEM222-202",
                                  "GABARAPL1-208",
                                  "SRSF2-208",
                                  "SRSF2-204",
                                  "RPL3-204")

SMG89_Rev_1_F5_DESeq2_DTE_combined_forPlot %>%  
  filter(transcript_name %in% standard_heatmap_transcripts) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, transcript_name) %>% 
  mutate(transcript_name = fct_relevel(transcript_name,
                                       "TMEM222-202",
                                       "SRSF2-208",
                                       "SRSF2-204",
                                       "RPL3-204",
                                       "GABARAPL1-208"
  )) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=transcript_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="",
       fill="DTE log2FC") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                # barwidth = unit(8, "mm"),
                                # barheight = unit(2, "mm"),
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 0, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             keyheight = unit(2, "mm"),
                             label.position = "bottom")) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(unique(SMG89_Rev_1_F5_DESeq2_DTE_combined_forPlot$condition_2))*3+0.6, "mm"),
                   cols = unit(5*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure5_DTE_Dotplot_Examples", Sys.Date(), ".pdf")),
       width = cw2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

### Overlap of multiple Cores -----------------------------------------------


#### SMG1i 1µM core ----------------------------------------------------------


# Define common SMG1i gene set
SMG1i_DGE_core_set <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 %in% c("control_1uM",
                                   "SMG8_KO_1uM",
                                   "SMG9_KO_1uM")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_SMG1i = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up_SMG1i",
                                  log2FoldChange < -1 & padj < 0.0001 ~ "down_SMG1i",
                                  TRUE ~ "ns_SMG1i")) %>% 
  dplyr::select(gene_id, condition_2, UpDown_SMG1i) %>%
  group_by(condition_2) %>% 
  mutate(gene_id = as.character(gene_id)) %>% 
  group_by(gene_id) %>% 
  mutate(n_up_SMG1i = sum(UpDown_SMG1i == "up_SMG1i"),
         n_down_SMG1i = sum(UpDown_SMG1i == "down_SMG1i"),
         n_ns_SMG1i = sum(UpDown_SMG1i == "ns_SMG1i")) %>% 
  mutate(gene_set_SMG1i = case_when(n_up_SMG1i == 3 ~ "core_up_SMG1i",
                                    n_up_SMG1i == 2 ~ "shell_up_SMG1i",
                                    n_up_SMG1i == 1 ~ "cloud_up_SMG1i",
                                    n_down_SMG1i == 3 ~ "core_down_SMG1i",
                                    n_down_SMG1i == 2 ~ "shell_down_SMG1i",
                                    n_down_SMG1i == 1  ~ "cloud_down_SMG1i",
                                    TRUE ~ "ns_SMG1i")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_SMG1i") %>% 
  ungroup() %>% 
  mutate(gene_set_SMG1i = fct_relevel(gene_set_SMG1i,
                                      "core_up_SMG1i",
                                      "shell_up_SMG1i",
                                      "cloud_up_SMG1i",
                                      "ns_SMG1i",
                                      "cloud_down_SMG1i",
                                      "shell_down_SMG1i",
                                      "core_down_SMG1i"))

SMG1i_DGE_core_set %>% 
  ungroup() %>% 
  dplyr::count(gene_set_SMG1i)  


##### Plot --------------------------------------------------------------------

SMG1i_DGE_core_set %>% 
  ungroup() %>% 
  dplyr::count(gene_set_SMG1i) %>% 
  ggplot(aes(x="SMG1i core\nDGE",
             y=n,
             fill=gene_set_SMG1i)) +
  geom_col() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_blank(),
        # panel.grid.major.y = element_line(colour = 'darkgray', linewidth = 0.1),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = "black", linewidth = 0.1),
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        strip.text = element_text(size = 6),
        text=element_text(family="Arial")) +
  scale_fill_manual(values=c("core_up_SMG1i" = "#b2182b",
                             "shell_up_SMG1i" = "#d6604d",
                             "cloud_up_SMG1i" = "#f4a582",
                             "ns_SMG1i" = "darkgray",
                             "cloud_down_SMG1i" = "#92c5de",
                             "shell_down_SMG1i" = "#4393c3",
                             "core_down_SMG1i" = "#2166ac")) +
  geom_text(aes(label = n),
            color="black",
            size = 5*0.36,
            position = position_stack(vjust = 0.5))  +
  labs(x="",
       y="number of genes",
       fill="SMG1i\n1µM DGE") +
  force_panelsizes(rows = unit(40, "mm"),
                   cols = unit(10, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_SMG1i_core_definition.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### UpSet plot of cores -----------------------------------------------------

##### UPF1 core -------------------------------------------------------------
UPF1_diffDegrons_DGE_wide <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 %in% c("UPF1_Nter_12h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_FKBP_HEK_12h")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up",
                            log2FoldChange < -1 & padj < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  dplyr::select(gene_id, condition_2, UpDown) %>%
  group_by(condition_2) %>% 
  mutate(gene_id = as.character(gene_id)) %>% 
  group_by(gene_id) %>% 
  mutate(n_up = sum(UpDown == "up"),
         n_down = sum(UpDown == "down"),
         n_ns = sum(UpDown == "ns")) %>% 
  mutate(gene_set = case_when(n_up == 3 ~ "core_up",
                              n_up == 2 ~ "shell_up",
                              n_up == 1 ~ "cloud_up",
                              n_down == 3 ~ "core_down",
                              n_down == 2 ~ "shell_down",
                              n_down == 1 ~ "cloud_down",
                              TRUE ~ "ns")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown") %>% 
  ungroup() %>% 
  mutate(gene_set = fct_relevel(gene_set,
                                "core_up",
                                "shell_up",
                                "cloud_up",
                                "ns",
                                "cloud_down",
                                "shell_down",
                                "core_down"))

##### SMG567 core -------------------------------------------------------------

# Define common SMG567 gene set
SMG567_DGE_core_set <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 %in% c("SMG7_KO_2_SMG5_KD",
                                   "SMG7_KO_2_SMG6_KD",
                                   "SMG7_KO_34_SMG5_KD",
                                   "SMG7_KO_34_SMG6_KD")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_SMG567 = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up_SMG567",
                                   log2FoldChange < -1 & padj < 0.0001 ~ "down_SMG567",
                                   TRUE ~ "ns_SMG567")) %>% 
  dplyr::select(gene_id, condition_2, UpDown_SMG567) %>%
  group_by(condition_2) %>% 
  mutate(gene_id = as.character(gene_id)) %>% 
  group_by(gene_id) %>% 
  mutate(n_up_SMG567 = sum(UpDown_SMG567 == "up_SMG567"),
         n_down_SMG567 = sum(UpDown_SMG567 == "down_SMG567"),
         n_ns_SMG567 = sum(UpDown_SMG567 == "ns_SMG567")) %>% 
  mutate(gene_set_SMG567 = case_when(n_up_SMG567 == 4 ~ "core_up_SMG567",
                                     n_up_SMG567 == 3 ~ "shell_up_SMG567",
                                     n_up_SMG567 == 2 | n_up_SMG567 == 1 ~ "cloud_up_SMG567",
                                     n_down_SMG567 == 4 ~ "core_down_SMG567",
                                     n_down_SMG567 == 3 ~ "shell_down_SMG567",
                                     n_down_SMG567 == 2 | n_down_SMG567 == 1 ~ "cloud_down_SMG567",
                                     TRUE ~ "ns_SMG567")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_SMG567") %>% 
  ungroup() %>% 
  mutate(gene_set_SMG567 = fct_relevel(gene_set_SMG567,
                                       "core_up_SMG567",
                                       "shell_up_SMG567",
                                       "cloud_up_SMG567",
                                       "ns_SMG567",
                                       "cloud_down_SMG567",
                                       "shell_down_SMG567",
                                       "core_down_SMG567"))


##### UPF3 core ---------------------------------------------------------------

# Define common UPF3 gene set
UPF3_DGE_core_set <- SMG89_DESeq2_DGE_combined %>% 
  dplyr::filter(condition_2 %in% c("UPF3dKO1",
                                   "UPF3dKO2",
                                   "UPF3dKO1KD",
                                   "UPF3dKO2KD")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_UPF3 = case_when(log2FoldChange > 1 & padj < 0.0001 ~ "up_UPF3",
                                 log2FoldChange < -1 & padj < 0.0001 ~ "down_UPF3",
                                 TRUE ~ "ns_UPF3")) %>% 
  dplyr::select(gene_id, condition_2, UpDown_UPF3) %>%
  group_by(condition_2) %>% 
  mutate(gene_id = as.character(gene_id)) %>% 
  group_by(gene_id) %>% 
  mutate(n_up_UPF3 = sum(UpDown_UPF3 == "up_UPF3"),
         n_down_UPF3 = sum(UpDown_UPF3 == "down_UPF3"),
         n_ns_UPF3 = sum(UpDown_UPF3 == "ns_UPF3")) %>% 
  mutate(gene_set_UPF3 = case_when(n_up_UPF3 == 4 ~ "core_up_UPF3",
                                   n_up_UPF3 == 3 ~ "shell_up_UPF3",
                                   n_up_UPF3 == 2 | n_up_UPF3 == 1 ~ "cloud_up_UPF3",
                                   n_down_UPF3 == 4 ~ "core_down_UPF3",
                                   n_down_UPF3 == 3 ~ "shell_down_UPF3",
                                   n_down_UPF3 == 2 | n_down_UPF3 == 1 ~ "cloud_down_UPF3",
                                   TRUE ~ "ns_UPF3")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_UPF3") %>% 
  ungroup() %>% 
  mutate(gene_set_UPF3 = fct_relevel(gene_set_UPF3,
                                     "core_up_UPF3",
                                     "shell_up_UPF3",
                                     "cloud_up_UPF3",
                                     "ns_UPF3",
                                     "cloud_down_UPF3",
                                     "shell_down_UPF3",
                                     "core_down_UPF3"))

###### plot - UP ----------------------------------------------------------------------


# Filter for cutoffs and PTC == TRUE, get isoform id and condition2
all_SMG89_core_dat <- SMG1i_DGE_core_set %>% 
  dplyr::filter(gene_set_SMG1i == "core_up_SMG1i") %>% 
  dplyr::select(gene_id) %>% 
  mutate(condition="SMG1i core") %>% 
  bind_rows(UPF1_diffDegrons_DGE_wide %>% 
              dplyr::filter(gene_set == "core_up") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="UPF1 core")) %>% 
  bind_rows(SMG567_DGE_core_set %>% 
              dplyr::filter(gene_set_SMG567 == "core_up_SMG567") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="SMG567 core")) %>% 
  bind_rows(UPF3_DGE_core_set %>% 
              dplyr::filter(gene_set_UPF3 == "core_up_UPF3") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="UPF3 core")) %>% 
  mutate(gene_id = as.character(gene_id))

# Make lists
all_SMG89_core_UpSet_List <- split(all_SMG89_core_dat$gene_id,all_SMG89_core_dat$condition)

# Make combination matrix
Mat_Up = make_comb_mat(all_SMG89_core_UpSet_List, top_n_sets = 10)

Mat_Up <- Mat_Up[comb_degree(Mat_Up) >= 2]

# Filter for top 10 combinations
Mat_Up <- Mat_Up[order(comb_size(Mat_Up), decreasing = TRUE)[1:10]]


# Make UpSet plot for upregulated
# Make lists
all_SMG89_core_UpSet_ht_up <- UpSet(t(Mat_Up),
                                    name = "log2FC",
                                    column_names_gp = gpar(fontsize=6),
                                    pt_size = unit(3, "mm"),
                                    lwd = 2, comb_order = order(comb_size(Mat_Up)),
                                    # comb_col = brewer.pal(length(set_name(Mat_Up))+3, "Blues")[comb_degree(Mat_Up)+3],
                                    width = unit((length(set_name(Mat_Up))/2), "cm"),
                                    height = unit(2.5, "cm"),
                                    top_annotation = upset_top_annotation(t(Mat_Up), 
                                                                          add_numbers = TRUE,
                                                                          annotation_name_gp = gpar(fontsize = 6),
                                                                          axis_param = list(gp = gpar(fontsize = 6)),
                                                                          gp = gpar(fill = "black"),
                                                                          numbers_gp = gpar(fontsize = 6),
                                                                          numbers_rot = 90,
                                                                          height = unit(1, "cm")), 
                                    right_annotation = upset_right_annotation(t(Mat_Up),
                                                                              add_numbers = TRUE,
                                                                              annotation_name_gp = gpar(fontsize = 6),
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              gp = gpar(fill = "black"), 
                                                                              numbers_gp = gpar(fontsize = 6),
                                                                              numbers_rot = 0,
                                                                              width = unit(1, "cm"))
)


pdf(file = file.path("/home/volker/2023_SMG89_project", "Revision1_all_core_DGE_UpSet_up.pdf"))


draw(all_SMG89_core_UpSet_ht_up)


dev.off()

###### plot - DOWN ----------------------------------------------------------------------


# Filter for cutoffs and PTC == TRUE, get isoform id and condition2
all_SMG89_core_dat_down <- SMG1i_DGE_core_set %>% 
  dplyr::filter(gene_set_SMG1i == "core_down_SMG1i") %>% 
  dplyr::select(gene_id) %>% 
  mutate(condition="SMG1i core") %>% 
  bind_rows(UPF1_diffDegrons_DGE_wide %>% 
              dplyr::filter(gene_set == "core_down") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="UPF1 core")) %>% 
  bind_rows(SMG567_DGE_core_set %>% 
              dplyr::filter(gene_set_SMG567 == "core_down_SMG567") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="SMG567 core")) %>% 
  bind_rows(UPF3_DGE_core_set %>% 
              dplyr::filter(gene_set_UPF3 == "core_down_UPF3") %>% 
              dplyr::select(gene_id) %>% 
              mutate(condition="UPF3 core")) %>% 
  mutate(gene_id = as.character(gene_id))

# Make lists
all_SMG89_core_UpSet_List_down <- split(all_SMG89_core_dat_down$gene_id,all_SMG89_core_dat_down$condition)

# Make combination matrix
Mat_Down = make_comb_mat(all_SMG89_core_UpSet_List_down, top_n_sets = 10)

Mat_Down <- Mat_Down[comb_degree(Mat_Down) >= 2]

# Filter for top 10 combinations
Mat_Down <- Mat_Down[order(comb_size(Mat_Down), decreasing = TRUE)[1:10]]


# Make UpSet plot for upregulated
# Make lists
all_SMG89_core_UpSet_ht_down <- UpSet(t(Mat_Down),
                                      name = "log2FC",
                                      column_names_gp = gpar(fontsize=6),
                                      pt_size = unit(3, "mm"),
                                      lwd = 2, comb_order = order(comb_size(Mat_Down)),
                                      # comb_col = brewer.pal(length(set_name(Mat_Up))+3, "Blues")[comb_degree(Mat_Up)+3],
                                      width = unit((length(set_name(Mat_Down))/2), "cm"),
                                      height = unit(2.5, "cm"),
                                      top_annotation = upset_top_annotation(t(Mat_Down), 
                                                                            add_numbers = TRUE,
                                                                            annotation_name_gp = gpar(fontsize = 6),
                                                                            axis_param = list(gp = gpar(fontsize = 6)),
                                                                            gp = gpar(fill = "black"),
                                                                            numbers_gp = gpar(fontsize = 6),
                                                                            numbers_rot = 90,
                                                                            height = unit(1, "cm")), 
                                      right_annotation = upset_right_annotation(t(Mat_Down),
                                                                                add_numbers = TRUE,
                                                                                annotation_name_gp = gpar(fontsize = 6),
                                                                                axis_param = list(gp = gpar(fontsize = 6)),
                                                                                gp = gpar(fill = "black"), 
                                                                                numbers_gp = gpar(fontsize = 6),
                                                                                numbers_rot = 0,
                                                                                width = unit(1, "cm"))
)


pdf(file = file.path("/home/volker/2023_SMG89_project", "Revision1_all_all_core_DGE_UpSet_down.pdf"))


draw(all_SMG89_core_UpSet_ht_down)


dev.off()

####  NMD factor expression   -----------------------------------------------------------

# Define NMD factors and few representative NMD targets
A0_heatmap_genes <- c("UPF1",
                      "UPF2",
                      "UPF3A",
                      "UPF3B",
                      "SMG1",
                      "SMG5",
                      "SMG6",
                      "SMG7",
                      "SMG8",
                      "SMG9",
                      "EIF4A3",
                      "RBM8A",
                      "MAGOH",
                      "MAGOHB",
                      "CASC3")

SMG89_A0_DESeq2_DGE_selectedNMD_heatmap_length <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "control_01uM",
                            "control_1uM",
                            "SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  distinct(condition_2) %>% 
  pull(condition_2)

SMG89_A0_DESeq2_DGE_selectedNMD_heatmap <- SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_12h",
                            "control_01uM",
                            "control_1uM",
                            "SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM")) %>% 
  arrange(experimentSet) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  filter(gene_name %in% A0_heatmap_genes) %>% 
  dplyr::select(condition_2, experimentSet, log2FoldChange, padj, gene_name) %>% 
  mutate(gene_name = fct_relevel(gene_name,
                                 "UPF1",
                                 "UPF2",
                                 "UPF3A",
                                 "UPF3B",
                                 "SMG1",
                                 "SMG5",
                                 "SMG6",
                                 "SMG7",
                                 "SMG8",
                                 "SMG9",
                                 "EIF4A3",
                                 "RBM8A",
                                 "MAGOH",
                                 "MAGOHB",
                                 "CASC3",
                                 "ZFAS1",
                                 "SNHG12",
                                 "GAS5",
                                 "GADD45B")) %>% 
  mutate(padj = replace(padj, padj == 0, 1e-320)) %>% 
  ggplot(aes(x=gene_name,
             y=condition_2,
             fill=log2FoldChange
  )) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FoldChange,
    size=-log10(padj)),
    shape=21,
    stroke=0.1) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  theme_minimal() + 	
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="",
       y="") +
  coord_fixed(ratio=1) +
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 1,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(aspect.ratio = 1) +
  force_panelsizes(rows = unit(length(SMG89_A0_DESeq2_DGE_selectedNMD_heatmap_length)*3+0.6, "mm"),
                   cols = unit(15*3+0.6, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_AA0_NMD_overview_", Sys.Date(), ".pdf")),
       SMG89_A0_DESeq2_DGE_selectedNMD_heatmap,
       width = cw3,
       height = 40,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### delKID plus SMG1i -------------------------------------------------

SMG89_delKID_DGE_for_plot_sig_woNS <- SMG89_delKID_DGE_for_plot_sig %>% 
  filter(DGE_class != "n.s.")

SMG89_DESeq2_DGE_combined %>% 
  filter(condition_2 %in% c("SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM")) %>% 
  filter(gene_id %in% SMG89_delKID_DGE_for_plot_sig_woNS$gene_id) %>% 
  left_join(SMG89_delKID_DGE_for_plot_sig_woNS %>% dplyr::select(gene_id, DGE_class)) %>% 
  mutate(condition_2 = fct_rev(condition_2)) %>% 
  ggplot(aes(x=log2FoldChange,
             y=condition_2,
             fill=DGE_class)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = 'darkgray', linewidth = 0.1),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6), 
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        strip.text = element_text(size = 6),
        text=element_text(family="Arial")) +
  scale_fill_manual(values=c("up" = "#b2182b",
                             "down" = "#2166ac")) +
  scale_x_continuous(limits=c(-10,10)) +
  geom_boxplot(outlier.shape = NA,
               linewidth = 0.2) +
  labs(y="",
       fill="Sig. regulated\nin delKID-0µM") +
  force_panelsizes(rows = unit(20, "mm"),
                   cols = unit(20, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_delKID_SMG1i_up_down_boxplot_new.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### Expression 8/9KO 0.1 µM in Core -----------------------------------------

all_SMG89_core_set1 <- SMG1i_DGE_core_set %>% 
  dplyr::filter(gene_set_SMG1i == "core_up_SMG1i") %>% 
  pull(gene_id)

SMG89_DESeq2_DGE_combined %>% 
  filter(gene_id %in% all_SMG89_core_set1) %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF3dKO1",
                            "UPF1_Nter_0h",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_0h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_0h",
                            "UPF1_FKBP_HEK_12h",
                            "control_01uM",
                            "control_1uM",
                            "SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM"
  )) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                   "SMG7_KO_34_SMG6_KD",
                                   "UPF3dKO1",
                                   "UPF1_FKBP_HEK_0h",
                                   "UPF1_FKBP_HEK_12h",
                                   "UPF1_FKBP_HCT_0h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_Nter_0h",
                                   "UPF1_Nter_12h",
                                   "SMG8_KO_0uM",
                                   "SMG9_KO_0uM",
                                   "SMG8_delKID_0uM",
                                   "control_01uM",
                                   "SMG8_KO_01uM",
                                   "SMG9_KO_01uM",
                                   "SMG8_delKID_01uM",
                                   "control_1uM",
                                   "SMG8_KO_1uM",
                                   "SMG9_KO_1uM",
                                   "SMG8_delKID_1uM"))) %>% 
  ggplot(aes(x=log2FoldChange,
             y=condition_2,
             fill=condition_2)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6),
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="solid") +
  geom_vline(xintercept=1,
             linewidth=0.25,
             linetype="dashed") +
  stat_density_ridges(aes(height = after_stat(ndensity)),
                      scale = 0.9,
                      size=0.1,
                      # quantiles = c(0.25,0.5,0.75),
                      # calc_ecdf = TRUE,
                      geom = "density_ridges_gradient",
                      rel_min_height = 0.01) +
  geom_point(data=SMG89_DESeq2_DGE_combined %>% 
               filter(gene_name %in% c("ZFAS1", "GAS5")) %>% 
               filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                                         "UPF3dKO1",
                                         "UPF1_Nter_0h",
                                         "UPF1_Nter_12h",
                                         "UPF1_FKBP_HCT_0h",
                                         "UPF1_FKBP_HCT_12h",
                                         "UPF1_FKBP_HEK_0h",
                                         "UPF1_FKBP_HEK_12h",
                                         "control_01uM",
                                         "control_1uM",
                                         "SMG8_KO_0uM",
                                         "SMG8_KO_01uM",
                                         "SMG8_KO_1uM",
                                         "SMG9_KO_0uM",
                                         "SMG9_KO_01uM",
                                         "SMG9_KO_1uM",
                                         "SMG8_delKID_0uM",
                                         "SMG8_delKID_01uM",
                                         "SMG8_delKID_1uM"
               )) %>% 
               mutate(condition_2 = fct_rev(condition_2)),
             aes(color=gene_name),
             size=1) +
  scale_color_manual(values=c("ZFAS1" = "#3B8B79",
                              "GAS5" = "#F08143")) +
  scale_fill_manual(values=c("SMG7_KO_34_SMG6_KD" = "#5B5B5B",
                             "UPF3dKO1" = "#5B5B5B",
                             "UPF1_Nter_0h" = "#5B5B5B",
                             "UPF1_Nter_12h" = "#5B5B5B",
                             "UPF1_FKBP_HCT_0h" = "#5B5B5B",
                             "UPF1_FKBP_HCT_12h" = "#5B5B5B",
                             "UPF1_FKBP_HEK_0h" = "#5B5B5B",
                             "UPF1_FKBP_HEK_12h" = "#5B5B5B",
                             "control_01uM" = "#D9DADA",
                             "control_1uM"  = "#D9DADA",
                             "SMG8_KO_0uM" = "#084B94",
                             "SMG8_KO_01uM" = "#084B94",
                             "SMG8_KO_1uM" = "#084B94",
                             "SMG9_KO_0uM" = "#645A91",
                             "SMG9_KO_01uM" = "#645A91",
                             "SMG9_KO_1uM" = "#645A91",
                             "SMG8_delKID_0uM" = "#4599D1",
                             "SMG8_delKID_01uM" = "#4599D1",
                             "SMG8_delKID_1uM" = "#4599D1")) +
  scale_x_continuous(breaks = c(-5,-2.5,-1,0,1, 2.5, 5),
                     labels = c(-5,-2.5,-1,0,1, 2.5, 5), 
                     limits = c(-5,5)) +
  labs(y="") +
  guides(fill = "none") +
  force_panelsizes(rows = unit(19*3.5, "mm"),
                   cols = unit(40, "mm")) 

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_SMG1i_core_expressionLevels_allCond_ZFAS1_GAS5.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

#### DTE cor SMG1i -----------------------------------------
SMG1i_DTE_core_set <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 %in% c("control_1uM",
                                   "SMG8_KO_1uM",
                                   "SMG9_KO_1uM")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  dplyr::rename("log2FC" = "logFC",
                "qvalue" = "FDR") %>% 
  mutate(UpDown_SMG1i = case_when(log2FC > 1 & qvalue < 0.0001 ~ "up_SMG1i",
                                  log2FC < -1 & qvalue < 0.0001 ~ "down_SMG1i",
                                  TRUE ~ "ns_SMG1i")) %>% 
  dplyr::select(transcript_id, condition_2, UpDown_SMG1i) %>%
  group_by(condition_2) %>% 
  mutate(transcript_id = as.character(transcript_id)) %>% 
  group_by(transcript_id) %>% 
  mutate(n_up_SMG1i = sum(UpDown_SMG1i == "up_SMG1i"),
         n_down_SMG1i = sum(UpDown_SMG1i == "down_SMG1i"),
         n_ns_SMG1i = sum(UpDown_SMG1i == "ns_SMG1i")) %>% 
  mutate(transcript_set_SMG1i = case_when(n_up_SMG1i == 3 ~ "core_up_SMG1i",
                                          n_up_SMG1i == 2 ~ "shell_up_SMG1i",
                                          n_up_SMG1i == 1 ~ "cloud_up_SMG1i",
                                          n_down_SMG1i == 3 ~ "core_down_SMG1i",
                                          n_down_SMG1i == 2 ~ "shell_down_SMG1i",
                                          n_down_SMG1i == 1  ~ "cloud_down_SMG1i",
                                          TRUE ~ "ns_SMG1i")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_SMG1i") %>% 
  ungroup() %>% 
  mutate(transcript_set_SMG1i = fct_relevel(transcript_set_SMG1i,
                                            "core_up_SMG1i",
                                            "shell_up_SMG1i",
                                            "cloud_up_SMG1i",
                                            "ns_SMG1i",
                                            "cloud_down_SMG1i",
                                            "shell_down_SMG1i",
                                            "core_down_SMG1i"))

SMG1i_DTE_core_set %>% 
  ungroup() %>% 
  dplyr::count(transcript_set_SMG1i) 

#### DTE Expression 8/9KO 0.1 µM in Core -----------------------------------------

DTE_all_SMG89_core_set1 <- SMG1i_DTE_core_set %>% 
  dplyr::filter(transcript_set_SMG1i == "core_up_SMG1i") %>% 
  pull(transcript_id)

SMG89_edgeR_DTE_combined %>% 
  filter(transcript_id %in% DTE_all_SMG89_core_set1) %>% 
  filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                            "UPF3dKO1",
                            "UPF1_Nter_0h",
                            "UPF1_Nter_12h",
                            "UPF1_FKBP_HCT_0h",
                            "UPF1_FKBP_HCT_12h",
                            "UPF1_FKBP_HEK_0h",
                            "UPF1_FKBP_HEK_12h",
                            "control_01uM",
                            "control_1uM",
                            "SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM",
                            "SMG8_delKID_0uM",
                            "SMG8_delKID_01uM",
                            "SMG8_delKID_1uM"
  )) %>% 
  mutate(condition_2 = fct_rev(fct_relevel(condition_2,
                                           "SMG7_KO_34_SMG6_KD",
                                           "UPF3dKO1",
                                           "UPF1_FKBP_HEK_0h",
                                           "UPF1_FKBP_HEK_12h",
                                           "UPF1_FKBP_HCT_0h",
                                           "UPF1_FKBP_HCT_12h",
                                           "UPF1_Nter_0h",
                                           "UPF1_Nter_12h",
                                           "SMG8_KO_0uM",
                                           "SMG9_KO_0uM",
                                           "SMG8_delKID_0uM",
                                           "control_01uM",
                                           "SMG8_KO_01uM",
                                           "SMG9_KO_01uM",
                                           "SMG8_delKID_01uM",
                                           "control_1uM",
                                           "SMG8_KO_1uM",
                                           "SMG9_KO_1uM",
                                           "SMG8_delKID_1uM"))) %>% 
  ggplot(aes(x=logFC,
             y=condition_2,
             fill=condition_2)) +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.y = element_line(color = "black",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_text(size=6),
        axis.text=element_text(size=6), 
        axis.title=element_text(size=6),
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial")) +
  geom_vline(xintercept=0,
             linewidth=0.25,
             linetype="solid") +
  geom_vline(xintercept=1,
             linewidth=0.25,
             linetype="dashed") +
  stat_density_ridges(aes(height = after_stat(ndensity)),
                      scale = 0.9,
                      size=0.1,
                      # quantiles = c(0.25,0.5,0.75),
                      # calc_ecdf = TRUE,
                      geom = "density_ridges_gradient",
                      rel_min_height = 0.01) +
  geom_point(data=SMG89_edgeR_DTE_combined %>% 
               filter(transcript_name %in% c("SRSF2-202",
                                             "SRSF2-204",
                                             "SRSF2-208")) %>% 
               filter(condition_2 %in% c("SMG7_KO_34_SMG6_KD",
                                         "UPF3dKO1",
                                         "UPF1_Nter_0h",
                                         "UPF1_Nter_12h",
                                         "UPF1_FKBP_HCT_0h",
                                         "UPF1_FKBP_HCT_12h",
                                         "UPF1_FKBP_HEK_0h",
                                         "UPF1_FKBP_HEK_12h",
                                         "control_01uM",
                                         "control_1uM",
                                         "SMG8_KO_0uM",
                                         "SMG8_KO_01uM",
                                         "SMG8_KO_1uM",
                                         "SMG9_KO_0uM",
                                         "SMG9_KO_01uM",
                                         "SMG9_KO_1uM",
                                         "SMG8_delKID_0uM",
                                         "SMG8_delKID_01uM",
                                         "SMG8_delKID_1uM")) %>% 
               mutate(condition_2 = fct_rev(condition_2)),
             aes(color=transcript_name),
             size=1) +
  scale_color_manual(values=c("SRSF2-204" = "#F08143",
                              "SRSF2-208" = "#FFED00",
                              "SRSF2-202" = "#3B8B79")) +
  scale_fill_manual(values=c("SMG7_KO_34_SMG6_KD" = "#5B5B5B",
                             "UPF3dKO1" = "#5B5B5B",
                             "UPF1_Nter_0h" = "#5B5B5B",
                             "UPF1_Nter_12h" = "#5B5B5B",
                             "UPF1_FKBP_HCT_0h" = "#5B5B5B",
                             "UPF1_FKBP_HCT_12h" = "#5B5B5B",
                             "UPF1_FKBP_HEK_0h" = "#5B5B5B",
                             "UPF1_FKBP_HEK_12h" = "#5B5B5B",
                             "control_01uM" = "#D9DADA",
                             "control_1uM"  = "#D9DADA",
                             "SMG8_KO_0uM" = "#084B94",
                             "SMG8_KO_01uM" = "#084B94",
                             "SMG8_KO_1uM" = "#084B94",
                             "SMG9_KO_0uM" = "#645A91",
                             "SMG9_KO_01uM" = "#645A91",
                             "SMG9_KO_1uM" = "#645A91",
                             "SMG8_delKID_0uM" = "#4599D1",
                             "SMG8_delKID_01uM" = "#4599D1",
                             "SMG8_delKID_1uM" = "#4599D1")) +
  scale_x_continuous(breaks = c(-10,0,1, 5, 10),
                     labels = c(-10,0,1, 5,  10), 
                     limits = c(-15,15)) +
  labs(y="") +
  guides(fill = "none") +
  force_panelsizes(rows = unit(19*3.5, "mm"),
                   cols = unit(40, "mm")) 

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_SMG1i_core_DTE_expressionLevels_allCond_SRSF2.pdf"),
       width = cw3-2,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


## Other analyses -----------------------------------------------------


### CPM_DTE -----------------------------------------------------------------

SMG89_UPF1_AID_DTE_mydir="/home/volker/gencode.v42.datasets/2023_UPF1_AID_DW"
SMG89_DTE_mydir="/home/volker/gencode.v42.datasets/2023_11_SMG89_SK"

# Get samples and check if all files are present
SMG89_UPF1_AID_DTE_samples <- read.table(file.path(SMG89_UPF1_AID_DTE_mydir, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "control_0h",
                               TRUE ~ condition))

SMG89_DTE_samples <- read.table(file.path(SMG89_DTE_mydir, "Samples.txt"), header = TRUE) %>% 
  mutate(sample = as.character(sample)) %>% 
  mutate(condition = case_when(condition == "control" ~ "control_0uM",
                               TRUE ~ condition))

CPM_DTE_samples <- as_tibble(bind_rows(SMG89_UPF1_AID_DTE_samples,
                                       SMG89_DTE_samples))

# Get unique conditions
CPM_DTE_condition <- CPM_DTE_samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate data frame used for tximeta
CPM_DTE_coldata <- tibble(files = file.path(SMG89_UPF1_AID_DTE_mydir, "Salmon", SMG89_UPF1_AID_DTE_samples$sample)) %>% 
  bind_rows(tibble(files = file.path(SMG89_DTE_mydir, "Salmon", SMG89_DTE_samples$sample))) 

# Supplement with sample IDs as "names"
CPM_DTE_coldata$names <- CPM_DTE_samples$sample

# Join with samples to get condition
CPM_DTE_coldata <- CPM_DTE_coldata %>% 
  left_join(CPM_DTE_samples,
            by = c("names" = "sample"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(CPM_DTE_coldata$files)))

# Import transcript-level counts 
catch <- catchSalmon(paths = CPM_DTE_coldata$files)

# Account for the mapping ambiguity
scaled.counts <- catch$counts/catch$annotation$Overdispersion

# Create DGEList objec
DGEList <- DGEList(counts = scaled.counts,
                   samples = CPM_DTE_samples,
                   group = CPM_DTE_samples$condition,
                   genes = catch$annotation)

CPM_DTE_cpm <- as_tibble(cpmByGroup(DGEList),rownames="transcript_id") %>% 
  left_join(gtf_gencode_df_short %>% filter(type=="transcript") %>% dplyr::select(transcript_id, transcript_name))

CPM_DTE_cpm_long <- CPM_DTE_cpm %>% 
  pivot_longer(cols=-c(transcript_id, transcript_name),
               values_to = "cpm",
               names_to = "condition_2")

CPM_DTE_cpm_long %>% 
  filter(transcript_name %in% c("SRSF2-204", "SRSF2-208", "ZFAS1-206", "ZFAS1-204", "ZFAS1-205", "GAS5-245", "GAS5-255", "GAS5-259")) %>% 
  mutate(transcript_name = fct_relevel(transcript_name,
                                       "ZFAS1-204", "ZFAS1-205", "ZFAS1-206", "GAS5-245", "GAS5-255", "GAS5-259")) %>% 
  mutate(condition_2 = fct_relevel(condition_2,
                                   "control_0uM",
                                   "control_01uM",
                                   "control_1uM",
                                   "SMG8_KO_0uM",
                                   "SMG8_KO_01uM",
                                   "SMG8_KO_1uM",
                                   "SMG9_KO_0uM",
                                   "SMG9_KO_01uM",
                                   "SMG9_KO_1uM",
                                   "control_0h",
                                   "control_48h",
                                   "UPF1_Nter_0h", 
                                   "UPF1_Nter_2h", 
                                   "UPF1_Nter_4h", 
                                   "UPF1_Nter_8h",
                                   "UPF1_Nter_12h",
                                   "UPF1_Nter_24h", 
                                   "UPF1_Nter_48h")) %>% 
    mutate(concentration = case_when(condition_2 %in% c("control_0uM","SMG8_KO_0uM","SMG9_KO_0uM") ~ "0 µM SMG1i",
                                     condition_2 %in% c("control_01uM","SMG8_KO_01uM","SMG9_KO_01uM") ~ "0.1 µM SMG1i",
                                     condition_2 %in% c("control_1uM","SMG8_KO_1uM","SMG9_KO_1uM") ~ "1 µM SMG1i",
                                     TRUE ~ "500 µM IAA")) %>% 
  mutate(timepoint = case_when(condition_2 == "control_0h" ~ -12,
                               condition_2 == "control_48h" ~ -6,
                               condition_2 == "UPF1_Nter_0h" ~ 0,
                               condition_2 == "UPF1_Nter_2h" ~ 2,
                               condition_2 == "UPF1_Nter_4h" ~ 4,
                               condition_2 == "UPF1_Nter_8h" ~ 8,
                               condition_2 == "UPF1_Nter_12h" ~ 12,
                               condition_2 == "UPF1_Nter_24h" ~ 24,
                               condition_2 == "UPF1_Nter_48h" ~ 48,
                               condition_2 %in% c("control_0uM","control_01uM","control_1uM") ~ 54,
                               condition_2 %in% c("SMG8_KO_0uM","SMG8_KO_01uM","SMG8_KO_1uM") ~ 58,
                               condition_2 %in% c("SMG9_KO_0uM","SMG9_KO_01uM","SMG9_KO_1uM") ~ 62)) %>% 
  mutate(group=case_when(condition_2 %in% c("control_0uM","control_01uM","control_1uM") ~ "control",
                         condition_2 %in% c("SMG8_KO_0uM","SMG8_KO_01uM","SMG8_KO_1uM") ~ "SMG8 KO",
                         condition_2 %in% c("SMG9_KO_0uM","SMG9_KO_01uM","SMG9_KO_1uM") ~ "SMG9 KO",
                         condition_2 %in% c("UPF1_Nter_0h","UPF1_Nter_2h","UPF1_Nter_4h","UPF1_Nter_8h",
                                            "UPF1_Nter_12h",
                                            "UPF1_Nter_24h", 
                                            "UPF1_Nter_48h") ~ "N-AID-UPF1",
                         condition_2 %in% c("control_0h","control_48h") ~ "control_UPF1")) %>% 
  filter(!condition_2 %in% c("control_0h","control_48h")) %>% 
    ggplot(aes(x=timepoint,
               y=cpm)) +
  geom_line(aes(group=group,
                color=group),
            alpha=0.5,
            size=0.5) +
  geom_point(aes(fill=concentration),
             shape=21,
             color="black",
             stroke=0.2) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(colour = 'darkgray', linewidth = 0.1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("N-AID-UPF1" = "#DACB47",
                              "control" = "darkgray",
                              "SMG8 KO" = "#084B94",
                              "SMG9 KO" = "#645A91")) +
  scale_fill_manual(values=c("500 µM IAA" = "#DACB47",
                              "0 µM SMG1i" = "gray90",
                              "0.1 µM SMG1i" = "gray50",
                              "1 µM SMG1i" = "gray10")) +
  scale_x_continuous(breaks = c(0,12,24,36,48,54,58,62),
                     labels = c("0",
                                "12",
                                "24",
                                "36",
                                "48",
                                "control",
                                "SMG8 KO",
                                "SMG9 KO")) +
  facet_wrap(~transcript_name,
             scales = "free_y") +
  force_panelsizes(total_width = unit(160, "mm"),
                   total_height = unit(100, "mm"))

ggsave(file.path("/home/volker/2023_SMG89_project", "Revision1_SRSF2_ZFAS1_GAS5.pdf"),
       width = 30,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")


#### Top /DTE targets -----------------------------------------------------

# Top 10 SMG8 and SMG9 0.1 µM NMD-annotated DTE targets
SMG8_01uM_DTE_NMD_top10_transcripts_FDR <- SMG89_edgeR_DTE_combined %>% 
  filter(type == "NMD") %>% 
  filter(condition_2 %in% c("SMG8_KO_01uM")) %>% 
  arrange(FDR) %>% 
  slice_min(FDR, n=10) %>% 
  dplyr::select(transcript_name) %>% 
  mutate(reason="top_FDR")

SMG8_01uM_DTE_NMD_top10_transcripts_logFC <- SMG89_edgeR_DTE_combined %>% 
  filter(type == "NMD") %>% 
  filter(condition_2 %in% c("SMG8_KO_01uM")) %>% 
  arrange(desc(logFC)) %>% 
  slice_max(logFC, n=10) %>% 
  dplyr::select(transcript_name) %>% 
  mutate(reason="top_logFC")

SMG9_01uM_DTE_NMD_top10_transcripts_FDR <- SMG89_edgeR_DTE_combined %>% 
  filter(type == "NMD") %>% 
  filter(condition_2 %in% c("SMG9_KO_01uM")) %>% 
  arrange(FDR) %>% 
  slice_min(FDR, n=10) %>% 
  dplyr::select(transcript_name) %>% 
  mutate(reason="top_FDR")

SMG9_01uM_DTE_NMD_top10_transcripts_logFC <- SMG89_edgeR_DTE_combined %>% 
  filter(type == "NMD") %>% 
  filter(condition_2 %in% c("SMG9_KO_01uM")) %>% 
  arrange(desc(logFC)) %>% 
  slice_max(logFC, n=10) %>% 
  dplyr::select(transcript_name) %>% 
  mutate(reason="top_logFC")

Other_top_transcripts <- tibble(transcript_name = c("SRSF2-204", "SRSF2-208", "ZFAS1-206", "ZFAS1-204", "ZFAS1-205"),
                                reason="standard")

SMG8_SMG9_01uM_DTE_NMD_top_transcripts <- bind_rows(SMG8_01uM_DTE_NMD_top10_transcripts_FDR, SMG8_01uM_DTE_NMD_top10_transcripts_logFC,
                                                    SMG9_01uM_DTE_NMD_top10_transcripts_FDR, SMG9_01uM_DTE_NMD_top10_transcripts_logFC, Other_top_transcripts) %>% 
  distinct(transcript_name, .keep_all = TRUE)

# DTE
SMG89_edgeR_DTE_combined %>% 
  filter(transcript_name %in% SMG8_SMG9_01uM_DTE_NMD_top_transcripts$transcript_name) %>%
  filter(condition_2 %in% c("SMG8_KO_0uM",
                            "SMG8_KO_01uM",
                            "SMG8_KO_1uM",
                            "SMG9_KO_0uM",
                            "SMG9_KO_01uM",
                            "SMG9_KO_1uM"
                            )) %>% 
  mutate(cell_line = case_when(condition_2 %in% c("SMG8_KO_0uM",
                                                  "SMG8_KO_01uM",
                                                  "SMG8_KO_1uM") ~ "SMG8",
                               condition_2 %in% c("SMG9_KO_0uM",
                                                  "SMG9_KO_01uM",
                                                  "SMG9_KO_1uM") ~ "SMG9")) %>% 
  left_join(SMG8_SMG9_01uM_DTE_NMD_top_transcripts) %>% 
  arrange(reason, desc(logFC)) %>% 
  mutate(transcript_name = fct_rev(fct_inorder(transcript_name))) %>%
  # mutate(transcript_name = (fct_reorder2(transcript_name,
  #                                       reason,
  #                                       logFC))) %>% 
  # mutate(transcript_name = fct_rev(fct_relevel(transcript_name,
  #                                              "ZFAS1-206", "ZFAS1-205", "ZFAS1-204", "SRSF2-208", "SRSF2-204"))) %>%
  ggplot(aes(x=logFC,
             y=transcript_name,
             fill=condition_2,
             size=-log10(FDR))) +
  geom_linerange(aes(xmin=0, xmax=logFC, color=condition_2),
                 linewidth=0.25,
                 alpha=0.5,
                 show_guide = FALSE) +
  geom_point(shape=21,
             stroke=0.1,
             alpha=0.75
             ) +
  scale_size(range = c(1, 3)) +
  theme_minimal() +
  theme(legend.position="right", 
        strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid.major.x = element_line(color = "gray80",
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(colour = 'black',size=6), 
        axis.title=element_text(colour = 'black',size=6), 
        axis.line.x = element_line(colour = 'black', linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        legend.title = element_text(colour = 'black',size = 6), 
        legend.text = element_text(colour = 'black',size = 6), 
        legend.key.size = unit(0.25, "cm"),
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        legend.box="vertical", 
        plot.title = element_text(size = 6), 
        plot.subtitle = element_text(size = 6), 
        plot.caption = element_text(size = 6),
        text=element_text(family="Arial", color="black")) +
  geom_vline(xintercept = 0,
             color="darkgray",
             linewidth=0.25) +
  scale_color_manual(values=c("SMG8_KO_0uM" = "gray70",
                             "SMG8_KO_01uM" = "#6491CA",
                             "SMG8_KO_1uM" = "#063363",
                             "SMG9_KO_0uM" = "gray70",
                             "SMG9_KO_01uM" = "#9088B9",
                             "SMG9_KO_1uM" = "#372E5D")) +
  scale_fill_manual(values=c("SMG8_KO_0uM" = "gray70",
                             "SMG8_KO_01uM" = "#6491CA",
                             "SMG8_KO_1uM" = "#063363",
                             "SMG9_KO_0uM" = "gray70",
                             "SMG9_KO_01uM" = "#9088B9",
                             "SMG9_KO_1uM" = "#372E5D")) +
  facet_wrap(~cell_line,ncol=2) +
  force_panelsizes(rows = unit(length(unique(SMG8_SMG9_01uM_DTE_NMD_top_transcripts$transcript_name))*2+0.4, "mm"),
                   cols = unit(60, "mm")) +
  labs(y="",
       x="DTE log2FC",
       fill="",
       size="-log10(FDR)") 

ggsave(file.path("/home/volker/2023_SMG89_project", paste0("Revision1_SMG89_Figure5_DTE_Top_transcripts", Sys.Date(), ".pdf")),
       width = cw3,
       height = 20,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")
  

## Tables EV ------------------------------------------------------------------


### Table EV2 - DGE/DTE ----------------------------------------------

# DGE
SMG89_DESeq2_DGE_combined %>% 
  filter(padj < 0.0001 & experimentSet == "HCT116_SMG89KO_SMG1i_this_Study") %>%  
  dplyr::select(-c(experimentSet,publicationName)) %>% 
  dplyr::rename("GENCODE_gene_biotype" = "gene_type",
                "gene_biotype_short" = "type") %>% 
  relocate(c(gene_name, GENCODE_gene_biotype, gene_biotype_short), .after = gene_id) %>% 
  write_excel_csv(file.path("/home/volker/2023_SMG89_project", "Revision1_Table_S2_part1.csv"))

# DTE
SMG89_edgeR_DTE_combined %>% 
  filter(FDR < 0.0001 & experimentSet == "HCT116_SMG89KO_SMG1i_this_Study") %>%  
  dplyr::rename("GENCODE_gene_biotype" = "gene_type",
                "GENCODE_transcript_biotype" = "transcript_type",
                "transcript_biotype_short" = "type") %>% 
  relocate(c(transcript_biotype_short), .after = GENCODE_transcript_biotype) %>% 
  relocate(c(transcript_name), .after = transcript_id) %>% 
  write_excel_csv(file.path("/home/volker/2023_SMG89_project", "Revision1_Table_S2_part2.csv"))

### Table EV3 - core/shell/cloud ----------------------------------------------

# DGE cores
DGE_all_cores_table <- SMG1i_DGE_core_set %>% 
  dplyr::select(gene_id,gene_set_SMG1i) %>% 
  full_join(UPF1_diffDegrons_DGE_wide %>% 
              dplyr::select(gene_id,gene_set) %>% 
              mutate(gene_set = paste0(gene_set,"_UPF1")) %>% 
              dplyr::rename("gene_set_UPF1" = "gene_set")) %>% 
  full_join(SMG567_DGE_core_set %>% 
              dplyr::select(gene_id,gene_set_SMG567)) %>% 
  full_join(UPF3_DGE_core_set %>% 
              dplyr::select(gene_id,gene_set_UPF3)) %>% 
  left_join(gtf_gencode_df_short %>% 
              filter(type=="gene") %>%
              dplyr::select(gene_id, gene_name, gene_type)
  ) %>% 
  relocate(c(gene_name,gene_type), .after = gene_id) 

DGE_all_cores_table %>% 
  write_excel_csv(file.path("/home/volker/2023_SMG89_project", "Revision1_Table_S3_part1.csv"))

# DTE cores
SMG1i_DTE_core_set <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 %in% c("control_1uM",
                                   "SMG8_KO_1uM",
                                   "SMG9_KO_1uM")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_SMG1i = case_when(logFC > 1 & FDR < 0.0001 ~ "up_SMG1i",
                                  logFC < -1 & FDR < 0.0001 ~ "down_SMG1i",
                                  TRUE ~ "ns_SMG1i")) %>% 
  dplyr::select(transcript_id, condition_2, UpDown_SMG1i) %>%
  group_by(condition_2) %>% 
  mutate(transcript_id = as.character(transcript_id)) %>% 
  group_by(transcript_id) %>% 
  mutate(n_up_SMG1i = sum(UpDown_SMG1i == "up_SMG1i"),
         n_down_SMG1i = sum(UpDown_SMG1i == "down_SMG1i"),
         n_ns_SMG1i = sum(UpDown_SMG1i == "ns_SMG1i")) %>% 
  mutate(transcript_set_SMG1i = case_when(n_up_SMG1i == 3 ~ "core_up_SMG1i",
                                          n_up_SMG1i == 2 ~ "shell_up_SMG1i",
                                          n_up_SMG1i == 1 ~ "cloud_up_SMG1i",
                                          n_down_SMG1i == 3 ~ "core_down_SMG1i",
                                          n_down_SMG1i == 2 ~ "shell_down_SMG1i",
                                          n_down_SMG1i == 1  ~ "cloud_down_SMG1i",
                                          TRUE ~ "ns_SMG1i")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_SMG1i") %>% 
  ungroup() %>% 
  mutate(transcript_set_SMG1i = fct_relevel(transcript_set_SMG1i,
                                            "core_up_SMG1i",
                                            "shell_up_SMG1i",
                                            "cloud_up_SMG1i",
                                            "ns_SMG1i",
                                            "cloud_down_SMG1i",
                                            "shell_down_SMG1i",
                                            "core_down_SMG1i"))

SMG1i_DTE_core_set %>% 
  ungroup() %>% 
  dplyr::count(transcript_set_SMG1i)  
#### Core Definition DTE UPF1 -----------------------------------------------------

UPF1_diffDegrons_DTE_wide <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 %in% c("UPF1_Nter_12h",
                                   "UPF1_FKBP_HCT_12h",
                                   "UPF1_FKBP_HEK_12h")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown = case_when(logFC > 1 & FDR < 0.0001 ~ "up",
                            logFC < -1 & FDR < 0.0001 ~ "down",
                            TRUE ~ "ns")) %>% 
  dplyr::select(transcript_id, condition_2, UpDown) %>%
  group_by(condition_2) %>% 
  mutate(transcript_id = as.character(transcript_id)) %>% 
  group_by(transcript_id) %>% 
  mutate(n_up = sum(UpDown == "up"),
         n_down = sum(UpDown == "down"),
         n_ns = sum(UpDown == "ns")) %>% 
  mutate(transcript_set = case_when(n_up == 3 ~ "core_up",
                                    n_up == 2 ~ "shell_up",
                                    n_up == 1 ~ "cloud_up",
                                    n_down == 3 ~ "core_down",
                                    n_down == 2 ~ "shell_down",
                                    n_down == 1 ~ "cloud_down",
                                    TRUE ~ "ns")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown") %>% 
  ungroup() %>% 
  mutate(transcript_set = fct_relevel(transcript_set,
                                      "core_up",
                                      "shell_up",
                                      "cloud_up",
                                      "ns",
                                      "cloud_down",
                                      "shell_down",
                                      "core_down"))

# Count gene_set for global overview
UPF1_diffDegrons_DTE_wide %>% 
  ungroup() %>% 
  dplyr::count(transcript_set)

#### Core Definition DTE SMG567  -----------------------------------------------------
SMG567_DTE_core_set <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 %in% c("SMG7_KO_2_SMG5_KD",
                                   "SMG7_KO_2_SMG6_KD",
                                   "SMG7_KO_34_SMG5_KD",
                                   "SMG7_KO_34_SMG6_KD")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_SMG567 = case_when(logFC > 1 & FDR < 0.0001 ~ "up_SMG567",
                                   logFC < -1 & FDR < 0.0001 ~ "down_SMG567",
                                   TRUE ~ "ns_SMG567")) %>% 
  dplyr::select(transcript_id, condition_2, UpDown_SMG567) %>%
  group_by(condition_2) %>% 
  mutate(transcript_id = as.character(transcript_id)) %>% 
  group_by(transcript_id) %>% 
  mutate(n_up_SMG567 = sum(UpDown_SMG567 == "up_SMG567"),
         n_down_SMG567 = sum(UpDown_SMG567 == "down_SMG567"),
         n_ns_SMG567 = sum(UpDown_SMG567 == "ns_SMG567")) %>% 
  mutate(transcript_set_SMG567 = case_when(n_up_SMG567 == 4 ~ "core_up_SMG567",
                                           n_up_SMG567 == 3 ~ "shell_up_SMG567",
                                           n_up_SMG567 == 2 | n_up_SMG567 == 1 ~ "cloud_up_SMG567",
                                           n_down_SMG567 == 4 ~ "core_down_SMG567",
                                           n_down_SMG567 == 3 ~ "shell_down_SMG567",
                                           n_down_SMG567 == 2 | n_down_SMG567 == 1 ~ "cloud_down_SMG567",
                                           TRUE ~ "ns_SMG567")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_SMG567") %>% 
  ungroup() %>% 
  mutate(transcript_set_SMG567 = fct_relevel(transcript_set_SMG567,
                                             "core_up_SMG567",
                                             "shell_up_SMG567",
                                             "cloud_up_SMG567",
                                             "ns_SMG567",
                                             "cloud_down_SMG567",
                                             "shell_down_SMG567",
                                             "core_down_SMG567"))

SMG567_DTE_core_set %>% 
  ungroup() %>% 
  dplyr::count(transcript_set_SMG567)

#### Core Definition DTE UPF3  -----------------------------------------------------
UPF3_DTE_core_set <- SMG89_edgeR_DTE_combined %>% 
  dplyr::filter(condition_2 %in% c("UPF3dKO1",
                                   "UPF3dKO2",
                                   "UPF3dKO1KD",
                                   "UPF3dKO2KD")) %>% 
  mutate(condition_2 = fct_drop(condition_2)) %>%
  mutate(UpDown_UPF3 = case_when(logFC > 1 & FDR < 0.0001 ~ "up_UPF3",
                                 logFC < -1 & FDR < 0.0001 ~ "down_UPF3",
                                 TRUE ~ "ns_UPF3")) %>% 
  dplyr::select(transcript_id, condition_2, UpDown_UPF3) %>%
  group_by(condition_2) %>% 
  mutate(transcript_id = as.character(transcript_id)) %>% 
  group_by(transcript_id) %>% 
  mutate(n_up_UPF3 = sum(UpDown_UPF3 == "up_UPF3"),
         n_down_UPF3 = sum(UpDown_UPF3 == "down_UPF3"),
         n_ns_UPF3 = sum(UpDown_UPF3 == "ns_UPF3")) %>% 
  mutate(transcript_set_UPF3 = case_when(n_up_UPF3 == 4 ~ "core_up_UPF3",
                                         n_up_UPF3 == 3 ~ "shell_up_UPF3",
                                         n_up_UPF3 == 2 | n_up_UPF3 == 1 ~ "cloud_up_UPF3",
                                         n_down_UPF3 == 4 ~ "core_down_UPF3",
                                         n_down_UPF3 == 3 ~ "shell_down_UPF3",
                                         n_down_UPF3 == 2 | n_down_UPF3 == 1 ~ "cloud_down_UPF3",
                                         TRUE ~ "ns_UPF3")) %>% 
  pivot_wider(names_from = "condition_2",
              values_from = "UpDown_UPF3") %>% 
  ungroup() %>% 
  mutate(transcript_set_UPF3 = fct_relevel(transcript_set_UPF3,
                                           "core_up_UPF3",
                                           "shell_up_UPF3",
                                           "cloud_up_UPF3",
                                           "ns_UPF3",
                                           "cloud_down_UPF3",
                                           "shell_down_UPF3",
                                           "core_down_UPF3"))

UPF3_DTE_core_set %>% 
  ungroup() %>% 
  dplyr::count(transcript_set_UPF3)


DTE_all_cores_table <- SMG1i_DTE_core_set %>% 
  dplyr::select(transcript_id, transcript_set_SMG1i) %>% 
  full_join(UPF1_diffDegrons_DTE_wide %>% 
              dplyr::select(transcript_id,transcript_set) %>% 
              mutate(transcript_set = paste0(transcript_set,"_UPF1")) %>% 
              dplyr::rename("transcript_set_UPF1" = "transcript_set")) %>% 
  full_join(SMG567_DTE_core_set %>% 
              dplyr::select(transcript_id, transcript_set_SMG567)) %>% 
  full_join(UPF3_DTE_core_set %>% 
              dplyr::select(transcript_id, transcript_set_UPF3)) %>% 
  left_join(gtf_gencode_df_short %>% 
              filter(type=="transcript") %>%
              dplyr::select(transcript_id, transcript_name, transcript_type, gene_name, gene_type)
  ) %>% 
  relocate(c(transcript_id, transcript_name, transcript_type, gene_name, gene_type), .after = transcript_id) 

DTE_all_cores_table %>% 
  dplyr::count(transcript_set_SMG1i)

DTE_all_cores_table %>% 
  write_excel_csv(file.path("/home/volker/2023_SMG89_project", "Revision1_Table_S3_part2.csv"))
