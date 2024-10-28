#!/usr/bin/env Rscript

# Title: swish_DTE_gencode.v42.SIRVomeERCCome.R
# Objective: Main script for running Swish - DTE mode
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(RColorBrewer)
  library(tximeta)
  library(fishpond)
  library(SummarizedExperiment)
  library(ggrepel)
  library(ggridges)
  library(cowplot)
  library(limma)
  library(scales)
  library(extrafont)
})

# Load gencode annotation
gtf_gencode_df_short <- read_csv("/home/volker/reference/Gencode/gencode.v42.gtf_df_short.csv")

# # Initial setup - run only once
# indexDir <- file.path("/home/volker/reference/Transcriptome/gencode.v42.SIRVomeERCCome")
# fastaFTP <- c("/home/volker/reference/Gencode/v42/Raw_Files/gencode.v42.transcripts.fa.gz",
#               "home/volker/reference/spike_ins/SIRV_ERCC_2.fa.gz")
# gtfPath <- file.path("/home/volker/reference/gencode.v42.SIRVomeERCCome.annotation.gtf.gz")
# 
# makeLinkedTxome(indexDir=indexDir,
#                 source="GENCODE",
#                 organism="Homo sapiens",
#                 release="42",
#                 genome="GRCh38",
#                 fasta=fastaFTP,
#                 gtf=gtfPath,
#                 write=TRUE,
#                 jsonFile="/home/volker/reference/Transcriptome/gencode.v42.SIRVomeERCCome.json")


# Preparation steps -------------------------------------------------------

# Define and get arguments from bash input
arguments <- parse_args(OptionParser(), positional_arguments = 1)

# mydir is the first argument = working directory
mydir=arguments$args[1]

# Create folders if not exists
dir.create(file.path(mydir, "Swish", "DTE"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(mydir, "Plots", "Swish", "DTE"), showWarnings = FALSE, recursive = TRUE)

# Get samples and check if all files are present
samples <- read.table(file.path(mydir, "Samples.txt"), header = TRUE)

# Get unique conditions
condition <- samples %>% 
  distinct(condition) %>% 
  pull(condition)

# Generate data frame used for tximeta
coldata <- tibble(files = file.path(mydir, "Salmon", samples$sample, "quant.sf"))

# Supplement with sample IDs as "names"
coldata$names <- samples$sample

# Join with samples to get condition
coldata <- coldata %>% 
  left_join(samples,
            by = c("names" = "sample")) %>% 
  mutate(condition = fct_relevel(as_factor(condition),
                                 "control"))

# Check if all files are present - otherwise stop here
stopifnot("*** Not all salmon files are present! ***" = all(file.exists(coldata$files)))

##
# Individual analyses ----------------------------------------------------------
##

# Build for-loop to go through all conditions vs control
mycount=2
for (i in 2:length(condition)){
  
  cond <- paste(condition[mycount])
  
  coldata_2 <- coldata %>%
    filter(condition %in% c("control",cond)) %>% 
    mutate(condition = fct_drop(condition))
  
  # Create condition-specific subfolder
  dir.create(file.path(mydir,
                       "Swish",
                       "DTE",
                       paste(cond)),
             showWarnings = FALSE)
  
  # Use tximeta to import salmon data
  y <- tximeta(coldata_2,
               type = "salmon",
               txOut = TRUE,
               useHub = FALSE) # reads in counts and inf reps
  
  # Scale inferential replicate counts
  y <- scaleInfReps(y) # scales counts
  
  # Apply filtering step (default parameters for one condition)
  y <- labelKeep(y, x = "condition") # labels features to keep
  
  # Label spike-ins (ERCC and SIRV) as FALSE in keep column -> do not analyze them!
  mcols(y)$keep <- mcols(y)$keep & !str_detect(rownames(assays(y)[["counts"]]), 'ERCC') & !str_detect(rownames(assays(y)[["counts"]]), 'SIRV')
  
  set.seed(1)
  
  # Run the "real" swish analysis
  y <- swish(y, x="condition") # simplest Swish case
  
  # Obtain vector to be used as column names for count table
  rep_count_names <- coldata_2 %>% 
    group_by(condition) %>% 
    mutate(condition = case_when(condition == "control" ~ "control",
                                 TRUE ~ "condition2")) %>% 
    mutate(count_name = paste0("rep_",row_number(),"_",condition)) %>% 
    pull(count_name)
  
  rep_count_names <- c("transcript_id", rep_count_names)
  
  # From: https://support.bioconductor.org/p/p134531/
  # Obtain column median of Gibbs sample counts
  infReps <- assays(y)[ grep("infRep", assayNames(y)) ]
  infArray <- abind::abind( as.list(infReps), along=3 )
  
  infMed <- apply(infArray, 1:2, median)
  count_table <- as_tibble(infMed, rownames="transcript_id") 
  colnames(count_table) <- rep_count_names
  
  # Generate annotated final table
  y_df_gencode <- as_tibble(mcols(y)) %>% 
    dplyr::select(-gene_id) %>% 
    left_join(gtf_gencode_df_short %>% filter(type == "transcript") %>% dplyr::select(-type),
              by=c("tx_name" = "transcript_id")) %>% 
    dplyr::rename("transcript_id" = "tx_name") %>% 
    mutate(condition_2 = cond) %>% 
    left_join(count_table)
  
  # Generate output csv file
  message("Generate Swish DTE csv output file for ", paste(cond))
  write.csv(y_df_gencode, 
            file=(file.path(mydir,
                            "Swish",
                            "DTE",
                            cond, 
                            paste0(paste(cond),
                                   "_vs_control_Swish_results",
                                   ".csv"))))
  
  # Stats printing
  message("Swish DTE results significant qvalue < 0.0001 for  ", paste(cond))
  print(table(mcols(y)$qvalue < 0.0001))
  
  if (mycount == 2) {
    DTE_cat <- y_df_gencode
  } else {
    DTE_cat <- bind_rows(DTE_cat,
                         y_df_gencode)
  }
  
  # Counter increment
  mycount = mycount + 1
  
}

write.csv(DTE_cat, file=file.path(mydir, "Swish", "DTE", paste0("Swish_DTE_cat",
                                                             ".csv")))

## DTE NMD Ridgeplot -------------------------------------------------------

# Ridgeplot
DTE_ridgeplot <- DTE_cat %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other", 
                          transcript_type == "protein_coding" ~ "coding",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "lncRNA" ~ "lncRNA")) %>% 
  filter(type == "NMD") %>% 
  ggplot(aes(x=log2FC,
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

## DTE NMD Counts ----------------------------------------------------------

DTE_cat_counts_absolute <- DTE_cat %>% 
  mutate(condition_2 = fct_rev(fct_inorder(as_factor(condition_2)))) %>% 
  mutate(type = case_when(!transcript_type %in% c("protein_coding", "nonsense_mediated_decay", "lncRNA") ~ "other", 
                          transcript_type == "protein_coding" ~ "coding",
                          transcript_type == "nonsense_mediated_decay" ~ "NMD",
                          transcript_type == "lncRNA" ~ "lncRNA")) %>% 
  filter(type == "NMD") %>% 
  filter(keep == TRUE) %>% 
  mutate(UpDown = case_when(log2FC > 1 & qvalue < 0.0001 ~ "up",
                            log2FC < -1 & qvalue < 0.0001 ~ "down",
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
  labs(x = "Number of identified transcripts (n)",
       y = "",
       color="Significantly Up-/Down-regulated",
       caption=paste0("\nSwish-DTE - Filter: |log2FoldChange| > 1 & qvalue < 0.0001"))


DTE_ridgeplot_counts <- plot_grid(DTE_ridgeplot + 
                                    labs(title="Swish DTE - NMD-annotated-only! - gencode.v42.SIRVomeERCCome") +
                                    guides(fill = guide_legend(title.position = "top",
                                                               title.hjust = 0.5,
                                                               label.position = "bottom",
                                                               label.hjust = 0.5)),
                                  DTE_cat_counts_absolute + 
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

ggsave(file.path(mydir, "Plots", "Swish", "DTE", "Swish_DTE_summary.pdf"),
       width = 12,
       height=12,
       units = "cm",
       device=cairo_pdf,
       bg = "transparent")

writeLines(capture.output(sessionInfo()), paste0(mydir, "/Swish/DTE/Swish_DTE_session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
# Based on: http://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
