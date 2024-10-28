#!/usr/bin/env Rscript

# Title: SMG89_Proteomics_Analysis
# Objective: Proteomics analysis for "SMG1:SMG8:SMG9-complex integrity maintains robustness of nonsense-mediated mRNA decay"
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

##
# Load libraries ----------------------------------------------------------
##

library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggh4x)
library(readxl)
library(janitor)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)

##
# Setup -------------------------------------------------------------------
##

setwd("W:/Projects/02_SMG8-9/GitHub/2024_SMG189/data")

# Read in data and clean-up
data <- read_excel("SMG89_proteomics.xlsx", 
                   sheet = "matrix17 - PGs",
                   skip = 1) %>% 
  clean_names() %>% 
  dplyr::select(-c("log_students_t_test_p_value_8ko_01_control_2",
                   "students_t_test_q_value_8ko_01_control_2",
                   "students_t_test_difference_8ko_01_control_2",
                   "students_t_test_test_statistic_8ko_01_control_2"
                   ))

# List columns
colnames(data)

# Define genes of interest
special_NMD <- c("SMG1", "SMG5", "SMG6", "SMG7", "SMG8", "SMG9", "UPF1", "UPF2", "UPF3A", "UPF3B")
special_EJC <- c("EIF4A3", "RBM8A", "MAGOH", "MAGOHB", "RNPS1", "ACIN1", "PNN", "SAP18", "CASC3")
special_other <- c("STAU1", "STAU2")
labels <- c(special_NMD, special_EJC, special_other)   

# Prepare data for plotting
SMG89_degron_proteomics_all <- data %>% 
  dplyr::select(majority_protein_i_ds,
                gene_names,
                starts_with("students_t_test_difference")) %>% 
  pivot_longer(cols = c(starts_with("students_t_test_difference")),
               names_to = c(".value", "comparison"),
               names_sep = "difference_") %>% 
  filter(comparison != "8ko_01_control_2") %>% 
  rename("log2FC" = "students_t_test_") %>% 
  left_join(data %>% 
              dplyr::select(majority_protein_i_ds,
                            gene_names,
                            starts_with("students_t_test_q_value")) %>% 
              pivot_longer(cols = c(starts_with("students_t_test_q_value")),
                           names_to = c(".value", "comparison"),
                           names_sep = "q_value_") %>% 
              filter(comparison != "8ko_01_control_2") %>% 
              rename("padj" = "students_t_test_")) %>% 
  mutate(log2FC = case_when(comparison == "control_wt_01" ~ -log2FC,
                            comparison == "wt_8ko" ~ -log2FC,
                            comparison == "wt_9ko" ~ -log2FC,
                            comparison == "8ko_8ko_01" ~ -log2FC,
                            comparison == "9ko_9ko_01" ~ -log2FC,
                            TRUE ~ log2FC)) %>% 
  mutate(comparison = case_when(comparison == "control_wt_01" ~ "wt_01_control",
                                comparison == "wt_8ko" ~ "8ko_wt",
                                comparison == "wt_9ko" ~ "9ko_wt",
                                comparison == "8ko_8ko_01" ~ "8ko_01_8ko",
                                comparison == "9ko_9ko_01" ~ "9ko_01_9ko",
                                TRUE ~ comparison)) %>% 
  mutate(comparison = case_when(comparison == "wt_control" ~ "WT 0µM",
                                comparison == "wt_01_control" ~ "WT 0.1µM",
                                comparison == "wt_1_control" ~ "WT 1µM",
                                comparison == "8ko_control" ~ "SMG8 KO 0µM",
                                comparison == "8ko_01_control" ~ "SMG8 KO 0.1µM",
                                comparison == "9ko_control" ~ "SMG9 KO 0µM",
                                comparison == "9ko_01_control" ~ "SMG9 KO 0.1µM",
                                comparison == "wt_01_wt" ~ "WT 0.1µM/0µM",
                                comparison == "wt_1_wt" ~ "WT 1µM/0µM",
                                comparison == "8ko_wt" ~ "SMG8 KO/WT 0µM",
                                comparison == "9ko_wt" ~ "SMG9 KO/WT 0µM",
                                comparison == "8ko_01_wt_01" ~ "SMG8 KO/WT 0.1µM",
                                comparison == "9ko_01_wt_01" ~ "SMG9 KO/WT 0.1µM",
                                comparison == "wt_1_wt_01" ~ "WT 1µM/0.1µM",
                                comparison == "8ko_01_8ko" ~ "SMG8 KO 0.1/0µM",
                                comparison == "9ko_01_9ko" ~ "SMG9 KO 0.1/0µM",
                                ))
##
# Lollipop plots ----------------------------------------------------------------
##

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(gene_names %in% special_NMD) %>% 
  dplyr::filter(comparison %in% c("WT 1µM/0µM",
                                  "SMG8 KO/WT 0µM",
                                  "SMG9 KO/WT 0µM")) %>% 
  mutate(gene_names = fct_rev(fct_relevel(gene_names,
                                   "UPF1", "UPF2", "UPF3B", "SMG1", "SMG5", "SMG6", "SMG7", "SMG8", "SMG9" 
  ))) %>% 
  mutate(comparison = (fct_relevel(comparison,
                                          "WT 1µM/0µM",
                                          "SMG8 KO/WT 0µM",
                                          "SMG9 KO/WT 0µM"
  ))) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(x=log2FC,
             y=gene_names)) +
  theme(strip.background = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=6, color="black"),
        axis.title=element_text(size=6),
        strip.text = element_text(size=6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line.y = element_blank(), 
        axis.line.x = element_line(colour = "black", linewidth = 0.1), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  geom_vline(xintercept = 0, linewidth=0.2) +
  geom_linerange(aes(xmin = 0, xmax = log2FC),
                 alpha=0.5,
                 color="darkgray") +
  geom_point(aes(fill=log2FC,
                 shape=significant,
                 stroke=stroke,
                 size=-log10(padj))) +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
  xlim(c(-4,4)) +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       limits = c(-4,4),
                       midpoint = 0,
                       na.value = "grey80") +
  facet_wrap(~comparison) +
  labs(y="") +
  force_panelsizes(cols = unit(25, "mm"),
                   rows = unit(25, "mm"))

ggsave(filename = "Revision1_Lollipop_Tannenbaums_2ndEdition.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

##
# Heat-Dot-Maps -----------------------------------------------------------
##

NMD_names <- c("UPF1",
               "UPF2",
               "UPF3A",
               "UPF3B",
               "SMG1",
               "SMG5",
               "SMG6",
               "SMG7",
               "SMG8",
               "SMG9")

EJC_names <- c("EIF4A3",
               "MAGOH;MAGOHB",
               "MAGOH",
               "MAGOHB",
               "RBM8A",
               "CASC3")

STAU_names <- c("STAU1",
                "STAU2")

Other_names <- c("DHX34",
                 "MOV10",
                 "G3BP1",
                 "G3BP2",
                 "GSPT1",
                 "GSPT2",
                 "RUVBL1",
                 "RUVBL2",
                 "PABPC1",
                 "PABPC4",
                 "YTHDF1",
                 "ETF1",
                 "NBAS",
                 "ATM",
                 "SLBP",
                 "PNRC2",
                 "ATR",
                 "ZC3H12A",
                 "TRIM71",
                 "GW182",
                 "AGO2",
                 "DCP1A",
                 "DCP1B",
                 "XRN1"
)

Names_combined <- c(NMD_names,
                    EJC_names,
                    STAU_names)

Names_combined_other <- c("UPF1",
                          Other_names)


## NMD_EJC_SMD -------------------------------------------------------------

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(comparison %in% c("WT 0µM",
                                  "WT 0.1µM",
                                  "WT 1µM",
                                  "SMG8 KO 0µM",
                                  "SMG8 KO 0.1µM",
                                  "SMG9 KO 0µM",
                                  "SMG9 KO 0.1µM")) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "WT 0µM",
                                          "WT 0.1µM",
                                          "WT 1µM",
                                          "SMG8 KO 0µM",
                                          "SMG8 KO 0.1µM",
                                          "SMG9 KO 0µM",
                                          "SMG9 KO 0.1µM"))) %>% 
  filter(gene_names %in% Names_combined) %>% 
  mutate(gene_names = (fct_relevel(gene_names,
                              "UPF1",
                              "UPF2",
                              "UPF3B",
                              "SMG1",
                              "SMG5",
                              "SMG6",
                              "SMG7",
                              "SMG8",
                              "SMG9",
                              "EIF4A3",
                              "MAGOH",
                              "MAGOHB",
                              "RBM8A",
                              "CASC3",
                              "STAU1",
                              "STAU2"))) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(y=comparison,
             x=gene_names,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    shape=significant,
    stroke=stroke,
    size=-log10(padj)),
    color="black") +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
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
        legend.position = "right",
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
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(15*3+0.6, "mm"),
                   rows = unit(7*3+0.6, "mm"))

ggsave(filename = "Revision1_HeatDotMap_NMD_control.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")


## Other interactors -------------------------------------------------------

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(comparison %in% c("WT 0µM",
                                  "WT 0.1µM",
                                  "WT 1µM",
                                  "SMG8 KO 0µM",
                                  "SMG8 KO 0.1µM",
                                  "SMG9 KO 0µM",
                                  "SMG9 KO 0.1µM")) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "WT 0µM",
                                          "WT 0.1µM",
                                          "WT 1µM",
                                          "SMG8 KO 0µM",
                                          "SMG8 KO 0.1µM",
                                          "SMG9 KO 0µM",
                                          "SMG9 KO 0.1µM"))) %>% 
  filter(gene_names %in% Names_combined_other) %>% 
  mutate(gene_names = (fct_relevel(gene_names,
                              "UPF1",
                              "GSPT1",
                              "GSPT2",
                              "ETF1",
                              "PABPC1",
                              "PABPC4",
                              "DHX34",
                              "NBAS",
                              "RUVBL1",
                              "RUVBL2",
                              "MOV10",
                              "G3BP1",
                              "G3BP2",
                              "DCP1A",
                              "DCP1B",
                              "XRN1",
                              "YTHDF1",
                              "ATM",
                              "SLBP",
                              "PNRC2",
                              "ATR",
                              "ZC3H12A",
                              "TRIM71",
                              "GW182",
                              "AGO2"))) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(y=comparison,
             x=gene_names,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    shape=significant,
    stroke=stroke,
    size=-log10(padj)),
    color="black") +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
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
        legend.position = "right",
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
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(13*3+0.6, "mm"),
                   rows = unit(7*3+0.6, "mm"))

ggsave(filename = "Revision1_HeatDotMap_other_control.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

# Heat-Dot-Maps all -----------------------------------------------------------

## NMD_EJC_SMD -------------------------------------------------------------

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(comparison %in% c("WT 0.1µM/0µM",
                                  "WT 1µM/0µM",
                                  "SMG8 KO/WT 0µM",
                                  "SMG9 KO/WT 0µM",
                                  "SMG8 KO/WT 0.1µM",
                                  "SMG9 KO/WT 0.1µM",
                                  "WT 1µM/0.1µM",
                                  "SMG8 KO 0.1/0µM",
                                  "SMG9 KO 0.1/0µM")) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "WT 0.1µM/0µM",
                                          "WT 1µM/0µM",
                                          "SMG8 KO/WT 0µM",
                                          "SMG9 KO/WT 0µM",
                                          "SMG8 KO/WT 0.1µM",
                                          "SMG9 KO/WT 0.1µM",
                                          "WT 1µM/0.1µM",
                                          "SMG8 KO 0.1/0µM",
                                          "SMG9 KO 0.1/0µM"))) %>% 
  filter(gene_names %in% Names_combined) %>% 
  mutate(gene_names = (fct_relevel(gene_names,
                                   "UPF1",
                                   "UPF2",
                                   "UPF3B",
                                   "SMG1",
                                   "SMG5",
                                   "SMG6",
                                   "SMG7",
                                   "SMG8",
                                   "SMG9",
                                   "EIF4A3",
                                   "MAGOH",
                                   "MAGOHB",
                                   "RBM8A",
                                   "CASC3",
                                   "STAU1",
                                   "STAU2"))) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(y=comparison,
             x=gene_names,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    shape=significant,
    stroke=stroke,
    size=-log10(padj)),
    color="black") +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
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
        legend.position = "right",
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
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(15*3+0.6, "mm"),
                   rows = unit(9*3+0.6, "mm"))

ggsave(filename = "Revision1_HeatDotMap_NMD_all_2nd.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")

## Other interactors -------------------------------------------------------

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(comparison %in% c("WT 0.1µM/0µM",
                                  "WT 1µM/0µM",
                                  "SMG8 KO/WT 0µM",
                                  "SMG9 KO/WT 0µM",
                                  "SMG8 KO/WT 0.1µM",
                                  "SMG9 KO/WT 0.1µM",
                                  "WT 1µM/0.1µM")) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "WT 0.1µM/0µM",
                                          "WT 1µM/0µM",
                                          "SMG8 KO/WT 0µM",
                                          "SMG9 KO/WT 0µM",
                                          "SMG8 KO/WT 0.1µM",
                                          "SMG9 KO/WT 0.1µM",
                                          "WT 1µM/0.1µM"))) %>% 
  filter(gene_names %in% Names_combined_other) %>% 
  mutate(gene_names = (fct_relevel(gene_names,
                                   "UPF1",
                                   "GSPT1",
                                   "GSPT2",
                                   "ETF1",
                                   "PABPC1",
                                   "PABPC4",
                                   "DHX34",
                                   "NBAS",
                                   "RUVBL1",
                                   "RUVBL2",
                                   "MOV10",
                                   "G3BP1",
                                   "G3BP2",
                                   "DCP1A",
                                   "DCP1B",
                                   "XRN1",
                                   "YTHDF1",
                                   "ATM",
                                   "SLBP",
                                   "PNRC2",
                                   "ATR",
                                   "ZC3H12A",
                                   "TRIM71",
                                   "GW182",
                                   "AGO2"))) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(y=comparison,
             x=gene_names,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    shape=significant,
    stroke=stroke,
    size=-log10(padj)),
    color="black") +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
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
        legend.position = "right",
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
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(13*3+0.6, "mm"),
                   rows = unit(7*3+0.6, "mm"))

ggsave(filename = "Revision1_HeatDotMap_other_all.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")


# For guide scaling -------------------------------------------------------

SMG89_degron_proteomics_all %>% 
  mutate(padj = replace(padj, padj == 0, 0.0001)) %>% 
  dplyr::filter(comparison %in% c("WT 1µM/0µM")) %>% 
  mutate(comparison = fct_rev(fct_relevel(comparison,
                                          "WT 1µM/0µM"))) %>% 
  filter(gene_names %in% c("UPF1",
                           "UPF2",
                           "UPF3B",
                           "SMG1",
                           "SMG5")) %>% 
  mutate(gene_names = (fct_relevel(gene_names,
                                   "UPF1",
                                   "UPF2",
                                   "UPF3B",
                                   "SMG1",
                                   "SMG5"))) %>% 
  mutate(padj = case_when(gene_names == "UPF1" ~ 1,
                          gene_names == "UPF2" ~ 0.1,
                          gene_names == "UPF3B" ~ 0.01,
                          gene_names == "SMG1" ~ 0.001,
                          gene_names == "SMG5" ~ 0.0001)) %>% 
  mutate(significant = case_when(padj < 0.01 ~ 21,
                                 TRUE ~ 22)) %>% 
  mutate(stroke = case_when(significant == 21 ~ 0.5,
                            significant == 22 ~ 0.2)) %>% 
  ggplot(aes(y=comparison,
             x=gene_names,
             fill=log2FC)) +
  geom_tile(color = "black",
            fill="white",
            lwd = 0.1,
            linetype = 1) +
  geom_point(aes(
    fill=log2FC,
    shape=significant,
    stroke=stroke,
    size=-log10(padj)),
    color="black") +
  scale_fill_gradient2(low = ("#2166AC"),
                       mid = "white",
                       high = ("#B2182B"),
                       midpoint = 0,
                       na.value = "grey80") +
  scale_size(range = c(1, 3)) +
  scale_shape_identity() +
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
        legend.position = "right",
        legend.margin=margin(c(0.25,0.25,0.25,0.25)),
        plot.title = element_text(size = 6),
        plot.subtitle = element_text(size = 6),
        plot.caption = element_text(size = 5),
        text=element_text(family="Arial"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(colour = "black", linewidth = 0.1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="TEST!",
       y="TEST!",
       fill="log2FC") +
  coord_fixed(ratio=1) + 
  theme(legend.direction = "horizontal", legend.box = "vertical") +
  guides(fill = guide_colourbar(order = 1,
                                title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                label.vjust = 0.5,
                                label.theme = element_text(angle = 90, size = 6)),
         size = guide_legend(order = 2,
                             title.position = "top",
                             label.position = "bottom")) +
  theme(axis.ticks.y = element_blank()) +
  theme(aspect.ratio = 1) +
  force_panelsizes(cols = unit(5*3+0.6, "mm"),
                   rows = unit(1*3+0.6, "mm"))

ggsave(filename = "Revision1_HeatDotMap_Scale.pdf", 
       width = 20,
       height = 20,
       units = "cm", 
       device= cairo_pdf, 
       bg = "transparent")