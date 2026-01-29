#!/usr/bin/env Rscript
## ============================================================================
## 05_Annotation_gtdb_suspects.R
## Purpose: Analyse and visualise metabolic annotations for GTDB-representative
##          methanotroph genomes, stratified by taxonomic groups and metabolic
##          gene presence. Produces publication-quality heatmaps showing methane
##          oxidation, carbon assimilation, and electron metabolism pathways.
## ============================================================================

library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(ggtree)
library(ape)
library(ggh4x)
library(viridis)

setwd("path/to/your/repo/MFD_methanotrophs_DK/")

## ---- Load KEGG annotation table and create tree labels ----
# KEGG_25_12_02.rds: merged functional annotations with GTDB taxonomy (from 04_Annotation.R)
# tree_label: formatted display label combining Family/Genus/Species with genome ID for heatmap y-axis
KEGG<-readRDS("./output/KEGG_25_12_02.rds")%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))

## Load color palette and gene metadata
# methane_colors: predefined RGB color map for enzyme types (pMMO, sMMO, H4MPT, etc.)
methane_colors<-readRDS("palette_methane_colors.rds")%>%
  unlist()
names(methane_colors)[names(methane_colors) == "-"] <- "Formate oxidation"

## ---- Load and standardize KO metadata (gene labels, pathways) ----
# KO_KSK_25_08_12.xlsx: KEGG Orthology metadata with pathway mappings and gene names
# Create abbreviated pathway names and gene labels for heatmap visualization; split Metabolic_step for multi-line facet labels
KOs<-read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

## ---- Prepare KEGG data for visualization ----
# Sort by genus for consistency; create combined gene+step label for heatmap x-axis annotation
KEGG_sort<-KEGG%>%
  arrange(Genus)

KEGG<-KEGG%>%
  mutate (gene_label = paste0(Gene_collapsed , ' ',  Pathway_step))

## ---- Create taxonomic labels for heatmap visualization ----
# class_label: combined Class;Family;Species|genome for detailed sample identification in heatmaps
KEGG<-KEGG%>%
  mutate(class_label = paste0(Class, ";", Family,";", tree_label))

## ---- Define exclusion lists for contaminants and non-methanotrophs ----
# genes_remove_string: gene families from non-methane-oxidizers (nitrifiers, amoA, non-specific MOX)
# genus_remove_string: non-methanotroph genera to exclude (nitrifiers, nitrospirae, etc.)
genes_remove_string<-c("Nitrosomonadaceae (1/3)", "Homologous_MO (1/3)", "Nitrosomonas (1/3)", "Propane_MO_Actino_cluster (1/3)",
                       "Homologous_pmoA (1/3)", "Cycloclasticus (1/3)", "Betaproteobacteria_amoA (1/3)",
                       "Nitrosococcus (1/3)", "Nitrospira_clade_B (1/3)", "Actinobacteria (1/3)", "Homologous_Rhodopila (1/3)")

genus_remove_string<-c("g__Methylocella", "g__Methylomonas", "g__Methylocystis", "g__Methylomirabilis", "g__Methylosinus", "g__Methylobacter_C",
                       "g__Methylovulum", "g__Methylobacter_A", "g__Crenothrix", "g__Methylosarcina", "g__Methylobacter", "g__Methyloferula", "g__Methylocaldum", "g__Nitrospira_D",
                       "g__Bradyrhizobium", "g__Nitrosospira", "g__Methylobacter_B", "g__Nitrosomonas", "g__Methyloglobulus", "g__Methylomicrobium", "g__Methylotuvimicrobium", 
                       "g__Methylicorpusculum", "g__Methylomarinum", "g__Methylococcus", "g__Methylohalobius", "g__Methylacidimicrobium", "g__Methylacidiphilum", "g__Methylacidithermus",
                       "g__Methyloprofundus")

## ---- Filter genomes  ----
# genomes_keep: GTDB representatives with Module C1 (methane oxidation) + mmoX/pmoA detection
# Exclude: MFD MAGs, non-methane oxidizers (nitrifiers), low-quality gene assignments
genomes_keep<-KEGG %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(!is.na(Genus))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(!KO %in% c("Root; Homologous_pmoA; Mycobacterium"))%>%
  filter(!gene_label %in% genes_remove_string)%>%
  filter(presence==1)%>%
  filter(!Genus %in% genus_remove_string)%>%
  distinct()

## ---- Extract unique methane-oxidation genes from GTDB methanotrophs ----
# genes_keep: unique mmoX/pmoA gene types present in Module C1 methanotrophs
genes_keep<-genomes_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()


unique(CM$data$genome) 

## ---- Filter and visualize Binatia/Gemmatimonadetes (candidate lineages) ----

genomes_Bin_keep<-KEGG %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(!is.na(Genus))%>%
  filter(Class%in%c("c__Binatia", "c__Gemmatimonadetes"))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(!KO %in% c("Root; Homologous_pmoA; Mycobacterium"))%>%
  filter(!gene_label %in% genes_remove_string)%>%
  filter(presence==1)%>%
  filter(!Genus %in% genus_remove_string)%>%
  distinct()

## ---- Extract genes unique to Binatia/Gemmatimonadetes for visualization ----
genes_Bin_keep<-genomes_Bin_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()


genes_remove<-KEGG%>%  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(!gene_label %in% genes_Bin_keep$gene_label)%>%
  select(gene_label)%>%distinct()

## ---- Create Binatia/Gemmatimonadetes heatmap showing methane oxidation genes ----
pl_bin_tusc <- KEGG %>% 
  filter(genome %in% genomes_Bin_keep$genome)%>%
  mutate(fam_label=paste0(Family, " ", Species))%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = fam_label))+
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "pMMO"="#b3943c",
    "sMMO"="darkgreen",
    "calcium (mxa)"="#5e4fa2",
    "lanthanide (xoxF)"="#c46ca1",
    "putative"="grey75",
    "H4MPT"="#d97512",
    "GSH"="turquoise4",
    "H4F"="darkred",
    "Serine cycle"="darkred",
    "RuMP cycle"="darkblue",
    "CBB cycle"="#679b60",
    "Formate oxidation"="#679b60"
  )) +
  facet_nested(Order~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_text(size=5),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size=5, face="bold", angle=0, hjust=0),
        #  strip.text.y = element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0.06, "cm", data = NULL),
        panel.spacing.y = unit(0.00, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )

pl_bin_tusc

## ---- Export Binatia/Gemmatimonadetes heatmap ----
# Save high-resolution PNG and SVG for publication
ggsave("./output/GTDB_Bin_TUSC_25_12_10.png",pl_bin_tusc,
       units = c("mm"),
       height = 90,
       width = 190,
       dpi=300)



ggsave("./output/GTDB_Bin_TUSC_25_12_10.svg",pl_bin_tusc,
       units = c("mm"),
       height = 90,
       width = 190,
       dpi=300)






## ---- Visualize methanotrophs (excluding Binatia and Methylomonadaceae) ----
pl <- KEGG %>% 
  filter(genome %in% genomes_keep$genome)%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  filter(!Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae"))%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = Species))+
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "pMMO"="#b3943c",
    "sMMO"="darkgreen",
    "calcium (mxa)"="#5e4fa2",
    "lanthanide (xoxF)"="#c46ca1",
    "putative"="grey75",
    "H4MPT"="#d97512",
    "GSH"="turquoise4",
    "H4F"="darkred",
    "Serine cycle"="darkred",
    "RuMP cycle"="darkblue",
    "CBB cycle"="#679b60",
    "Formate oxidation"="#679b60"
  )) +
  facet_nested(.~Metabolic_step, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size=8, face="bold", angle=0),
        #  strip.text.y = element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )

pl

## ---- Filter genes for broader methanotroph visualization (excluding Methylococcales) ----
# genes_keep: genes present in all GTDB methanotrophs except Binatia and Methylococcales families
genes_keep<-genomes_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  filter(!Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae",  "f__JACCXJ01", "f__Desulfobaccaceae", "f__JAFGDC01"))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  filter(!gene_label=="Homologous_Binatales (1/3)")%>%
  # mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  distinct()

genes_remove<-KEGG%>%  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(!gene_label %in% genes_keep$gene_label)%>%
  select(gene_label)%>%distinct()

pl_2 <- KEGG %>% 
  filter(genome %in% genomes_keep$genome)%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  mutate(fam_label=paste0(gsub("f__", "", Family), " ", Species))%>%
  filter(!Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae", "f__JACCXJ01", "f__Desulfobaccaceae", "f__JAFGDC01"))%>%
  mutate(Class=gsub("c__", "", Class))%>%  mutate(Class=gsub("proteobacteria", "-\nproteobacteria", Class))%>%
  mutate(Class=gsub("cocco", "-\ncocco", Class))%>%
  mutate(Class=gsub("mycetes", "-\nmycetes", Class))%>%
  mutate(gene_label=if_else(gene_label=="Nevskiales (1/3)", paste0("Nevskiales_Macondimonas_put_pmoA (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Burkholderiaceae (1/3)", paste0("Nevskiales_Macondimonas_put_pmoA (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Methylocystis_pmoA1_3 (1/3)", paste0("Methylosinus_Methylocystis_pmoA1 (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Methylobacter_clade_2 (1/3)", paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(grepl("Methylotenera|JABFRO01", gene_label), paste0("Methylococcales_like_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Methyloferula_Methylovirgula (1/3)", paste0("Rhizobiales_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Rhodopila_Put_mmoX (1/3)", paste0("Rhizobiales_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Mycobacterium (1/3)", paste0("Mycobacterium_Methanotrophicium_cluster (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label=="Methylococcales_unknown (1/3)", paste0("Methylococcaceae_pmoA (1/3)"), gene_label))%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = fam_label))+
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "pMMO"="#b3943c",
    "sMMO"="darkgreen",
    "calcium (mxa)"="#5e4fa2",
    "lanthanide (xoxF)"="#c46ca1",
    "putative"="grey75",
    "H4MPT"="#d97512",
    "GSH"="turquoise4",
    "H4F"="darkred",
    "Serine cycle"="darkred",
    "RuMP cycle"="darkblue",
    "CBB cycle"="#679b60",
    "Formate oxidation"="#679b60"
  )) +
  facet_nested(Class~Metabolic_step, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5),
        axis.text.y = element_text(size=5, margin = margin(l = -1)),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin = margin(b = -0.5)),  # Size for x-axis facet labels
        strip.text.y = element_text(size=5, angle=0, hjust=0, margin = margin(r = -2)),
        #  strip.text.y = element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0.03, "cm", data = NULL),
        panel.spacing.y = unit(0.08, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white"), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )

## ---- Export detailed multi-family heatmap ----
ggsave("./output/metabolism_Extended_figs/GTDB_gtdb_sub_met_25_12_10.png",pl_2,
       units = c("mm"),
       height = 210,
       width = 185,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/GTDB_gtdb_sub_met_25_12_10.svg",pl_2,
       units = c("mm"),
       height = 210,
       width = 185,
       dpi=300)

## ---- Filter for Methylococcales methanotrophs and create focused heatmap ----
# Subset genomes to Methylococcales order + JACCXJ01 family 

genes_keep<-genomes_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  filter(Order %in% c("o__Methylococcales") | Family %in% c("f__JACCXJ01"))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()

## ---- Identify Methylococcales-specific genes  ----
print(n=26, KEGG%>%filter(gene_label %in% genes_keep$gene_label)%>%
        select(KO, gene_label)%>%distinct())

## ---- Methylococcales-focused heatmap  ----
pl_3 <- KEGG %>% 
  filter(genome %in% genomes_keep$genome)%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  mutate(fam_label=paste0(gsub("f__", "", Family), " ", Species))%>%
  #filter(Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae"))%>%
  filter(Order %in% c("o__Methylococcales") | Family %in% c("f__JACCXJ01"))%>%
  mutate(Order=gsub("o__", "", Order))%>% 
  mutate(gene_label=if_else(grepl("Methylobacter_clade_1|CAIQWF01|JANEMG01|SXIZ01", gene_label), paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label %in% c("Methylogaea (1/3)", "Methylococcales_unknown (1/3)"), paste0("Methylococcaceae_pmoA (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label %in% c("Methylobacter_C_mmoX (1/3)", "Methyloprofundus_WTBX01_clade (1/3)", "UBA10906_mmoX (1/3)", "JABFRC01 (1/3)", "UBA10906_mmoX (1/3)", "UBA4132 (1/3)"), paste0("Methylomonadaceae_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(gene_label %in% c("Methylococcales_mmoX (1/3)"), paste0("Methyloccocaceae_mmoX (1/3)"), gene_label))%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = fam_label))+
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "pMMO"="#b3943c",
    "sMMO"="darkgreen",
    "calcium (mxa)"="#5e4fa2",
    "lanthanide (xoxF)"="#c46ca1",
    "putative"="grey75",
    "H4MPT"="#d97512",
    "GSH"="turquoise4",
    "H4F"="darkred",
    "Serine cycle"="darkred",
    "RuMP cycle"="darkblue",
    "CBB cycle"="#679b60",
    "Formate oxidation"="#679b60"
  )) +
  facet_nested(Order~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_text(size=5, margin = margin(l = -1)),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin = margin(b = -0.5)),  # Size for x-axis facet labels
        strip.text.y = element_text(size=5, angle=0, hjust=0, margin = margin(r = -2)),
        #  strip.text.y = element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0.03, "cm", data = NULL),
        panel.spacing.y = unit(0.08, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white"), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )

pl_3


ggsave("output/GTDB_Methylococcales_gtdb_sub_met_25_12_10.png",pl_3,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)



ggsave("output/GTDB_Methylococcales_gtdb_sub_met_25_12_10.svg",pl_3,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)



## ---- MFD Binatia/Gemmatimonadetes heatmap (MFD-TUSC/Bin candidates) ----
# Subset MFD MAGs to Binatia and Gemmatimonadetes classes for visualization
genomes_Bin_keep<-KEGG %>% 
  filter(grepl("MFD|LIB", genome))%>%
  filter(!is.na(Genus))%>%
  filter(Class%in%c("c__Binatia", "c__Gemmatimonadetes"))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(!KO %in% c("Root; Homologous_pmoA; Mycobacterium"))%>%
  filter(!gene_label %in% genes_remove_string)%>%
  filter(presence==1)%>%
  filter(!Genus %in% genus_remove_string)%>%
  distinct()


genes_Bin_keep<-genomes_Bin_keep %>% 
  filter(grepl("MFD|LIB", genome))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()

## Identify genes for exclusion in other visualizations
genes_remove<-KEGG%>%  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(!gene_label %in% genes_Bin_keep$gene_label)%>%
  select(gene_label)%>%distinct()

## ---- MFD Binatia/Gemmatimonadetes metabolic heatmap ----
pl_bin_tusc <- KEGG %>% 
  mutate(Order=gsub("o__", "", Order))%>% 
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = label))+
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "pMMO"="#b3943c",
    "sMMO"="darkgreen",
    "calcium (mxa)"="#5e4fa2",
    "lanthanide (xoxF)"="#c46ca1",
    "putative"="grey75",
    "H4MPT"="#d97512",
    "GSH"="turquoise4",
    "H4F"="darkred",
    "Serine cycle"="darkred",
    "RuMP cycle"="darkblue",
    "CBB cycle"="#679b60",
    "Formate oxidation"="#679b60"
  )) +
  facet_nested(Order~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_text(size=5, margin = margin(l = -1)),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin = margin(b = -0.5)),  # Size for x-axis facet labels
        strip.text.y = element_text(size=5, angle=0, hjust=0, margin = margin(r = -2)),
        #  strip.text.y = element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0.03, "cm", data = NULL),
        panel.spacing.y = unit(0.08, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white"), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )

pl_bin_tusc

## ---- Export MFD Binatia/Gemmatimonadetes heatmap ----

ggsave("output/MFD_Bin_TUSC_sub_met_25_12_10.png",pl_bin_tusc,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

ggsave("output/MFD_Bin_TUSC_sub_met_25_12_10.svg",pl_bin_tusc,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

## ---- Electron donor/acceptor metabolism heatmap (MFD Binatia/Gemmatimonadetes) ----
# Focus on nitrogen, sulfur, energy (H2/CO oxidation) pathways; exclude methane oxidation
NS <- KEGG %>% 
  filter(genome %in% genomes_Bin_keep$genome)%>%
  mutate(Order=gsub("o__", "", Order))%>% 
  filter(!grepl("Dsr|dsr|apr", Type))%>%
  filter(!grepl("hox|hup|anfG|SOX|sox|nxrA\\/|nxrB\\/", Gene_collapsed))%>%
  filter(Module %in% c("TCA cycle", "PHB", "Nitrogen", "Sulphur", "PPP", "CO and H2", "Glycolysis")) %>% 
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Gene_collapsed = fct_relevel(Gene_collapsed, KOs$Gene_collapsed)) %>%
  mutate(Type = factor(Type, levels = unique(KOs$Type))) %>%
  ggplot(aes(x = Gene_collapsed, y = label))+ 
  geom_tile(aes(fill = if_else(is.na(presence), "white", Type)), color="grey90", linewidth=0.1 )+
  #  geom_tile(aes(fill = if_else(is.na(presence),  "white", Pathway)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "Carbon monoxide oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "Assimilatory sulphate reduction"="#0DB7A8",
    "sat"="#6C961E",
    "aprAB" = "#18B862",
    "DsrABL" = "#77AA79",
    "DsrC"= "#40B645",
    "dsrMKJOP"= "lightblue4",
    "dsrEFH"="#5B815B",
    "SOX"="#92C76C",
    "TCA cycle"="#C7BB6D",
    "PHB"="#BC8615",
    "acetate uptake"="#C7966D",
    "EMC"="#BD5531",
    "PPP"="#AD8ACC",
    "Glycolysis"="#1237B7"))+
  facet_nested(Order ~factor(Module, levels=c("CO and H2", "Nitrogen", "Sulphur", "TCA cycle", "PHB", "EMC / acetate","PPP", "Glycolysis")), scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_text(size=5, margin = margin(l = -1)),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin = margin(b = -0.5)),  # Size for x-axis facet labels
        strip.text.y = element_text(size=5, angle=0, hjust=0, margin = margin(r = -2)),
        strip.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0.03, "cm", data = NULL),
        panel.spacing.y = unit(0.08, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white"), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
  )


NS

ggsave("output/MFD_Bin_TUSC_Nitro_sulphur_25_12_10.png",NS,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

ggsave("output/MFD_Bin_TUSC_Nitro_sulphur_25_12_10.svg",NS,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)


