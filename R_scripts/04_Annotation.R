#!/usr/bin/env Rscript

## 04_Annotation.R
## Purpose: parse and merge KEGG annotations, custom HMMs (pmoA/mmoX), hydrogenase
##          and DsrMKJOP hits; create presence/absence tables and heatmap plots.
## Usage: load KEGG and custom HMM tables, merge with taxonomy, and export
##        cleaned tables for genome-level metabolic analyses.

library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)

setwd("path/to/your/repo/MFD_methanotrophs_DK/")

## ---- Load HMM-derived KOs (pmoA/mmoX) and merge with KEGG annotations ----
# HMM_ko contains custom HMM hits for pmoA and mmoX genes from various datasets
HMM_ko<-readRDS("data/HMM_KOs_25_06_06.rds")

## Clean whitespace and merge with KEGG annotation table
HMM_ko$ko_id <- gsub("\u00A0", " ", HMM_ko$ko_id)


KEGG<-vroom("data/DRAM_output.tsv",delim="\t", col_select=c("fasta", "ko_id")) %>%
  mutate(fasta=gsub("_genomic","",fasta))%>%
  filter(!is.na(ko_id))%>%
  rbind(HMM_ko)



genomes_KEGG <- unique(KEGG$fasta)

## ---- Load KO-to-pathway mapping table ----
# KO_KSK file contains metabolic pathway metadata for each KEGG KO
KOs<-read_excel("data/KO_KSK_methanotroph_paper.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed, ' - ',  Pathway_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)

## Clean up pathway names for compact plotting
KOs<-KOs%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))


## ---- Create presence/absence table for all KOs across all genomes ----
# Cartesian join of KO metadata with genomes, then left join with actual hits
df <- data.frame()
for (i in 1:length(genomes_KEGG)) {
  d <- KOs %>% mutate(genome = genomes_KEGG[i])
  df = rbind(df, d)
}
# Mark presence as 1 where hit found, then join and fill NA with 0 (absence)
KEGG <- KEGG %>% mutate(presence = 1)
KEGG <- left_join(df, KEGG, by = c('KO'='ko_id','genome' ='fasta'), relationship = "many-to-many")
KEGG$presence[is.na(KEGG$presence)] <- 0


## ---- Resolve conflicts when duplicate KO-gene pairs have mixed presence (0,1) ----
# Keep rows with presence==1 (they represent confirmed hits) when duplicates exist
cleaned_KEGG <- KEGG %>%
  filter(presence %in% c(0, 1)) %>%
  group_by(genome, gene_label) %>%
  filter(all(c(0, 1) %in% presence)) %>%   # Keep groups with both 0 and 1
  filter(presence == 1) %>%                # Keep only the rows where presence is 1
  ungroup()


non_mixed <- KEGG %>%
  anti_join(
    KEGG %>%
      filter(presence %in% c(0, 1)) %>%
      group_by(genome, gene_label) %>%
      filter(all(c(0, 1) %in% presence)) %>%
      ungroup(),
    by = c("genome", "gene_label", "presence")
  )

# Combine cleaned + untouched data
final_KEGG <- bind_rows(non_mixed, cleaned_KEGG)


tax<-readRDS("data/tax_MAGs_gtdb.rds")
drep<-read_lines("data/drep_list.txt")


KEGG <- left_join(final_KEGG, tax, by =c('genome'='user_genome')) %>%
  mutate(label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))%>%
  mutate(drep=if_else(genome %in% drep, "yes", "no"))%>
  mutate(drep = if_else(grepl("GCA|GCF", genome), "yes", drep))

## ---- Load and process hydrogenase annotations (hydDB) ----
# hydDB provides additional subtyping of hydrogenase (K00436, K23549, K06281) from hidden Markov models
hydDB<-vroom("data/hyddb-results_24_09_19_v2.csv", delim=";", col_names = c("sequence", "db_hit"))

hydDB$KO <- str_extract(hydDB$sequence, "ko:K\\d+")
hydDB$genome <- str_extract(hydDB$sequence, "^[^_]+(?:_[^_]+)*(?=_(k127|contig|genomic|NZ|NC))")

## Extract, filter, and clean hydrogenase KO hits
# Keep only KOs for Group 1, Group 2, and Group 3d hydrogenases
hydDB<-hydDB%>%
  mutate(KO=gsub("ko:", "", KO))%>%
  mutate(genome=gsub("_genomic", "", genome))%>%
  filter(KO %in% c("K00436", "K23549", "K06281"))%>%
  filter(!db_hit %in% c("NONHYDROGENASE"))%>%
  select(!sequence)%>%
  distinct()

KEGG2 <- KEGG %>%
  distinct()%>%
  left_join(hydDB, by = c("KO", "genome"))%>%
  mutate(Type=if_else(is.na(db_hit), Type, db_hit))%>%
  select(!db_hit)

## Assign hydrogenase subtypes based on gene name patterns
# hox genes → [NiFe] Group 3d, hya genes → [NiFe] Group 1, hupU/V → [NiFe] Group 2
KEGG2<-KEGG2%>%
  mutate(Type=if_else(grepl("hox", Gene), "[NiFe] Group 3d", Type))%>%
  mutate(Type=if_else(grepl("hya", Gene), "[NiFe] Group 1", Type))%>%
  mutate(Type=if_else(grepl("hupU|hupV", Gene), "[NiFe] Group 2", Type))

## ---- Load and process sulfur cycling genes (DsrMKJOP) ----
# DsrMKJOP genes involved in dissimilatory sulfur reduction; map db_hit to standardized KO codes
DsrMKJOP<-vroom("data/dsrMKJOPe10.tsv", delim="\t", col_names = c("sequence", "db_hit"))

DsrMKJOP$genome <- str_extract(DsrMKJOP$sequence, "^[^_]+(?:_[^_]+)*(?=_(k127|contig|genomic|NZ|NC))")

## Map DsrMKJOP protein identifiers to standardized KO codes, clean names
DsrMKJOP<-DsrMKJOP%>%
  select(!sequence)%>%
  mutate(KO = case_when(
    db_hit == "DsrM" ~ "K27187",
    db_hit == "DsrK" ~ "K27188",
    db_hit == "DsrJ" ~ "K27189",
    db_hit == "DsrO" ~ "K27190",
    db_hit == "DsrP" ~ "K27191",
    TRUE ~ "not"  # default case when none of the conditions match
  ))%>%
  mutate(db_hit=gsub("Dsr", "dsr", db_hit))%>%
  mutate(Gene=db_hit)%>%
  distinct()

## Merge DsrMKJOP hits into main KEGG annotation table
KEGG2 <- KEGG2 %>%
  left_join(DsrMKJOP, by = c("genome", "KO", "Gene"))%>%
  mutate(presence=if_else(is.na(db_hit), presence, 1))%>%
  select(!db_hit)

## ---- Export final KEGG annotation table ----
# KEGG2 contains all KEGG annotations with metabolic classifications and custom HMM assignments
saveRDS(KEGG2, "output/KEGG_25_12_02.rds")

## Verify methanotroph presence in dataset
# Filter to genomes with KEGG Module C1 (methane oxidation core) and mmoX/pmoA genes
count<-KEGG%>%filter(Metabolic_step=="Custom\nHMM")%>%filter(presence>0)%>%
  filter(grepl("MFD|LIB", genome))%>%select(genome, Species, Family, Genus, Gene_collapsed, drep)%>%
  filter(!grepl("Burkholderiaceae|Homologous*|Binatales_JAKAVN01|Nitrosomonadaceae|Nevskiales|Mycobacterium|Nevskiales_Macondimonas_put_pmoA", Gene_collapsed))%>%
  filter(!grepl("f__Acetobacteraceae|f__Acidobacteriaceae|f__Nitrosomonadaceae|f__Gallionellaceae|f__Nitrospiraceae", Family))%>%
  select(genome, Species, Genus, Family, drep)%>%
  distinct()

count%>%count(Family)
count%>%count(drep)
count%>%count(Species, drep)
length(count$genome)


