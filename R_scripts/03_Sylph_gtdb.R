#!/usr/bin/env Rscript

########## Not needed to run, just documentation #############
#### Obs, sylph file too large for git repo, can be found at Zenodo repo ####
### parsing of output from Sylph genome quantification with GTDB taxonomy ###
## The output file sylph_gtdb_25_03_06 (uploaded to data folder) can be used directly for plotting genome-level abundances

library(ggplot2)

#!/usr/bin/env Rscript

## 03_Sylph_gtdb.R
## Purpose: link Sylph genome quantification output with GTDB taxonomy,
##          resolve naming issues, prepare genome labels and export a
##          joined dataframe used for plotting genome-level abundance.

library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)
library(ggh4x)
library(treeio)
library(ggtree)


## NOTE: update to your repository root if required
setwd("path/to/your/repo/MFD_methanotrophs_DK/")


##### Load Sylph selection and GTDB link files #####
## `taxonomy_filter_r220.txt` contains lineage fragments used to subset the GTDB output
tax_filter<-vroom("data/taxonomy_filter_r220.txt", delim = "\t", col_names = "tax")%>%
  mutate(tax=gsub("o__Methylococcales;f__UBA1147;g__UBA1147", "g__UBA1147", tax))%>%pull(tax)

## Load renamed taxonomy mapping for user genomes/samples
tax<-readRDS("data/MFD_renamed_tax_25_03_04.rds")


## Identify taxa with name conflicts (duplicate Species entries) to resolve later
tax_issues<-tax%>%group_by(Species)%>%
  summarise(count = n())%>%
  ungroup()

resolved<-tax%>%left_join(tax_issues)%>%filter(count>1)
resolved2<-tax%>%left_join(tax_issues)%>%mutate(user_genome=as.character(user_genome))%>%filter(count>2)%>%pull(user_genome)

## Build trial tax sets: prefer MAG representatives where duplicates exist
tax_try<-tax%>%filter(!user_genome %in% resolved$user_genome)
tax_try2<-tax%>%filter(user_genome %in% resolved$user_genome)%>%
  filter(MAG_Flag=="MAG_")%>%rbind(tax_try)%>%
  filter(!str_detect(user_genome, str_c(resolved2, collapse = "|")))


## ---- Read Sylph GTDB quantification table and subset to taxa of interest ----
# The path below is the Sylph quant output; filter to species-level lines matching `tax_filter`.
MFD_gtdbtest<-vroom("MFD_drep_gtdb_tax_relative_abundance.tsv", delim = "\t")%>% ############# OBS, file too big for git repo, can be found at Zenodo repo
  filter(str_detect(clade_name, str_c(tax_filter, collapse = "|")))%>%
  filter(grepl("s__", clade_name))%>%
  filter(!(grepl("t__", clade_name) & !str_detect(clade_name, str_c(resolved2, collapse = "|"))))


## Fix a specific pair of related entries by summing their numeric columns and
## replacing with a single corrected clade entry (small manual fix for known issue)
cella_fix<-MFD_gtdbtest %>%
  filter(grepl("GCA_003162995.1|MFD02134.bin.3.45", clade_name)) %>%
  summarise(across(where(is.numeric), ~ sum(. , na.rm = TRUE))) %>%
  mutate(clade_name="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995_t__MFD02134.bin.3.45.fa")%>%
  relocate(clade_name)
  
## Remove the broken two entries, add the fixed row, and harmonise a specific species label
MFD_gtdbtest2<-MFD_gtdbtest%>%
  filter(!grepl("GCA_003162995.1|MFD02134.bin.3.45", clade_name)) %>%
  rbind(cella_fix) %>%
  mutate(clade_name=if_else(clade_name=="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995|t__MFD05580.bin.2.10.fa",
                 "d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995_t__MFD05580.bin.2.10.fa", clade_name))%>%
  filter(!clade_name=="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995|t__MFD05580.bin.2.10.fa")%>%  
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")



## Split the GTDB clade string into taxonomic columns and clean SeqId
MFD_gtdb<-MFD_gtdbtest2%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|")%>%
  distinct()


## For taxa with special resolution needs, patch the `tax` table and create `tax_try3`
tax_try3<-tax%>%filter(user_genome %in% resolved2)%>%
  filter(!user_genome=="GCA_003162995.1")%>%
  mutate(Species=if_else(user_genome=="MFD02134.bin.3.45", "s__Methylocella sp003162995_t__MFD02134.bin.3.45.fa",
"s__Methylocella sp003162995_t__MFD05580.bin.2.10.fa"))%>%
  rbind(tax_try2)
  
  
## Join the curated tax mapping with GTDB quant entries by Species
MFD_gtdb_linked<-tax_try3%>%select(Species, user_genome, MAG_Flag, genome_number, label_2, label_3)%>%
  merge(MFD_gtdb, by="Species")  


## If a prepared selection exists, load it and proceed to label/genome ordering
MFD_gtdb<-readRDS("dataframes/sylph_selection_25_03_06.rds")
#### Getting the genomes of interest ####


##########################################################################
## Prepare genome list and phylogenetic ordering for plotting ###########
##########################################################################

genomes<-as.data.frame(tax_try3$user_genome)

## Read an existing tree of representative genomes and clean tip labels
tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/sp_reps_sylph_mfd_25_03_06/MSA.faa.treefile")
tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

## Root tree at midpoint and extract plotting order from ggtree
tr <- phytools::midpoint.root(tree_drep)
pr <- ggtree(tr)
pr

order<-pr$data%>%
  filter(isTip==T)%>%
  arrange(desc(y)) %>%
  pull(label)

## Check for genomes present in GTDB link but missing from tree order
setdiff(unique(MFD_gtdb$user_genome), order)


## Build a small `tax` mapping for labels and genome numbering using the tree order
tax <- MFD_gtdb%>%select(!c(SeqId, Taxonomic_abundance))%>%distinct()%>%
  filter(!is.na(user_genome))%>%
  mutate(MAG_Flag_2 = if_else(grepl("LIB|MFD", user_genome), "MAG_", "")) 

tax <- tax %>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(Genus) %>%
  mutate(genome_number = row_number()) %>%
  ungroup() %>%
  mutate(label_2_new = case_when(
    str_detect(Genus, "drep") ~ Family,
    str_detect(Species, "drep") ~ gsub("g__", "", Genus),
    TRUE ~ gsub("s__", "", Species)))%>%
      mutate(label_3_new = if_else(grepl("MAG", MAG_Flag_2), 
                           paste0(label_2, " ", MAG_Flag_2, genome_number), 
                           label_2_new))%>%
  select(label_2_new, label_3_new, user_genome)%>%
  mutate(user_genome=as.character(user_genome))



## Join GTDB quant table to label mapping and add sample metadata
MFD_gtdb_link <- left_join(MFD_gtdb, tax) %>%
  distinct()

meta<-readRDS("./output/MFD_all_MAGs.rds")%>%
  select(SeqId, fieldsample_barcode, starts_with("mfd_"))%>%
  distinct()

df_full<-MFD_gtdb_link%>%left_join(meta)%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome)

saveRDS(df_full, "output/sylph_gtdb_25_03_06.rds")




