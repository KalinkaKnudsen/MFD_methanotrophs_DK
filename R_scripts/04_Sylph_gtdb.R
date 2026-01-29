#!/usr/bin/env Rscript


library(ggplot2)
library(vroom)
library(tidyverse)
library(readxl)
library(ggh4x)
library(treeio)
library(ggtree)


setwd("path/to/your/repo/MFD_methanotrophs_DK/")


##### Loading the sylph file


tax_filter<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/taxonomy_filter_r220.txt", delim = "\t", col_names = "tax")%>%
  mutate(tax=gsub("o__Methylococcales;f__UBA1147;g__UBA1147", "g__UBA1147", tax))%>%pull(tax)

tax<-readRDS("MFD_renamed_tax_25_03_04.rds")

tax_issues<-tax%>%group_by(Species)%>%
  summarise(count = n())%>%
  ungroup()

resolved<-tax%>%left_join(tax_issues)%>%filter(count>1)
resolved2<-tax%>%left_join(tax_issues)%>%mutate(user_genome=as.character(user_genome))%>%filter(count>2)%>%pull(user_genome)


tax_try<-tax%>%filter(!user_genome %in% resolved$user_genome)
tax_try2<-tax%>%filter(user_genome %in% resolved$user_genome)%>%
  filter(MAG_Flag=="MAG_")%>%rbind(tax_try)%>%
  filter(!str_detect(user_genome, str_c(resolved2, collapse = "|")))

MFD_gtdbtest<-vroom("/home/bio.aau.dk/wz65bi/mfd_sylph_quant/analysis/MFD_drep_gtdb_quant/MFD_drep_gtdb_tax_relative_abundance.tsv", delim = "\t")%>%
  #filter(grepl("g__Methylobacter_C", clade_name))
  filter(str_detect(clade_name, str_c(tax_filter, collapse = "|")))%>%
  filter(grepl("s__", clade_name))%>%
  filter(!(grepl("t__", clade_name) & !str_detect(clade_name, str_c(resolved2, collapse = "|"))))


cella_fix<-MFD_gtdbtest %>%
  filter(grepl("GCA_003162995.1|MFD02134.bin.3.45", clade_name)) %>%
  summarise(across(where(is.numeric), ~ sum(. , na.rm = TRUE))) %>%
  mutate(clade_name="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995_t__MFD02134.bin.3.45.fa")%>%
  relocate(clade_name)
  
MFD_gtdbtest2<-MFD_gtdbtest%>%
  filter(!grepl("GCA_003162995.1|MFD02134.bin.3.45", clade_name)) %>%
  rbind(cella_fix) %>%
  mutate(clade_name=if_else(clade_name=="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995|t__MFD05580.bin.2.10.fa",
                 "d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995_t__MFD05580.bin.2.10.fa", clade_name))%>%
  filter(!clade_name=="d__Bacteria|p__Pseudomonadota|c__Alphaproteobacteria|o__Rhizobiales|f__Beijerinckiaceae|g__Methylocella|s__Methylocella sp003162995|t__MFD05580.bin.2.10.fa")%>%  
  pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")



MFD_gtdb<-MFD_gtdbtest2%>%
  mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
  separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|")%>%
  distinct() #%>%
  #mutate(user_genome=gsub(".fa", "", user_genome))%>%
  #mutate(user_genome=gsub("t__", "", user_genome))



tax_try3<-tax%>%filter(user_genome %in% resolved2)%>%
  filter(!user_genome=="GCA_003162995.1")%>%
  mutate(Species=if_else(user_genome=="MFD02134.bin.3.45", "s__Methylocella sp003162995_t__MFD02134.bin.3.45.fa",
"s__Methylocella sp003162995_t__MFD05580.bin.2.10.fa"))%>%
  rbind(tax_try2)
  
  
  
MFD_gtdb_linked<-tax_try3%>%select(Species, user_genome, MAG_Flag, genome_number, label_2, label_3)%>%
  merge(MFD_gtdb, by="Species")  


#MFD_gtdb_linked<-MFD_gtdb%>%left_join(tax_try3, by="Species")
#### From ANI, it is evident that GCA_003162995.1 should be aggregated with MFD02134.bin.3.45
#saveRDS(MFD_gtdb_linked, "dataframes/sylph_selection_25_03_06.rds")

MFD_gtdb<-readRDS("dataframes/sylph_selection_25_03_06.rds")
#### Getting the genomes of interest ####


##########################################################################
## Exporting genomes for a genome tree which has the genomes from Sylph ##
##########################################################################

genomes<-as.data.frame(tax_try3$user_genome)

#readr::write_csv(genomes, "/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/sp_reps_sylph_mfd_25_03_06/genomes.txt", col_names = F)


#tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/no_gtdb/MSA.faa.treefile")
tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/sp_reps_sylph_mfd_25_03_06/MSA.faa.treefile")
# 
 tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
 tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

tr <- phytools::midpoint.root(tree_drep)
# plot rooted tree
pr <- ggtree(tr) 
pr

order<-pr$data%>%
  filter(isTip==T)%>%
  #filter(grepl("LIB|MFD", label))%>%
  arrange(desc(y)) %>%
  pull(label)

setdiff(unique(MFD_gtdb$user_genome), order)

tax <- MFD_gtdb%>%select(!c(SeqId, Taxonomic_abundance))%>%distinct()%>%
  filter(!is.na(user_genome))%>%
#  group_by(Species)%>%
  mutate(MAG_Flag_2 = if_else(grepl("LIB|MFD", user_genome), "MAG_", "")) 

#tax <- MFD_gtdb%>%select(!c(SeqId, Taxonomic_abundance))%>%distinct()%>%left_join(linkage)


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



# tax <- tax %>%
#   mutate(
#     label = case_when(
#       str_detect(Genus, "drep") ~ paste0(Family, ";", Genus),
#       str_detect(Species, "drep") ~ paste0(Genus, ";", Species),
#       TRUE ~ Species
#     ),
#     label_2 = case_when(
#       str_detect(Genus, "drep") ~ Family,
#       str_detect(Species, "drep") ~ gsub("g__", "", Genus),
#       TRUE ~ gsub("s__", "", Species)
#     ) %>% paste0(" MAG_")
#   )%>%
#   mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
#   arrange(user_genome) %>%
#   group_by(Genus) %>%
#   mutate(genome_number = row_number()) %>%
#   ungroup() %>%
#   mutate(label_3 = paste0(label_2, genome_number))%>%
#   select(!label_2)

#saveRDS(tax, "MFD_renamed_tax_25_03_06.rds")



# Display the modified dataframe

  
MFD_gtdb_link <- left_join(MFD_gtdb, tax) %>%
  #  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  distinct()




meta<-readRDS("./output/genome_quantification/MFD_all_MAGs.rds")%>%
  select(SeqId, fieldsample_barcode, starts_with("mfd_"))%>%
  distinct()

df_full<-MFD_gtdb_link%>%left_join(meta)%>%
mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome)        

saveRDS(df_full, "output/genome_quantification/sylph_gtdb_25_03_06.rds")


