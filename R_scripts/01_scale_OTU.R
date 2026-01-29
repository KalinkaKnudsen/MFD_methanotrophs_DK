#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)


setwd("path/to/your/repo/MFD_methanotrophs_DK/")


habitat <-read_excel("data/2025-02-19_mfd_db.xlsx")
seq_meta <- vroom("data/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)
com_meta<-merge(seq_meta, habitat, by="fieldsample_barcode")%>%
  relocate(SeqId)



# Merging the names and files to fit;

mmoX <- vroom("data/combined_count_table_e10.txt", delim = "\t") %>%
rename(OTU = ConsensusLineage)  %>%
mutate(OTU = if_else(OTU == "Root", "Root_mmoX", OTU))


mmoX_scale<-com_meta%>%
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>%
  mutate(kb_length=(527*3)/1000)%>% #length of seearch mmoX = 527
  select(SeqId, per_million_scale, kb_length)

otutable_mmoX_transposed <- mmoX %>%
  select(-OTU) %>%
  t() %>%
  as.data.frame()

otutable_mmoX_transposed<-otutable_mmoX_transposed%>%
  mutate(SeqId=rownames(otutable_mmoX_transposed))%>%
  relocate(SeqId)
row.names(otutable_mmoX_transposed)=NULL
colnames(otutable_mmoX_transposed)[2:ncol(otutable_mmoX_transposed)] <- mmoX$OTU

# Merge the transposed otutable dataframe with the mmoX_scale dataframe based on SeqId
merged_df_mmoX <- merge(mmoX_scale, otutable_mmoX_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_mmoX <- merged_df_mmoX %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)



##### Loading the pmoA data #####

pmoA <- vroom("data/combined_count_table_e10.txt", delim = "\t") %>%
  rename(OTU = ConsensusLineage)  %>%
  mutate(OTU = if_else(OTU == "Root", "Root_pmoA", OTU))


pmoA_scale<-com_meta%>%
  mutate(per_million_scale=(after_total_reads/2)/1000000)%>%
  mutate(kb_length=(309*3)/1000)%>% #length of seearch pmoA = 309
  select(SeqId, per_million_scale, kb_length)

otutable_pmoA_transposed <- pmoA %>%
  select(-OTU) %>%
  t() %>%
  as.data.frame()

otutable_pmoA_transposed<-otutable_pmoA_transposed%>%
  mutate(SeqId=rownames(otutable_pmoA_transposed))%>%
  relocate(SeqId)
row.names(otutable_pmoA_transposed)=NULL
colnames(otutable_pmoA_transposed)[2:ncol(otutable_pmoA_transposed)] <- pmoA$OTU

# Merge the transposed otutable dataframe with the pmoA_scale dataframe based on SeqId
merged_df_pmoA <- merge(pmoA_scale, otutable_pmoA_transposed, by = "SeqId", all.x = T)%>%
  replace(is.na(.), 0)

# Divide each column in the merged dataframe by the corresponding value in per_million_scale
scaled_otu_table_pmoA <- merged_df_pmoA %>%
  mutate(across(starts_with("Root"), ~./(per_million_scale*kb_length))) %>%
  select(-per_million_scale, -kb_length)



### Merging the dataframes

otu_meta_mmoX<-merge(scaled_otu_table_mmoX, com_meta, by="SeqId")
otu_merged<-merge(otu_meta_mmoX, scaled_otu_table_pmoA, by="SeqId")


#saveRDS(otu_merged, file = "output/merged_otu_table_24_06_24.rds")

########## Now filtering on samples with at least 500.000 reads after processing #####
png("filtering_of_reads.png", width=3000, height=1600, res=350)
jpeg("filtering_of_reads.jpeg", width=3000, height=1600, res=350)
hist(com_meta$after_total_reads, breaks=1500, xlim=c(0,50000000), 
     xlab="Total reads after filtering", main="Histogram of reads per sample")
abline(v=500000, col="red3", lwd=2)
# Add a label for the red line
text(500000, 150, "cutoff: 500 000 reads", col="grey20", pos=4, cex=0.8)
dev.off()


too_few_reads<-com_meta%>%filter(after_total_reads<500000)%>%pull(SeqId)


otu_merged<-otu_merged%>%filter(!SeqId %in% too_few_reads)


#### Cleaning up the tax. 


#otu_merged<-readRDS("output/merged_otu_table_24_06_24.rds")

otu_filtered <- otu_merged %>%
  select(
    matches("^Root; Likely_mmoX"),
    matches("^Root; TUSC"),
    matches("^Root; Putative_pmoA_Binatales"),
    matches("^Root; o_Rhizobiales_pmoA"),
    matches("^Root; pxmA"),
    matches("^Root; o_Methylococcales_pmoA"),
    matches("^Root; Nevskiales_Macondimonas_put_pmoA"),
    matches("^Root; USCg"),
    matches("^Root; Methylomirabilis_pmoA"),
    matches("^Root; Methylacidiphilaceae"),
    matches("Seq|mfd")
  )

# OTU_filtered_long<-pivot_longer(otu_filtered, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
#   mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
#   mutate(type=if_else(grepl("pxmA", Tax), "pxmA", "pmoA"))%>%
#   mutate(type=if_else(grepl("mmoX", Tax), "mmoX", type))

#saveRDS(OTU_filtered_long, "OTU_filtered_long_24_03_15.rds")
#otu_filtered<-readRDS("OTU_filtered_long_24_03_15.rds")

otu_filtered2 <- otu_filtered %>%
  mutate("Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; Methylomonas" = rowSums(select(., matches("Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; Methylomonas_")), na.rm = TRUE)) %>%
  mutate("Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; UBA10906" = rowSums(select(., matches("UBA10906_")), na.rm = TRUE)) %>%
  mutate("Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella" = rowSums(select(., matches("Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_")), na.rm = TRUE)) %>%
  mutate("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylosinus" = rowSums(select(., matches("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylosinus_")), na.rm = TRUE)) %>%
  mutate("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylocystis"=rowSums(select(., matches("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylocystis_")), na.rm = TRUE)) %>%
  mutate("Root; Likely_mmoX; Rhizobiales_mmoX; Rhodomicrobium" = rowSums(select(., matches("Root; Likely_mmoX; Rhizobiales_mmoX; Rhodomicrobium_umbrella")), na.rm = TRUE)) %>%
  mutate("Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylobacter" = rowSums(select(., matches("Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylobacter_")), na.rm = TRUE)) %>%
  mutate("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1; Methylocystis_pmoA1" = rowSums(select(., matches("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1; Methylocystis_pmoA1_")), na.rm = TRUE)) %>%
  mutate("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1; Methylosinus_pmoA1" = rowSums(select(., matches("Methylosinus_pmoA1_")), na.rm = TRUE)) %>%
  mutate("Root; o_Rhizobiales_pmoA; Methylosinus_Methylocystis_pmoA2; Methylocystis_pmoA2" = rowSums(select(., matches("Root; o_Rhizobiales_pmoA; Methylosinus_Methylocystis_pmoA2; Methylocystis_pmoA2_")), na.rm = TRUE)) %>%
  mutate("Root; o_Rhizobiales_pmoA; Rhodomicrobium_pmoA" = rowSums(select(., matches("Root; o_Rhizobiales_pmoA; Rhodomicrobium_pmoA")), na.rm = TRUE)) %>%  
  ## Specifically selecting the two clades beneath to aggregate the "Singleton" clade
  mutate("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1" = rowSums(select(.,c("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1; PmoA1_Singleton_clade", "Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1")), na.rm = TRUE)) %>%
  rename("Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylomonas"="Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylomonas_1")%>%
  rename("Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylovulum"="Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylovulum_2")%>%
  
  select(-matches("Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; Methylomonas_"), -matches("UBA10906_"), -matches("Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_"),
         -matches("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylosinus_"),
         -matches("Root; Likely_mmoX; Rhizobiales_mmoX; Methylocystis_Methylosinus; Methylocystis_"),         
         -matches("Root; Likely_mmoX; Rhizobiales_mmoX; Rhodomicrobium_umbrella"), 
         -matches("Root; o_Methylococcales_pmoA; Methylomonadaceae; Methylobacter_"), 
         -matches("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1; Methylocystis_pmoA1_"), 
         -matches("Methylosinus_pmoA1_"), 
         -matches("Root; o_Rhizobiales_pmoA; Methylosinus_Methylocystis_pmoA2; Methylocystis_pmoA2_"),
         -matches("Root; o_Rhizobiales_pmoA; Rhodomicrobium_pmoA; Rhodomicrobium_"),
         -matches("Methylosinus_Methylocystis_pmoA1; PmoA1_Singleton_clade"))


test3<-otu_filtered2%>%
  select(matches("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1|SeqId"))



##### After inspection of OTUs that have no representation and are not valuable to show, I am removing the following:

#### Now to collapsing some groups.
tax_remove=c("Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; UBA10906",
             "Root; Likely_mmoX; Methylococcales_mmoX; JABFRO01", 
             "Root; Likely_mmoX; Methylococcales_mmoX; Methylomonadaceae_mmoX; JABFRC01",
             "Root; o_Methylococcales_pmoA",
             "Root; pxmA; Methylumidiphilus_pxmA",
             "Root; o_Methylococcales_pmoA; Methylomonadaceae; CAIQWF01")


OTU_filtered_long<-pivot_longer(otu_filtered2, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  filter(!Tax %in% tax_remove)%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  mutate(type=if_else(grepl("pxmA", Tax), "pxmA", "pmoA"))%>%
  mutate(type=if_else(grepl("mmoX", Tax), "mmoX", type))%>%
  mutate(type=if_else(grepl("pmoA2", Tax), "pmoA2", type))%>%
  mutate(type=if_else(grepl("Nevskiales_Macondimonas|Binatales|Rhodopila", Tax), "Put. pmoA/mmoX", type))

sort(unique(OTU_filtered_long$Tax))




saveRDS(OTU_filtered_long, "./output/OTU_filtered_long_25_09_01.rds")


