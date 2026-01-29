#!/usr/bin/env Rscript


.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(patchwork)

setwd("~/scripts/MFD/methanotrophs/R_scripts/output/")




p_combine<-readRDS("../palette_mfd_hab2_ISME.rds") 
OTU_filtered_long<-readRDS("OTU_filtered_long_24_06_24.rds")



tax<-OTU_filtered_long%>%
  select(Tax)%>%
  distinct()

#tax_remove=c("Root; o_Nitrososphaerales; f_Nitrosopumilaceae; Nitrosotalea_TA20_cluster", "Root; o_Nitrososphaerales; f_Nitrososphaeraceae; TH5893_Nitrososphaera_cluster; JAJPHM01_Nitrososphaera_cluster", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster", "Root; periplasmic_umbrella; Nitrohelix_Nitrospina_Nitrospiraceae_cluster; f_Nitrospinaceae_unknown_cluster" )

tax_curated<-tax%>%
  mutate(Tax_curated=gsub("Root; Likely_mmoX; ", "", Tax), 
         Tax_curated = if_else(Tax=="Root; Likely_mmoX", "Likely_mmoX", Tax_curated),
         Tax_curated=gsub("Root; o_Methylococcales_pmoA; ", "", Tax_curated),
         Tax_curated=gsub("Root; o_Rhizobiales_pmoA; ", "", Tax_curated),
         Tax_curated=gsub("Root; pxmA; ", "", Tax_curated),
         Tax_curated=gsub("Root; ", "", Tax_curated),
         Tax_curated=gsub("Methylococcales_mmoX; Methylomonadaceae_mmoX", "Methylomonadaceae_mmoX", Tax_curated),
         Tax_curated=gsub("Beijerinckiaceae_pmoA1; Methylosinus_Methylocystis_pmoA1;", "Beijerinckiaceae_pmoA1;", Tax_curated),
         Tax_curated=gsub("Methylosinus_Methylocystis_pmoA2;", "", Tax_curated),
         Tax_curated = if_else(Tax_curated=="Methylococcales_mmoX", "Methylococcales_mmoX_umbrella", Tax_curated),
         Tax_curated=gsub("Methylococcales_mmoX; ", "", Tax_curated),
         Tax_curated = if_else(Tax_curated=="Methylomonadaceae_mmoX", "Methylomonadaceae_mmoX_umbrella", Tax_curated),
         Tax_curated = if_else(Tax_curated=="Methyloccocaceae_mmoX", "Methyloccocaceae_mmoX_umbrella", Tax_curated),
         Tax_curated=gsub("Rhizobiales_mmoX; Methylocystis_Methylosinus; ", "Rhizobiales_mmoX; ", Tax_curated),
         Tax_curated = if_else(Tax_curated=="Rhizobiales_mmoX; Methylocystis_Methylosinus", "Rhizobiales_mmoX; Methylocystis_Methylosinus_umbrella", Tax_curated),
         Tax_curated = if_else(Tax_curated=="USCg", "USCg_umbrella", Tax_curated),
         Tax_curated = if_else(Tax_curated=="pxmA", "pxmA_umbrella", Tax_curated),
         Tax_curated=gsub("_pmoA1;", ";", Tax_curated),
         Tax_curated=gsub("_pmoA2;", ";", Tax_curated),
         Tax_curated=gsub("_mmoX;", ";", Tax_curated),
         Tax_curated=gsub("_pxmA", "", Tax_curated))






#Merging the curated tax onto the dataframe:


OTU_filtered_long<-merge(OTU_filtered_long, tax_curated, by= "Tax")
# tax2 <- readxl::read_excel("tax2.xlsx")

OTU_filtered_long<-OTU_filtered_long %>%
  mutate(mfd_hab2=gsub("\\s*\\(non-habitat type\\)\\s*", "", mfd_hab2))%>%
  mutate(across(mfd_hab1, ~str_replace(., "Sclerophyllous scrub", "Temperate heath and scrub")),
         across(mfd_hab2, ~str_replace(., "Scrub", "Sclerophyllous scrub"))) %>%
  mutate(complex_long = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", ")) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = "\n")) %>%
  mutate(mfd_hab2 = if_else(complex == "Soil\nNatural\nForests", mfd_hab3, mfd_hab2))%>%
  mutate(mfd_hab2 = if_else(grepl("Beech", mfd_hab2), "Beech", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Birch", mfd_hab2), "Birch", mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(grepl("Oak", mfd_hab2), "Oak", mfd_hab2)) %>%
  mutate(mfd_hab2= if_else(grepl("Lake", mfd_hab3), "Standing freshwater, lake", mfd_hab2))%>%
  mutate(mfd_hab2 = if_else(complex == "Sediment\nUrban\nFreshwater", paste0(mfd_hab3), paste0(mfd_hab2)))%>%
  mutate(mfd_hab2 = if_else(grepl("Rainwater basin", mfd_hab2), "Rainwater basin", mfd_hab2))%>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="Standing freshwater", "Standing freshwater, other", mfd_hab2))%>%
  filter(!mfd_hab2 %in% c("Enclosed water, Dried", "Birch", "Pine", "Mire"))%>%
  mutate(mfd_hab2 = if_else(grepl("Enclosed water", mfd_hab2), "Urban enclosed water", mfd_hab2))%>% 
  mutate(mfd_hab2 = if_else(is.na(mfd_hab2) & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " no MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="NA" & complex == "Soil\nNatural\nForests", paste0(mfd_hab1, " no MFDO2"), mfd_hab2)) %>%
  mutate(mfd_hab2 = if_else(mfd_hab2=="Semi-natural tall-herb humid meadows", "Semi-natural humid meadows", mfd_hab2))%>%
  filter(!is.na(mfd_hab2))%>%
  filter(!mfd_hab2=="NA")%>%
  mutate(mfd_hab2=gsub("\\s*\\(exotic\\)", "", mfd_hab2))%>%
  distinct()%>%
  mutate(hab1_label=gsub("Bogs, mires and fens", "Bogs, mires\nand fens", mfd_hab1),
         hab1_label=gsub("Grassland formations", "Grassland\nformations", hab1_label),
         hab1_label=gsub("Temperate heath and scrub","Heath and\nscrub", hab1_label),
         hab1_label=gsub("Greenspaces","Green-\nspaces", hab1_label),
         hab1_label=gsub("Coastal","Coast", hab1_label),
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))





seq<-vroom("../2023-10-11_samples_minimal_metadata_collapsed.csv", delim=",")%>%select(fieldsample_barcode, flat_name)
meta<-read_excel("../2024-02-13_mfd_db.xlsx")%>%select(fieldsample_barcode, latitude, longitude, project_id)%>%left_join(seq)%>%rename(SeqId=flat_name)%>%mutate(SeqId=gsub(".fastq.gz", "", SeqId))

OTU_filtered_long<-OTU_filtered_long%>%left_join(meta)

library(mapDK) 
library(gridExtra) #For plotting


colors <- c( "gold3","#33a08c","coral","#1a0060", "darkred", "#F8E622","darkgreen","#6a3d9a","#1f78b4", "firebrick3")


m <- NULL
m <- mapDK(detail = 'region') #Setting the boarders to regions, you can change this if you want


grasslands_calserous_USCa<-OTU_filtered_long%>%filter(mfd_hab3=="Calcareous grassland")%>%filter(Tax_curated=="Beijerinckiaceae; Methylocella; Methylocella_USCa")
grasslands_calserous_USCg<-OTU_filtered_long%>%filter(mfd_hab3=="Calcareous grassland")%>%filter(Tax_curated=="USCg; JACCX01")

MAP<- m+ 
  geom_jitter(data = grasslands_calserous_USCa, aes(x=longitude, y=latitude,
                                alpha=RPKM,
                                color="black",group = fieldsample_barcode),
              size=2.5, stroke = 1.5, height=0.001, width=0.001, shape=21, fill="darkred") +
  theme(legend.position = c(0.8,0.6),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.2))

MAP



MAP<-grasslands_calserous_USCg<- m+ 
  geom_jitter(data = grasslands_calserous_USCg, aes(x=longitude, y=latitude,
                                                    alpha=RPKM,
                                                    color=mfd_hab1,group = fieldsample_barcode),
              size=2.5, stroke = 1.5, height=0.001, width=0.001, shape=21, fill="darkred") +
  theme(legend.position = c(0.8,0.6),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.2))

MAP


MAPa <- m + 
  geom_jitter(data = grasslands_calserous_USCa, aes(x=longitude, y=latitude,
                                                    fill=RPKM,  # Use fill instead of alpha for mapping
                                                    group=fieldsample_barcode),
              size=2.5, stroke=1, height=0.001, width=0.001, 
              shape=21, color="black", alpha=1) +  # Set alpha=1 globally for outlines
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c("#c7c5c510", "#ad232340", "#ad2323", "darkred"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "darkred"  # Handling NA values
  ) +
  theme(legend.position = c(0.8,0.6),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.2))

MAPa



MAPg <- m + 
  geom_jitter(data = grasslands_calserous_USCg, aes(x=longitude, y=latitude,
                                                    fill=RPKM,  # Use fill instead of alpha for mapping
                                                    group=fieldsample_barcode),
              size=2, stroke=1, height=0.03, width=0.03, 
              shape=21, color="black", alpha=0.7) + 
  
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c("#c7c5c510", "#ad232340", "#ad2323", "darkred"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "darkred"  # Handling NA values
  ) +
  theme(legend.position = c(0.8,0.6),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.2))

MAPg




Methylocystis<-sylph_broad%>%left_join(meta)%>%filter(Genus=="g__Methylocystis")
MFD_map_cystis<- m+ 
  geom_jitter(data = Methylocystis, aes(x=longitude, y=latitude, shape=label_3,
                                        #   fill = mfd_hab1,
                                        color=mfd_hab1,group = fieldsample_barcode),
              size=2.5, alpha = 0.9, stroke = 1.8, height=0.03, width=0.03) +
  scale_color_manual(values=c(
    "Bogs, mires and fens"="#b3943c",
    "Forests"="darkgreen",
    "Grassland formations"="#5e4fa2",
    "Temperate heath and scrub"="#c46ca1",
    "Dunes"="#d97512",
    "Freshwater"="turquoise4",
    "Fields"="darkred",
    "Coastal"="darkblue",
    "Rocky habitats and caves"="coral",
    "Greenspaces"="#679b60"))+  
  scale_shape_manual(values=c(21, 22, 23, 24, 25))+
  labs(color = "Habitattype", shape="Species_cluster",title="Methylocystis | Tax. Rel. abund. > 1%")+
  theme(legend.position = c(0.8,0.7),
        legend.background = element_blank(),
        plot.title = element_text(hjust = 0.2))

MFD_map_cystis

