#!/usr/bin/env Rscript


library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(ggtree)
library(ape)
library(ggh4x)
library(viridis)

setwd("path/to/your/repo/MFD_methanotrophs_DK/")


KEGG<-readRDS("./output/metabolism/KEGG_25_12_02.rds")%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))


methane_colors<-readRDS("palette_methane_colors.rds")%>%
  unlist()
names(methane_colors)[names(methane_colors) == "-"] <- "Formate oxidation"


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

KEGG_sort<-KEGG%>%
  arrange(Genus)

KEGG<-KEGG%>%
  mutate (gene_label = paste0(Gene_collapsed , ' ',  Pathway_step))


####################################################################################################################################
###################### From here, we are working with all the genomes from GTDB that I have investigated ###########################
####################################################################################################################################


KEGG<-KEGG%>%
  mutate(class_label = paste0(Class, ";", Family,";", tree_label))


genes_remove_string<-c("Nitrosomonadaceae (1/3)", "Homologous_MO (1/3)", "Nitrosomonas (1/3)", "Propane_MO_Actino_cluster (1/3)", 
                       "Homologous_pmoA (1/3)", "Cycloclasticus (1/3)", "Betaproteobacteria_amoA (1/3)",
                       "Nitrosococcus (1/3)", "Nitrospira_clade_B (1/3)", "Actinobacteria (1/3)", "Homologous_Rhodopila (1/3)")

genus_remove_string<-c("g__Methylocella", "g__Methylomonas", "g__Methylocystis", "g__Methylomirabilis", "g__Methylosinus", "g__Methylobacter_C",
                       "g__Methylovulum", "g__Methylobacter_A", "g__Crenothrix", "g__Methylosarcina", "g__Methylobacter", "g__Methyloferula", "g__Methylocaldum", "g__Nitrospira_D",
                       "g__Bradyrhizobium", "g__Nitrosospira", "g__Methylobacter_B", "g__Nitrosomonas", "g__Methyloglobulus", "g__Methylomicrobium", "g__Methylotuvimicrobium", 
                       "g__Methylicorpusculum", "g__Methylomarinum", "g__Methylococcus", "g__Methylohalobius", "g__Methylacidimicrobium", "g__Methylacidiphilum", "g__Methylacidithermus",
                       "g__Methyloprofundus")

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


genes_keep<-genomes_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()

# pl <- KEGG %>% 
#   filter(genome %in% genomes_keep$genome)%>%
#   filter(Module == 'C1') %>% 
#   filter(!Metabolic_step=="Custom\nHMM")%>%
#   #filter(type %in% c('LR-MAG', 'GTDB'))%>%
#   mutate(presence=if_else(presence==0, NA, presence))%>%
#   mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
#   #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
#   mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
#   ggplot(aes(x = gene_label, y = class_label))+ 
#   geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
#   scale_fill_manual(na.value="transparent", values=c(
#     "pMMO"="#b3943c",
#     "sMMO"="darkgreen",
#     "calcium (mxa)"="#5e4fa2",
#     "lanthanide (xoxF)"="#c46ca1",
#     "putative"="grey75",
#     "H4MPT"="#d97512",
#     "GSH"="turquoise4",
#     "H4F"="darkred",
#     "Serine cycle"="darkred",
#     "RuMP cycle"="darkblue",
#     "CBB cycle"="#679b60",
#     "Formate oxidation"="#679b60"
#   )) +
#   facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
#         axis.text.y = element_blank(),
#         axis.ticks.x = element_line(linewidth = 0.1),
#         axis.ticks.y = element_blank(),
#         strip.clip = "off",
#         strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
#         strip.text.y = element_blank(), 
#         strip.background = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "none",
#         panel.spacing = unit(0.03, "cm", data = NULL),
#         legend.title=element_blank(),
#         panel.background = element_rect(fill="white")
#   )
# 
# pl
# 
# 
# CM <- KEGG %>% 
#   filter(genome %in% genomes_keep$genome)%>%
#   filter(Module == 'C1') %>% 
#   filter(Metabolic_step=="Custom\nHMM")%>%
#   filter(gene_label %in% genes_keep$gene_label)%>%
#   mutate(presence=if_else(presence==0, NA, presence))%>%
# #  filter(presence>0)%>%
#   mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
#   mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
#   ggplot(aes(x = gene_label, y = class_label))+ 
#   geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.15 ) + 
#   scale_fill_manual(values=c(methane_colors), na.value = "white") +
#   facet_nested(. ~ Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
#         axis.text.y = element_text(size=8),
#         axis.ticks = element_line(linewidth = 0.1),
#         axis.title = element_blank(),
#         legend.position = "none",
#         strip.clip = "off",
#         strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
#         strip.text.y=element_blank(),
#         strip.background = element_blank(),
#         text = element_text(family = "Arial"),
#         panel.spacing = unit(0.03, "cm", data = NULL),
#         legend.title=element_blank(),
#         panel.background = element_rect(fill="white"))+
#   scale_x_discrete(c(0,0))+
#   scale_y_discrete(c(0,0))
# 
# CM
# 
# options("aplot_guides" = "keep")
# met<-pl%>%aplot::insert_left(CM, width=0.32)
# #met
# 
# 
# ggsave("./output/metabolism/gtdb_met_25_08_11.png",met, height=20, width=22, dpi=400)


unique(CM$data$genome)



# 
# missing_DRAM<-setdiff(CM$data%>%filter(presence==1)%>%select(genome)%>%distinct(), 
#                       pl$data%>%filter(!Metabolic_step=="Custom\nHMM")%>%filter(presence==1)%>%
#                         select(genome)%>%distinct())
# 
# x<-KEGG%>%filter(genome %in% missing_DRAM$genome)%>%select(class_label, Genus)%>%distinct()
# 

#readr::write_csv(missing_DRAM, "~/data/MFD/annotation_combined/DRAM_25_03_21_gtdb/missing_genomes.txt", col_names = FALSE)



# Okay now, I can start with seperation into TUSC+Bin and the rest 



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


pl_bin_tusc <- KEGG %>% 
  filter(genome %in% genomes_Bin_keep$genome)%>%
  mutate(fam_label=paste0(Family, " ", Species))%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  #filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
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
#ggsave("./output/metabolism/gtdb_Bin_TUSC_met_25_08_13.png",pl_bin_tusc, height=6, width=14, dpi=400)

ggsave("./output/metabolism_Extended_figs/GTDB_Bin_TUSC_25_12_10.png",pl_bin_tusc,
       units = c("mm"),
       height = 90,
       width = 190,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/GTDB_Bin_TUSC_25_12_10.svg",pl_bin_tusc,
       units = c("mm"),
       height = 90,
       width = 190,
       dpi=300)





############ Great, that works nicely #########

### Now, I want to have a look at the genomes that are somewhat doubtfull ##
## First round, I do not want any of the genomes from the "Methylomonadaceae" or the "Methylococcaceae"


pl <- KEGG %>% 
  filter(genome %in% genomes_keep$genome)%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  filter(!Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae"))%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
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
  facet_nested(.~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
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




############# Instead of patching two, I can do like so: ###


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
  #filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  #mutate(gene_label=gsub(" (1/3)", "", gene_label))%>%
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
  facet_nested(Class~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
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

#pl_2


#ggsave("./output/metabolism/gtdb_sub_met_25_08_13.png",pl_2, height=13, width=16, dpi=400)

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




genes_keep<-genomes_keep %>% 
  filter(!grepl("MFD|LIB", genome))%>%
  filter(!genome %in% genomes_Bin_keep$genome)%>%
  # filter(Family %in% c("f__Methylomonadaceae", "f__Methylococcaceae"))%>%
  filter(Order %in% c("o__Methylococcales") | Family %in% c("f__JACCXJ01"))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  # mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  distinct()

genes_remove<-KEGG%>%  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(!gene_label %in% genes_keep$gene_label)%>%
  select(gene_label)%>%distinct()


print(n=26, KEGG%>%filter(gene_label %in% genes_keep$gene_label)%>%
        select(KO, gene_label)%>%distinct())



  #mutate(gene_label=if_else(grepl("Methylotenera|JABFRO01", gene_label), paste0("Methylococcales_like_mmoX (1/3)"), gene_label))%>%

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
  #filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
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


ggsave("./output/metabolism_Extended_figs/GTDB_Methylococcales_gtdb_sub_met_25_12_10.png",pl_3,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/GTDB_Methylococcales_gtdb_sub_met_25_12_10.svg",pl_3,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)













##################### Making a quick MFD-TUSC plot ######################


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


## Okay let us have a tiny look at cont and comp. 
mags_v4_low_comp<-vroom("/projects/microflora_danica/MFD-LR/mags_v4/MFD_mags_V4.tsv")%>%
  filter(bin %in% genomes_Bin_keep$genome)%>%
  filter(Completeness_CheckM2<80)

mags_SR_low_comp<-vroom("/projects/microflora_danica/MFD-SR/results/mags_shallow_all.tsv")%>%
  filter(bin %in% genomes_Bin_keep$genome)%>%filter(Completeness_CheckM2<80)

genomes_Bin_keep_filt<-genomes_Bin_keep%>%
  filter(!genome %in% mags_v4_low_comp$bin)%>%
  filter(!genome %in% mags_SR_low_comp$bin)


genes_remove<-KEGG%>%  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(!gene_label %in% genes_Bin_keep$gene_label)%>%
  select(gene_label)%>%distinct()



### counts 
genomes_Bin_keep%>%ungroup()%>%filter(Gene=="TUSC")%>%count(drep)
genomes_Bin_keep%>%ungroup()%>%group_by(Order)%>%filter(!Gene=="TUSC")%>%count(drep)

# 
# drep<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/mfd_all_drep_reps.tsv")%>%
#   rename(genome=bin)
# 
# genomes_Bin_keep2 <- genomes_Bin_keep %>%
#   left_join(drep) %>%
#   # 1) Genus/family text to use in labels (strip g__/f__)
#   mutate(
#     genus_or_family = if_else(
#       is.na(Genus),
#       str_replace_na(Family, "") %>% str_replace("^f__", ""),
#       str_replace(Genus, "^g__", "")
#     )
#   ) %>%
#   # 2) MAG numbering: row number within (genus/family, drep)
#   group_by(genus_or_family, drep) %>%
#   mutate(mag_number = row_number()) %>%
#   ungroup() %>%
#   # 3) Dup numbering: row number within species_rep cluster
#   group_by(species_rep) %>%
#   mutate(dup_number = row_number()) %>%
#   ungroup() %>%
#   # 4) Build a per-row MAG label (used by the representative)
#   mutate(mag_label = paste0(genus_or_family, " MAG_", mag_number)) %>%
#   # 5) Broadcast the representative MAG label within each species_rep group
#   group_by(species_rep) %>%
#   mutate(
#     rep_label_1 = mag_label[which(genome == species_rep)][1],
#     rep_label_1 = coalesce(rep_label_1, first(mag_label)),
#     # 6) Final label:
#     #    - representative gets the MAG label (genus-based numbering)
#     #    - non-reps get "sp. dup. <dup_number>" (species-cluster numbering)
#     rep_label_2 = if_else(
#       genome == species_rep,
#       rep_label_1,
#       paste0(rep_label_1, " sp. dup. ", dup_number)
#     )
#   ) %>%
#   ungroup()%>%select(rep_label_2, genome)




pl_bin_tusc <- KEGG %>% 
  # filter(genome %in% genomes_Bin_keep_filt$genome)%>%
  # left_join(genomes_Bin_keep2)%>%
  mutate(Order=gsub("o__", "", Order))%>% 
  #mutate(fam_label=paste0(Family, " ", Species))%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  #filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
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

ggsave("./output/metabolism_Extended_figs/MFD_Bin_TUSC_sub_met_25_12_10.png",pl_bin_tusc,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

ggsave("./output/metabolism_Extended_figs/MFD_Bin_TUSC_sub_met_25_12_10.svg",pl_bin_tusc,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)




################### Let us make an alternative electrons plot ########################


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


NS

ggsave("./output/metabolism_Extended_figs/MFD_Bin_TUSC_Nitro_sulphur_25_12_10.png",NS,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

ggsave("./output/metabolism_Extended_figs/MFD_Bin_TUSC_Nitro_sulphur_25_12_10.svg",NS,
       units = c("mm"),
       height = 170,
       width = 185,
       dpi=300)

#ggsave("./output/metabolism/TUSC_Binatales_MFD_met_Nitro_sulphur_C_25_08_21.png",NS, height=12, width=18, dpi=400)









#### Okay some DRAM were missing, I will figure out which:
x1<-KEGG %>% 
  filter(genome %in% genomes_Bin_keep$genome)%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(presence==1)%>%
  select(genome)%>%
  distinct()


x2<-KEGG %>% 
  filter(genome %in% genomes_Bin_keep$genome)%>%
  filter(Module == 'C1') %>% 
  filter(!gene_label %in% genes_remove$gene_label)%>%
  filter(!Metabolic_step=="Custom\nHMM")%>%
  filter(presence==1)%>%
  select(genome)%>%
  distinct()

missing<-setdiff(x1, x2)


readr::write_csv(missing, "~/data/MFD/annotation_combined/genomes_missing_25_08_13.txt", col_names = FALSE)



###3 Getting the TUSC genomes:
x3<-genomes_Bin_keep%>%ungroup()%>%filter(Gene=="TUSC")%>%select(genome)
readr::write_csv(x3, "/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/TUSC_genes/genomes.txt", col_names = FALSE)

