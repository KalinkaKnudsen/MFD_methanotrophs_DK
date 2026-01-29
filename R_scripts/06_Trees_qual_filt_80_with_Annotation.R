#!/usr/bin/env Rscript

#.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
library(ggplot2)
#library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
#library(patchwork)
library(ggtree)
library(ape)
library(ggh4x)
library(treeio)
library(viridis)

setwd("~/scripts/MFD/methanotrophs/R_scripts/")

# Tree file
tree_80 <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/filter_80/tree/MSA.faa.treefile")

tree_80$tip.label <- gsub("RS_", "", tree_80$tip.label)
tree_80$tip.label <- gsub("GB_", "", tree_80$tip.label)


tr <- phytools::midpoint.root(tree_80)
# plot rooted tree
pr <- ggtree(tr) 
pr

order<-pr$data%>%
  filter(isTip==T)%>%
  arrange(desc(y)) %>%
  pull(label)


tax_old<-readRDS("tax_MAGs_gtdb.rds")

gtdb<-vroom("/home/bio.aau.dk/vj52ou/data/gtdb/bac120_taxonomy_r220.tsv", delim="\t", col_names = c("user_genome", "classification"))%>%
  mutate(user_genome=gsub("RS_", "",user_genome),
         user_genome=gsub("GB_", "",user_genome))%>%
  filter(user_genome %in% tree_80$tip.label)

tax<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/MFD_SR_LR_taxonomy_r220.tsv", delim = "\t", col_select =c("user_genome", "classification"))%>%
  filter(user_genome %in% tree_80$tip.label)%>%
  rbind(gtdb)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'Microflora Danica','GTDB'),
         type=if_else(grepl("LIB", user_genome), "Microflora Danica", type))

tax <- tax %>%
  separate(classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ';', remove = FALSE) %>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", user_genome ))

MFD_renamed<-readRDS("MFD_renamed_tax_25_03_04.rds")%>%
  #rename(user_genome=bin)%>%
  select(user_genome, label_3)

drep<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/mfd_all_drep_reps.tsv")%>%
  rename(user_genome=bin)

tax2<-left_join(tax, MFD_renamed)%>%
  mutate(label_3=if_else(type=="GTDB", gsub("s__", "", Species), label_3))%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(Genus) %>%
  mutate(genome_number = row_number()) %>%
  ungroup() %>%
  mutate(label_3=if_else(is.na(label_3), paste0(gsub("g__","",Genus), " sp. dup. ", genome_number), label_3))%>%
  mutate(user_genome=as.character(user_genome))

# tax<-readRDS("MFD_renamed_tax_25_03_04.rds")%>%
#   mutate(type = if_else(grepl("MFD", user_genome), 'Microflora Danica','GTDB'),
#          type=if_else(grepl("LIB", user_genome), "Microflora Danica", type))%>%
#   mutate(tip.label=as.character(user_genome))%>%
#   relocate(tip.label)%>%
tax2<-tax2%>%
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))%>%
  mutate(label_3=gsub("Methylocella sp002890675", "Ca. Methyloaffinis lahnbergensis", label_3))%>%
  mutate(label_3=gsub("Methylocella sp004564215", "Methylocapsa gorgona", label_3))%>%
  mutate(label_3=gsub("Methylocella sp029855125", "Methylocapsa sp. D3K7", label_3))%>%
  filter(user_genome %in% tree_80$tip.label)


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


#################################################################################################
###################### From here, we are working with Rhodomicrobium ############################
#################################################################################################

# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Rhodomicrobium")%>%
  pull(user_genome)


tax3<-tax2%>%
  filter(user_genome %in% outgroup_tips)%>%
  left_join(drep)%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(species_rep) %>%
  mutate(genome_number2 = row_number()) %>%
  mutate(rep_label = label_3[user_genome == species_rep][1]) %>%
  mutate(rep_label = gsub("Methylocapsa |M. sp003162995|Methylocystis |sp003134075 ", "", rep_label))%>%
  mutate(rep_label = if_else(is.na(rep_label) & type=="Microflora Danica",paste0(gsub("g__", "", Genus), " MAG_", genome_number),rep_label)) %>%
  ungroup() %>%
  mutate(label_3=if_else(type=="Microflora Danica" & !user_genome==species_rep, paste0(rep_label, " sp. dup. ", genome_number2), label_3))%>%
  mutate(user_genome=as.character(user_genome))

outgroup_tips2 <- tax3 %>%
  filter(!grepl("Rhodomicrobium MAG_10|Rhodomicrobium MAG_11|Rhodomicrobium MAG_12|Rhodomicrobium MAG_13|Rhodomicrobium MAG_14|Rhodomicrobium MAG_15
                |Rhodomicrobium MAG_16|Rhodomicrobium MAG_17|Rhodomicrobium MAG_18", label_3))%>%
  filter(!grepl("Rhodomicrobium MAG_8 sp. dup. 1|Rhodomicrobium MAG_8 sp. dup. 4|Rhodomicrobium MAG_8 sp. dup. 5|Rhodomicrobium MAG_8 sp. dup. 7|Rhodomicrobium MAG_8 sp. dup. 9|Rhodomicrobium MAG_8 sp. dup. 15|Rhodomicrobium MAG_8 sp. dup. 17|Rhodomicrobium MAG_8 sp. dup. 19
                |Rhodomicrobium MAG_8 sp. dup. 20|Rhodomicrobium MAG_8 sp. dup. 22|Rhodomicrobium MAG_8 sp. dup. 27|Rhodomicrobium MAG_8 sp. dup. 26|Rhodomicrobium MAG_8 sp. dup. 28", label_3))%>%
  filter(!grepl("Rhodomicrobium MAG_38 sp. dup. 3|Rhodomicrobium MAG_38 sp. dup. 4|Rhodomicrobium MAG_38 sp. dup. 10
                |Rhodomicrobium MAG_38 sp. dup. 13|Rhodomicrobium MAG_38 sp. dup. 14|Rhodomicrobium MAG_38 sp. dup. 6|Rhodomicrobium MAG_38 sp. dup. 8
                |Rhodomicrobium MAG_38 sp. dup. 15", label_3))%>%
  pull(user_genome)

tax3<-tax3%>%
  mutate(label_3 = if_else(grepl("sp. dup", label_3),gsub("Rhodomicrobium ", "", label_3),label_3))

# remove outgroup,I use this instead
Rhodo_tree <- ape::keep.tip(tr, tip = outgroup_tips2)



Rhodo_plot<-ggtree(Rhodo_tree, linewidth=0.3) %<+% tax3 + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.2, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  geom_treescale(x=0.03, y=30, label="Tree scale", fontsize =5/.pt)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  theme(legend.position = c(0.1, 0.9)) +
  scale_x_discrete(expand = c(0,0))


Rhodo_plot

# put on bootstrap values
df1 <- Rhodo_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Rhodo_plot <- Rhodo_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  ))+
  theme(legend.position = c(0.2, 0.7),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 



ax<-Rhodo_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Rhodo_plot$data%>%filter(node==min(parent))%>%pull(y)

Rhodo_plot<-Rhodo_plot + geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
                                      aes(x = x, y = y, xend = xend, yend = yend),
                                      linewidth = 0.3,
                                      inherit.aes = FALSE) 

  

Rhodo_plot 




pl <- KEGG %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
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
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

#pl


CM <- KEGG %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(grepl("Rhodomicro", gene_label))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
#  filter(presence>0)%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.15 ) + 
  scale_fill_manual(values=c(methane_colors), na.value = "white") +
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "bottom",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

#CM

options("aplot_guides" = "keep")
met<-pl%>%aplot::insert_left(CM, width=0.08)
tree_met <-met  %>% aplot::insert_left(Rhodo_plot, width=0.4)
tree_met


#ggsave("./output/metabolism/Rhodo_tree_met_80_filt.png",tree_met_ANI, height=8.5, width=17)


ggsave("./output/metabolism_Extended_figs/Rhodomicrobium_tree_met_80_filt_25_12_11.png",tree_met,
       units = c("mm"),
       height = 145,
       width = 186,
       dpi=300)

ggsave("./output/metabolism_Extended_figs/Rhodomicrobium_tree_met_80_filt_25_12_11.svg",tree_met,
       units = c("mm"),
       height = 145,
       width = 186,
       dpi=300)





NS <- KEGG %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
 # filter(!Metabolic_step %in% c("Complex I", "Complex II", "Complex III", "Complex V"))%>%
  filter(Module %in% c("PHB","Nitrogen", "Sulphur", "CO and H2")) %>% 
 # filter(!grepl("Dsr|dsr|apr|sat|assim. nitrate", Type))%>%
  #filter(!grepl("amo|nap|ace|nrf|anfG|nasB|narB|nxrA\\/", Gene_collapsed))%>%
  # filter(grepl("nas", Gene_collapsed))%>%
  mutate(Module=gsub("Oxi. Phos.", "Complex IV", Module))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Gene_collapsed = fct_relevel(Gene_collapsed, KOs$Gene_collapsed)) %>%
  mutate(Type = factor(Type, levels = unique(KOs$Type))) %>%
  ggplot(aes(x = Gene_collapsed, y = genome))+ 
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
    "PHB"="#5B815B",
    "acetate uptake"="#77AA79",
    "EMC"="lightblue4",
    "PPP"="#AD4ACC",
    "Glycolysis"="#1237B7",
    "NADH:quinone oxidoreductase"="#b3943c",
    "Succinate dehydrogenase"="darkgreen",
    "Furmarate reductase"="#5e4fa2",
    "Ubiquinol cyt. c reductase"="#c46ca1",
    #  "menaquinol cyt. c reductase"="#679b60",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    # "Cyt. aa3-600 menaquinol"="darkred",
    "High affin. cyt c bb3-type"="darkblue",
    "F-type ATPase"="#7F9AA1"))+
  facet_nested(. ~factor(Module, levels=c("CO and H2", "Nitrogen", "Sulphur", "TCA cycle", "PHB", "EMC / acetate","PPP", "Glycolysis", "Complex IV")), scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.07, "cm"),      # space between rows
        legend.spacing.x = unit(0.07, "cm"),      # space between columns
        legend.text = element_text(size=5, margin = margin(l = 0.4)),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

NS


# ggsave("./output/metabolism_Extended_figs/Rhodomicrobium_NS_leg_25_12_11.svg",NS,
#        units = c("mm"),
#        height = 120,
#        width = 186,
#        dpi=300)

tree_met2 <-NS  %>% aplot::insert_left(Rhodo_plot, width=0.30)


ggsave("./output/metabolism_Extended_figs/Rhodomicrobium_NS_tree_met_80_filt_25_12_11.png",tree_met2,
       units = c("mm"),
       height = 120,
       width = 186,
       dpi=300)
ggsave("./output/metabolism_Extended_figs/Rhodomicrobium_NS_tree_met_80_filt_25_12_11.svg",tree_met2,
       units = c("mm"),
       height = 120,
       width = 186,
       dpi=300)







######## Now I will try to import the annotation ######
fast_ani<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Rhodomicrobium_24_07_09/ANI_Rhodo_Methylocella/fastANI.tsv", 
                col_names = c("genome", "genome_2", "ANI","aligned_contigs","total_contigs"))%>%
  mutate(genome = gsub("_genomic.fna", "", genome),
         genome = gsub(".*\\/", "", genome),
         genome = gsub(".fa","", genome),
         genome_2 = gsub("_genomic.fna", "", genome_2),
         genome_2 = gsub(".*\\/", "", genome_2),
         genome_2 = gsub(".fa","", genome_2))

order<-Rhodo_plot$data%>%
  filter(isTip==T)%>%
  arrange(desc(y)) %>%
  pull(label)

ANI_heat <- fast_ani %>%
  filter(genome %in% Rhodo_tree$tip.label,
         genome_2 %in% Rhodo_tree$tip.label)%>%
  mutate(genome=factor(genome, levels=order, ordered = TRUE))%>%
  # filter(match(genome, order) < match(genome_2, order)) %>% 
  mutate(label="Average Nucleotide\nIdentity")%>%
  ggplot(aes(genome, genome_2)) +
  #theme_minimal()+
  geom_tile(aes(fill = ANI)) +
  facet_nested(~ label, scales = "free", space = "free") +
  #geom_text(aes(label = round(ANI, 2)),size=2.2) +
  scale_fill_viridis(option="magma", begin = 1, end = 0.1)+
  labs(y=" ")+
  theme(#legend.position = 'bottom',
    axis.title.x = element_blank(),
    axis.title.y=element_text(size=3),
    strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
    strip.background = element_blank(),
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill="white"),
    legend.background = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = c(0.89, 0.88),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10))









#################################################################################################
###################### From here, we are working with Methylocella   ############################
#################################################################################################

# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax2 %>%
  filter(Genus=="g__Methylocella")%>%
  pull(user_genome)

# remove outgroup,I use this instead
Methylocella_tree <- ape::keep.tip(tr, tip = outgroup_tips)
tax2<-tax2%>%
  mutate(label_3=gsub("Methylocella", "Methylocapsa", label_3))%>%
  mutate(label_3=gsub("Ca. Methyloaffinis", "Ca. M.", label_3))%>%
  mutate(label_3=gsub("Methylocapsa sp003162995", "M. sp003162995", label_3))

order<-Methylocella_tree$tip.label

tax3<-tax2%>%
  filter(user_genome %in% order)%>%
  left_join(drep)%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(species_rep) %>%
  mutate(genome_number2 = row_number()) %>%
  mutate(rep_label = label_3[user_genome == species_rep][1]) %>%
  mutate(rep_label = gsub("Methylocapsa |M. sp003162995", "", rep_label))%>%
  mutate(rep_label = if_else(is.na(rep_label) & type=="Microflora Danica",paste0(gsub("g__", "", Genus), " MAG_", genome_number),rep_label)) %>%
  ungroup() %>%
  mutate(label_3=if_else(type=="Microflora Danica" & !user_genome==species_rep, paste0(rep_label, " sp. dup. ", genome_number2), label_3))%>%
  mutate(user_genome=as.character(user_genome))


readr::write_csv(tax3, "/home/bio.aau.dk/vj52ou/data/MFD/annotation_combined/pmoA_methylocapsa/pmoA_v2/mags_to_contigs_sp_dups.csv", col_names = FALSE)



Methylocella_plot<-ggtree(Methylocella_tree, linewidth=0.3) %<+% tax3+ 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.2, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  geom_treescale(x=0.03, y=40, label="Tree scale", fontsize =5/.pt)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  theme(legend.position = c(0.1, 0.9)) +
  scale_x_discrete(expand = c(0,0))


Methylocella_plot

# put on bootstrap values
df1 <- Methylocella_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Methylocella_plot <- Methylocella_plot %<+% df_boot + geom_nodepoint(aes(color=status),  size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  ))+
theme(legend.position = c(0.2, 0.7),
      legend.title = element_text(size=5),
      legend.text = element_text(size=5)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 



ax<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(y)

Methylocella_plot<-Methylocella_plot +   
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
                      aes(x = x, y = y, xend = xend, yend = yend),
                      linewidth = 0.3,
                      inherit.aes = FALSE)


Methylocella_plot 



######## Now I will try to import the annotation ######

# 
# fast_ani<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Rhodomicrobium_24_07_09/ANI_Rhodo_Methylocella/fastANI.tsv", 
#                 col_names = c("genome", "genome_2", "ANI","aligned_contigs","total_contigs"))%>%
#   mutate(genome = gsub("_genomic.fna", "", genome),
#          genome = gsub(".*\\/", "", genome),
#          genome = gsub(".fa","", genome),
#          genome_2 = gsub("_genomic.fna", "", genome_2),
#          genome_2 = gsub(".*\\/", "", genome_2),
#          genome_2 = gsub(".fa","", genome_2))
# 
# order<-Methylocella_plot$data%>%
#   filter(isTip==T)%>%
#   arrange(desc(y)) %>%
#   pull(label)
# 
# ANI_heat <- fast_ani %>%
#   filter(genome %in% Methylocella_tree$tip.label,
#          genome_2 %in% Methylocella_tree$tip.label)%>%
#   mutate(genome=factor(genome, levels=order, ordered = TRUE))%>%
#   # filter(match(genome, order) < match(genome_2, order)) %>% 
#   mutate(label="Average Nucleotide\nIdentity")%>%
#   ggplot(aes(genome, genome_2)) +
#   #theme_minimal()+
#   geom_tile(aes(fill = ANI)) +
#   facet_nested(~ label, scales = "free", space = "free") +
#   #geom_text(aes(label = round(ANI, 2)),size=2.2) +
#   scale_fill_viridis(option="magma", begin = 1, end = 0.1)+
#   labs(y=" ")+
#   theme(#legend.position = 'bottom',
#     axis.title.x = element_blank(),
#     axis.title.y=element_text(size=3),
#     strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
#     strip.background = element_blank(),
#     text = element_text(family = "Arial"),
#     panel.background = element_rect(fill="white"),
#     legend.background = element_blank(),
#     axis.text = element_blank(), axis.ticks = element_blank(),
#     legend.position = c(0.89, 0.88),
#     legend.text = element_text(size=10),
#     legend.title = element_text(size=10))


###3 For ISME #####



pl <- KEGG %>% 
  filter(genome %in% Methylocella_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
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
  ))+
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))


pl



pl_simple <- KEGG %>%
  filter(genome %in% Methylocella_tree$tip.label) %>%
  mutate(presence = if_else(presence == 0, NA, presence)) %>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  mutate(Gene = if_else(KO %in% c("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1", "Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylocella"), "Beijerinckiaceae_pmoA", Gene)) %>%
  mutate(Gene = if_else(KO %in% c("Root; Likely_mmoX; Rhizobiales_mmoX", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_2",
                                  "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_1", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methyloferula_Methylovirgula"), "Beijerinckiaceae mmoX", Gene)) %>%
  filter(Gene %in% c("Beijerinckiaceae_pmoA", "Beijerinckiaceae mmoX", "rbcL - RuBisCo", "rbcS - RuBisCo", "Methylocella_pxmA", "Methylocella_USCa", "Beijerinckiaceae_pxmA", "Homologous_Methylocella", "pmoA_amoA_pxmA",
                     "mxaF", "mxaI", "xoxF")) %>%
  mutate(Gene = factor(Gene, levels = c("Beijerinckiaceae_pmoA", "Beijerinckiaceae_pxmA", "Methylocella_USCa", "Methylocella_pxmA", "pmoA_amoA_pxmA", "Beijerinckiaceae mmoX","Homologous_Methylocella", "mxaF", "mxaI", "xoxF","rbcL - RuBisCo", "rbcS - RuBisCo")), ordered = TRUE) %>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test"))) %>%
  #mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  ggplot(., aes(x = Gene, y = genome)) +  
  #  geom_tile(aes(fill = Type), color="grey90", linewidth=0.09 ) +
  geom_tile(aes(fill = if_else(is.na(presence), "white", Type)), color="grey90", linewidth=0.09 ) +  
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuBis\nCO")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

options("aplot_guides" = "keep")
met<-pl%>%aplot::insert_left(pl_simple, width=0.12)
tree_met <-met  %>% aplot::insert_left(Methylocella_plot, width=0.5)
#tree_met

ggsave("./output/metabolism_Extended_figs/Methylocapsa_tree_met_80_filt.png",tree_met,
       units = c("mm"),
       height = 186,
       width = 186,
       dpi=300)











#################################################################################################
###################### From here, we are working with Methylocella  AND Methylocystis in one #####
#################################################################################################

methylocystis_drop<-tax2%>%filter(!label_3 %in% c("Methylocystis sp015709515", "Methylocapsa sp. D3K7", "Methylocapsa palsarum"))%>%
  #filter(Genus=="g__Methylocystis")%>%
  filter(type=="GTDB")%>%
  select(user_genome)

# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax2 %>%
  filter(Genus %in% c("g__Methylocella", "g__Methylocystis"))%>%
  filter(!user_genome %in% methylocystis_drop$user_genome)%>%
  pull(user_genome)

# remove outgroup,I use this instead
Methylocombo_tree <- ape::keep.tip(tr, tip = outgroup_tips)
tax2<-tax2%>%
  mutate(label_3=gsub("Methylocella", "Methylocapsa", label_3))%>%
  mutate(label_3=gsub("Ca. Methyloaffinis", "Ca. M.", label_3))%>%
  mutate(label_3=gsub("Methylocapsa sp003162995", "M. sp003162995", label_3))

order<-Methylocombo_tree$tip.label

tax3<-tax2%>%
  filter(user_genome %in% order)%>%
  left_join(drep)%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(species_rep) %>%
  mutate(genome_number2 = row_number()) %>%
  mutate(rep_label = label_3[user_genome == species_rep][1]) %>%
  mutate(rep_label = gsub("Methylocapsa |M. sp003162995|Methylocystis |sp003134075 ", "", rep_label))%>%
  mutate(rep_label = if_else(is.na(rep_label) & type=="Microflora Danica",paste0(gsub("g__", "", Genus), " MAG_", genome_number),rep_label)) %>%
  ungroup() %>%
  mutate(label_3=if_else(type=="Microflora Danica" & !user_genome==species_rep, paste0(rep_label, " sp. dup. ", genome_number2), label_3))%>%
  mutate(user_genome=as.character(user_genome))

Methylocombo_plot<-ggtree(Methylocombo_tree, linewidth=0.3) %<+% tax3+ 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.3, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  geom_treescale(x=0.03, y=40, label="Tree scale", fontsize =5/.pt)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  theme(legend.position = c(0.1, 0.9)) +
  scale_x_discrete(expand = c(0,0))


Methylocombo_plot

# put on bootstrap values
df1 <- Methylocombo_plot %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Methylocombo_plot <- Methylocombo_plot %<+% df_boot + geom_nodepoint(aes(color=status),  size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  ))+
  theme(legend.position = c(0.2, 0.7),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 



ax<-Methylocombo_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocombo_plot$data%>%filter(node==min(parent))%>%pull(y)

Methylocombo_plot<-Methylocombo_plot +   
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3,
               inherit.aes = FALSE)


Methylocombo_plot 





pl <- KEGG %>% 
  filter(genome %in% Methylocombo_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
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
  ))+
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))


pl


genes_keep<-KEGG%>%filter(genome %in% Methylocombo_tree$tip.label)%>%
  filter(Metabolic_step %in% c("Custom\nHMM"))%>%
  group_by(Gene)%>%
  filter(presence==1)%>%
  select(Gene)%>%
  distinct()%>%
  filter(!Gene%in%c("Homologous_MO", "Rhodomicrobium_umbrella", "Rhodomicrobium_likely_mmoX", "Homologous_Methylocella"))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))





pl_simple <- KEGG %>%
  filter(genome %in% Methylocombo_tree$tip.label) %>%
  mutate(presence = if_else(presence == 0, NA, presence)) %>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(Gene %in% genes_keep$Gene)%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pxmA", "Beijerinckiaceae_pxmA", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA1_1", "Methylosinus_Methylocystis_pmoA1", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA2_1", "Methylosinus_Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = if_else(KO %in% c("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1", "Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylocella"), "Beijerinckiaceae_pmoA", Gene)) %>%
  mutate(Gene = if_else(KO %in% c("Root; Likely_mmoX; Rhizobiales_mmoX", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_2",
                                  "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_1", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methyloferula_Methylovirgula"), "Beijerinckiaceae mmoX", Gene)) %>%
  mutate(Gene=gsub("Beijerinckiaceae mmoX", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=gsub("Beijerinckiaceae_pmoA", "Methylosinus_Methylocystis_pmoA1", Gene))%>%
  #filter(Gene %in% c("Beijerinckiaceae_pmoA", "Beijerinckiaceae mmoX", "rbcL - RuBisCo", "rbcS - RuBisCo", "Methylocella_pxmA", "Methylocella_USCa", "Beijerinckiaceae_pxmA", "Homologous_Methylocella", "pmoA_amoA_pxmA",                  "mxaF", "mxaI", "xoxF")) %>%
 # mutate(Gene = factor(Gene, levels = c("Beijerinckiaceae_pmoA", "Beijerinckiaceae_pxmA", "Methylocella_USCa", "Methylocella_pxmA", "pmoA_amoA_pxmA", "Beijerinckiaceae mmoX","Homologous_Methylocella", "mxaF", "mxaI", "xoxF","rbcL - RuBisCo", "rbcS - RuBisCo")), ordered = TRUE) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test"))) %>%
  #mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  ggplot(., aes(x = Gene, y = genome)) +  
  #  geom_tile(aes(fill = Type), color="grey90", linewidth=0.09 ) +
  geom_tile(aes(fill = if_else(is.na(presence), "white", Type)), color="grey90", linewidth=0.09 ) +  
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuBis\nCO")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

options("aplot_guides" = "keep")
met<-pl%>%aplot::insert_left(pl_simple, width=0.22)
tree_met <-met  %>% aplot::insert_left(Methylocombo_plot, width=0.5)
#tree_met

ggsave("./output/metabolism_Extended_figs/Methylocombo_tree_met_80_filt.png",tree_met,
       units = c("mm"),
       height = 250,
       width = 186,
     # width = 200,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/Methylocombo_tree_met_80_filt.svg",tree_met,
       units = c("mm"),
       height = 250,
       width = 186,
       # width = 200,
       dpi=300)




#################################################################################
######################### Adding on Nitrogen, Sulphur, and other metabolic stuff 
#################################################################################

cytochromes <-  KEGG %>%
  filter(genome %in% Methylocombo_tree$tip.label) %>%
  filter(Module == 'Oxi. Phos.') %>% 
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence), "white", Type)), color="grey90", linewidth=0.1 )+
  #  geom_tile(aes(fill = if_else(is.na(presence),  "white", Pathway)), color="grey90", linewidth=0.1 ) + 
  scale_fill_manual(na.value="transparent", values=c(
    "NADH:quinone oxidoreductase"="#b3943c",
    "Succinate dehydrogenase"="darkgreen",
    "Furmarate reductase"="#5e4fa2",
    "Ubiquinol cyt. c reductase"="#c46ca1",
    #  "menaquinol cyt. c reductase"="#679b60",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    # "Cyt. aa3-600 menaquinol"="darkred",
    "High affin. cyt c bb3-type"="darkblue",
    "F-type ATPase"="#7F9AA1"))+
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )

cytochromes


NS <- KEGG %>% 
  filter(genome %in% Methylocombo_tree$tip.label) %>%
  filter(!Metabolic_step %in% c("Complex I", "Complex II", "Complex III", "Complex V"))%>%
  filter(Module %in% c("PHB","Nitrogen", "Sulphur", "CO and H2", "Oxi. Phos.")) %>% 
  filter(!grepl("Dsr|dsr|apr|sat|assim. nitrate", Type))%>%
  filter(!grepl("amo|nap|ace|nrf|anfG|nasB|narB|nxrA\\/", Gene_collapsed))%>%
 # filter(grepl("nas", Gene_collapsed))%>%
  mutate(Module=gsub("Oxi. Phos.", "Complex IV", Module))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Gene_collapsed = fct_relevel(Gene_collapsed, KOs$Gene_collapsed)) %>%
  mutate(Type = factor(Type, levels = unique(KOs$Type))) %>%
  ggplot(aes(x = Gene_collapsed, y = genome))+ 
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
    "PHB"="#5B815B",
    "acetate uptake"="#77AA79",
    "EMC"="lightblue4",
    "PPP"="#AD4ACC",
    "Glycolysis"="#1237B7",
    "NADH:quinone oxidoreductase"="#b3943c",
    "Succinate dehydrogenase"="darkgreen",
    "Furmarate reductase"="#5e4fa2",
    "Ubiquinol cyt. c reductase"="#c46ca1",
    #  "menaquinol cyt. c reductase"="#679b60",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    # "Cyt. aa3-600 menaquinol"="darkred",
    "High affin. cyt c bb3-type"="darkblue",
    "F-type ATPase"="#7F9AA1"))+
  facet_nested(. ~factor(Module, levels=c("CO and H2", "Nitrogen", "Sulphur", "TCA cycle", "PHB", "EMC / acetate","PPP", "Glycolysis", "Complex IV")), scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.07, "cm"),      # space between rows
        legend.spacing.x = unit(0.07, "cm"),      # space between columns
        legend.text = element_text(size=5, margin = margin(l = 0.4)),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

NS




tree_met2 <-NS  %>% aplot::insert_left(Methylocombo_plot, width=0.45)


ggsave("./output/metabolism_Extended_figs/Methylocombo_tree_met_80_filt_NS_cyt.png",tree_met2,
       units = c("mm"),
       height = 285,
       width = 190,
       # width = 200,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/Methylocombo_tree_met_80_filt_NS_cyt.svg",tree_met2,
       units = c("mm"),
       height = 285,
       width = 190,
       # width = 200,
       dpi=300)


#################################################################################
########################### Methylocystis #######################################
#################################################################################


# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocystis")%>%
  pull(user_genome)

# remove outgroup,I use this instead
Methylocystis_tree <- ape::keep.tip(tr, tip = outgroup_tips)
Methylocystis_plot<-ggtree(Methylocystis_tree) %<+% tax2 + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.07, size = 2,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  geom_treescale(x=0.03, y=30, label="Tree scale", fontsize =3)+
  # scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#FCEA8C", "LR-MAG"="#AFEEEE", "SR-MAG"="#B8BEFC"))+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  # geom_polygon(aes(fill = type, x = 0, y = 0)) +
  theme(legend.position = c(0.1, 0.9)) +
  scale_x_discrete(expand = c(0,0))

Methylocystis_plot

# put on bootstrap values
df1 <- Methylocystis_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


Methylocystis_plot <- Methylocystis_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=1.4) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.9)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=10),
    label.theme = element_text(size=9),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 2)
  )) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=10),
    label.theme = element_text(size=9),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 3)
  ))


ax<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(y)

Methylocystis_plot<-Methylocystis_plot + geom_segment(aes(x = ax, y = ay, xend = -0.05, yend = ay),
                                                      arrow = arrow(length = unit(0.2, "cm"), type="closed"), lwd=0.2)+
  annotate("text", x = -0.025, y = ay+1.5, label = "Outgroup", size=3)+
  theme(legend.position = c(0.2, 0.9)) 

Methylocystis_plot 

fast_ani<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/ANI_Methylocystis/fastANI.tsv", 
                col_names = c("genome", "genome_2", "ANI","aligned_contigs","total_contigs"))%>%
  mutate(genome = gsub("_genomic.fna", "", genome),
         genome = gsub(".*\\/", "", genome),
         genome = gsub(".fa","", genome),
         genome_2 = gsub("_genomic.fna", "", genome_2),
         genome_2 = gsub(".*\\/", "", genome_2),
         genome_2 = gsub(".fa","", genome_2))

order<-Methylocystis_plot$data%>%
  filter(isTip==T)%>%
  arrange(desc(y)) %>%
  pull(label)

ANI_heat <- fast_ani %>%
  filter(genome %in% Methylocystis_tree$tip.label,
         genome_2 %in% Methylocystis_tree$tip.label)%>%
  mutate(genome=factor(genome, levels=order, ordered = TRUE))%>%
  # filter(match(genome, order) < match(genome_2, order)) %>% 
  mutate(label="Average Nucleotide\nIdentity")%>%
  ggplot(aes(genome, genome_2)) +
  #theme_minimal()+
  geom_tile(aes(fill = ANI)) +
  facet_nested(~ label, scales = "free", space = "free") +
  #geom_text(aes(label = round(ANI, 2)),size=2.2) +
  scale_fill_viridis(option="magma", begin = 1, end = 0.1)+
  labs(y=" ")+
  theme(#legend.position = 'bottom',
    axis.title.x = element_blank(),
    axis.title.y=element_text(size=3),
    strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
    strip.background = element_blank(),
    strip.clip = "off",
    panel.background = element_rect(fill="white"),
    legend.background = element_blank(),
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = c(0.83, 0.78),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10))
#ANI_heat




#tree_ANI <- ANI_heat %>% aplot::insert_left(Methylocystis_plot, width=0.7)
#tree_ANI


################################ Getting the metabolism to work ##################################################

#lc <- unlist(methane_colors)

pl <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  #  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
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
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )

#pl

genes_keep<-KEGG%>%filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()%>%
  filter(!gene_label%in%c("Homologous_MO (1/3)", "Rhodomicrobium_umbrella (1/3)"))%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA1_", gene_label), "Methylocystis_pmoA1 (1/3)", gene_label))%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA2_", gene_label), "Methylocystis_pmoA2 (1/3)", gene_label))%>%
  mutate(gene_label=gsub("Methylosinus_umbrella", "Methylosinus", gene_label))

KOs<-KOs%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA1_", gene_label), "Methylocystis_pmoA1 (1/3)", gene_label))%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA2_", gene_label), "Methylocystis_pmoA2 (1/3)", gene_label))%>%
  mutate(gene_label=gsub("Methylosinus_umbrella", "Methylosinus", gene_label))
  

CM <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA1_", gene_label), "Methylocystis_pmoA1 (1/3)", gene_label))%>%
  mutate(gene_label=if_else(grepl("Methylocystis_pmoA2_", gene_label), "Methylocystis_pmoA2 (1/3)", gene_label))%>%
  mutate(gene_label=gsub("Methylosinus_umbrella", "Methylosinus", gene_label))%>%
  filter(gene_label %in% genes_keep$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_blank(),
        axis.ticks.x =element_line(linewidth = 0.1),
        # axis.text.y=element_text(size=7),
        axis.ticks.y = element_line(linewidth = 0.1), 
        panel.spacing = unit(0.03, "cm", data = NULL),
        strip.text.x = element_text(size = 8, face="bold"),  # Size for x-axis facet labels
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),legend.position = "none",
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

CM

options("aplot_guides" = "keep")
met<-pl%>%aplot::insert_left(CM, width=0.2)
tree_met <-met  %>% aplot::insert_left(Methylocystis_plot, width=0.72)
#tree_met
ggsave("./output/metabolism/Methylocystis_tree_met_80_filt_no_ANI_25_06_17.png",tree_met, height=10, width=16)


tree_met_ANI <-tree_met  %>% aplot::insert_right(ANI_heat, width=0.6)
#tree_met_ANI

ggsave("./output/metabolism/Methylocystis_tree_met_80_filt_25_06_17.png",tree_met_ANI, height=9, width=18 )




#################################################################################################
###################### From here, we are working with Gammaproteobacteria #######################
#################################################################################################


#### I cannot handle Gamma as one big, I need some subdivisions here! 
Gamma_genera<-tax%>%
  #filter(type=="Microflora Danica")%>%
  filter(Class=="c__Gammaproteobacteria")%>%
  filter(!Genus%in% c("g__Methyloprofundus", "g__Methylumidiphilus"))%>%
  select(Genus)%>%
  pull(unique(Genus))


outgroup_tips <- tax %>% 
  filter(Class=="c__Gammaproteobacteria")%>%
  filter(Genus %in% Gamma_genera)%>%
  filter(user_genome %in% unique(KEGG$genome))%>%
  filter(!type=="GTDB"  | Species %in% c("s__Methylobacter_C polyspora_A", "s__Methylobacter marinus", "s__Methylobacter_C polyspora_A", "s__JAGXGJ01 sp028714175", "s__JAGXGJ01 sp024636875", "s__JAGXGJ01 sp028715725",
                                         "s__Methyloglobulus morosus", "s__Crenothrix polyspora_B", "s__Crenothrix polyspora", "s__Methylovulum miyakonense", "s__Methylovulum oryzae", "s__Methylobacter_B favarea",
                                         "s__Methylobacter_A oryzae_A", "s__Methylobacter_A tundripaludum_A", "s__Methylobacter_A titanis", "s__Methylomonas paludis", "s__Methylomonas methanica_A", "s__Methylomonas fluvii",
                                         "s__Methylomonas koyamae", "s__Methylomonas lenta", "s__Methylomarinum vadi", "s__CAIQWF01 sp903875455", "s__CAIQWF01 sp020048335", 
                                         "s__USCg-Taylor sp030859565", "s__USCg-Taylor sp002007425", "s__JACCXJ01 sp013697045", "s__JACCXJ01 sp030860485"))%>%
  pull(user_genome)%>%as.character()

tax3<-tax2%>%
  filter(user_genome %in% outgroup_tips)%>%
  left_join(drep)%>%
  mutate(user_genome = factor(user_genome, levels = order, ordered = TRUE)) %>%
  arrange(user_genome) %>%
  group_by(species_rep) %>%
  mutate(genome_number2 = row_number()) %>%
  mutate(rep_label = label_3[user_genome == species_rep][1]) %>%
  mutate(rep_label = if_else(is.na(rep_label) & type=="Microflora Danica",paste0(gsub("g__", "", Genus), " MAG_", genome_number),rep_label)) %>%
  ungroup() %>%
  mutate(label_3=if_else(type=="Microflora Danica" & !user_genome==species_rep, paste0(rep_label, " sp. dup. ", genome_number2), label_3))%>%
  mutate(user_genome=as.character(user_genome))

outgroup_tips2 <- tax3 %>%
  filter(!grepl("Crenothrix", label_3))%>%
  filter(!(type=="GTDB"&grepl("Methylobacter|Methylomonas|Methyloglobulus|Methylovulum|Methylomarinum|JAGXGJ01|CAIQ", label_3)))%>%
  pull(user_genome)

tax3<-tax3%>%
  mutate(label_3 = if_else(grepl("sp. dup", label_3),gsub("Rhodomicrobium ", "", label_3),label_3))%>%
  mutate(label_3=gsub("USCg-Taylor sp002007425 MAG_2 sp. dup. 1", "MAG_2 sp. dup. 1", label_3))%>%
  mutate(label_3=gsub("Methylobacter_A MAG_15 sp. dup. ", "MAG_15 sp. dup. ", label_3))

# remove outgroup,I use this instead
Gamma_tree <- ape::keep.tip(tr, tip = outgroup_tips2)



Gamma_plot<-ggtree(Gamma_tree, linewidth=0.3) %<+% tax3 + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.8, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  geom_treescale(x=0.03, y=30, label="Tree scale", fontsize =5/.pt)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  theme(legend.position = c(0.1, 0.9),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
  scale_x_discrete(expand = c(0,0))


Gamma_plot

# put on bootstrap values
df1 <- Gamma_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Gamma_plot <- Gamma_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  ))+
  theme(legend.position = c(0.2, 0.7),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 



ax<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(y)

Gamma_plot<-Gamma_plot + geom_segment(data = data.frame(x = ax, y = ay, xend = -0.04, yend = ay),
                                      aes(x = x, y = y, xend = xend, yend = yend),
                                      linewidth = 0.3,
                                      inherit.aes = FALSE) 



Gamma_plot 





# Identify the tip labels corresponding to the outgroup sequences

pl <- KEGG %>% 
  filter(genome %in% Gamma_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.07, "cm"),      # space between rows
        legend.spacing.x = unit(0.07, "cm"),      # space between columns
        legend.text = element_text(size=5, margin = margin(l = 0.4)),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

pl


genes_keep<-KEGG %>% 
  filter(genome %in% Gamma_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%distinct()


CM <- KEGG %>% 
  filter(genome %in% Gamma_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(gene_label %in% genes_keep$gene_label)%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methylovulum_mmoX (1/3)", paste0("Methylomonadaceae_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methylomonas_mmoX (1/3)", paste0("Methylomonadaceae_mmoX (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methylomonas (1/3)", paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methylobacter_clade_1 (1/3)", paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methylobacter_A (1/3)", paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  mutate(gene_label=if_else(presence==1&gene_label=="Methyloglobulus (1/3)", paste0("Methylomonadaceae (1/3)"), gene_label))%>%
  filter(!gene_label %in% c("Methylomonas_mmoX (1/3)", "Methylovulum_mmoX (1/3)", "Methylomonas (1/3)", 
                            "Methylobacter_clade_1 (1/3)", "Methylobacter_A (1/3)", "Methyloglobulus (1/3)"))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.15 ) + 
  scale_fill_manual(values=c(methane_colors), na.value = "white") +
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        legend.position = "none",
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.07, "cm"),      # space between rows
        legend.spacing.x = unit(0.07, "cm"),      # space between columns
        legend.text = element_text(size=5, margin = margin(l = 0.4)),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

#CM

options("aplot_guides" = "keep")
met<-pl%>%aplot::insert_left(CM, width=0.12)
tree_met <-met  %>% aplot::insert_left(Gamma_plot, width=0.35)
#tree_met_no_HMM <-pl  %>% aplot::insert_left(Gamma_plot, width=1)
#tree_met




ggsave("./output/metabolism_Extended_figs/Gamma_tree_met_80_25_12_11.png",tree_met,
       units = c("mm"),
       height = 100,
       width = 190,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/Gamma_tree_met_80_25_12_11.svg",tree_met,
       units = c("mm"),
       height = 100,
       width = 190,
       dpi=300)



#ggsave("./output/metabolism/Gamma_tree_met_80_filt_25_08_25.png",tree_met, height=6, width=15, dpi=400)




NS <- KEGG %>% 
  filter(genome %in% Gamma_tree$tip.label)%>%
#  filter(!Metabolic_step %in% c("Complex I", "Complex II", "Complex III", "Complex V"))%>%
  filter(Module %in% c("PHB","Nitrogen", "Sulphur", "CO and H2")) %>% 
 # filter(!grepl("Dsr|dsr|apr|sat|assim. nitrate", Type))%>%
  filter(!grepl("anfG", Gene_collapsed))%>%
  # filter(grepl("nas", Gene_collapsed))%>%
  mutate(Module=gsub("Oxi. Phos.", "Complex IV", Module))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Gene_collapsed = fct_relevel(Gene_collapsed, KOs$Gene_collapsed)) %>%
  mutate(Type = factor(Type, levels = unique(KOs$Type))) %>%
  ggplot(aes(x = Gene_collapsed, y = genome))+ 
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
    "PHB"="#5B815B",
    "acetate uptake"="#77AA79",
    "EMC"="lightblue4",
    "PPP"="#AD4ACC",
    "Glycolysis"="#1237B7",
    "NADH:quinone oxidoreductase"="#b3943c",
    "Succinate dehydrogenase"="darkgreen",
    "Furmarate reductase"="#5e4fa2",
    "Ubiquinol cyt. c reductase"="#c46ca1",
    #  "menaquinol cyt. c reductase"="#679b60",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    # "Cyt. aa3-600 menaquinol"="darkred",
    "High affin. cyt c bb3-type"="darkblue",
    "F-type ATPase"="#7F9AA1"))+
  facet_nested(. ~factor(Module, levels=c("CO and H2", "Nitrogen", "Sulphur", "TCA cycle", "PHB", "EMC / acetate","PPP", "Glycolysis", "Complex IV")), scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.y = unit(0.07, "cm"),      # space between rows
        legend.spacing.x = unit(0.07, "cm"),      # space between columns
        legend.text = element_text(size=5, margin = margin(l = 0.4)),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        ) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

NS




tree_met2 <-NS  %>% aplot::insert_left(Gamma_plot, width=0.4)

ggsave("./output/metabolism_Extended_figs/Gamma_tree_met_NS_80_25_12_11.png",tree_met2,
       units = c("mm"),
       height = 105,
       width = 190,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/Gamma_tree_met_NS_80_25_12_11.svg",tree_met2,
       units = c("mm"),
       height = 105,
       width = 190,
       dpi=300)




###################### Getting all of the USCg ########################


KEGG<-readRDS("./output/metabolism/KEGG_25_12_02.rds")%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))

tax_oxi_rate<-vroom("/projects/microflora_danica/mmlong2/extra/mmlong2_bins.tsv")%>%
  select(bin, gtdb_tax)%>%
  rename(classification=gtdb_tax)%>%
  rename(genome=bin)%>%
  separate(classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ';', remove = FALSE)

pmoA<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/methanotrophs/graftm_output/graft_oxi_rates/25_10_07_pmoA/combined_count_table_e10.txt", delim="\t")
pmoA<-pmoA%>%
  pivot_longer(cols=2:length(pmoA), names_to="genome", values_to="presence")%>%
  filter(presence>0)%>%
  distinct()%>%
  rename(KO = ConsensusLineage)




KEGG_fixed <- KEGG %>%
  left_join(tax_oxi_rate, by = "genome", suffix = c("", "_new")) %>%
  mutate(
    classification = coalesce(classification, classification_new),
    Domain        = coalesce(Domain, Domain_new),
    Phylum        = coalesce(Phylum, Phylum_new),
    Class         = coalesce(Class, Class_new),
    Order         = coalesce(Order, Order_new),
    Family        = coalesce(Family, Family_new),
    Genus         = coalesce(Genus, Genus_new),
    Species       = coalesce(Species, Species_new)
  ) %>%
  select(-ends_with("_new"))%>%
  left_join(pmoA, by = c("genome", "KO"), suffix = c("", "_pmoA")) %>%
  mutate(
    presence = if_else(
      presence == 0 & presence_pmoA > 0,
      1,
      presence
    )
  ) %>%
  select(-presence_pmoA)%>%
  mutate(label=if_else(genome %in% tax_oxi_rate$genome, paste0("Dokkedal ", Genus,";", Species, " | ", genome), label))




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

KEGG_sort<-KEGG_fixed%>%
  arrange(Genus)

KEGG<-KEGG_fixed%>%
  mutate (gene_label = paste0(Gene_collapsed , ' ',  Pathway_step))

pl <- KEGG %>% 
  filter(Family %in% c("f__JACCXJ01"))%>%
  filter(grepl("MFD|LIB", genome))%>%
  filter(Module == 'C1') %>% 
  filter(!Metabolic_step=="Custom\nHMM")%>%
  filter(!grepl("mmo", gene_label))%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Type=fct_relevel(Type, unique(KOs$Type)))%>%
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
  facet_nested(.~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
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
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"))

pl



genes_keep<-KEGG %>% 
  filter(Family %in% c("f__JACCXJ01"))%>%
  filter(grepl("MFD|LIB", genome))%>%
  filter(label %in% pl$data$label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%distinct()


CM <- KEGG %>% 
  filter(Family %in% c("f__JACCXJ01"))%>%
  filter(grepl("MFD|LIB", genome))%>%
  filter(label %in% pl$data$label)%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(gene_label %in% genes_keep$gene_label)%>%
  # filter(grepl("Rhodomicro", gene_label))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  #  filter(presence>0)%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(aes(x = gene_label, y = label))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.15 ) + 
  scale_fill_manual(values=c(methane_colors), na.value = "white") +
  facet_nested(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
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
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"))

CM


options("aplot_guides" = "keep")



uscg<-pl%>%aplot::insert_left(CM, width=0.09)
uscg

ggsave("./output/metabolism_Extended_figs/USCg_all_met_25_12_11.png",uscg,
       units = c("mm"),
       height = 60,
       width = 186,
       dpi=300)



ggsave("./output/metabolism_Extended_figs/USCg_all_met_25_12_11.svg",uscg,
       units = c("mm"),
       height = 60,
       width = 186,
       dpi=300)
























# Look here when I want to label the clades of the tree https://4va.github.io/biodatasci/r-ggtree.html 
fast_ani<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Rhodomicrobium_24_07_09/ANI_Rhodo_Methylocystis/fastANI.tsv", 
                col_names = c("genome", "genome_2", "ANI","aligned_contigs","total_contigs"))%>%
  mutate(genome = gsub("_genomic.fna", "", genome),
         genome = gsub(".*\\/", "", genome),
         genome = gsub(".fa","", genome),
         genome_2 = gsub("_genomic.fna", "", genome_2),
         genome_2 = gsub(".*\\/", "", genome_2),
         genome_2 = gsub(".fa","", genome_2))

order<-Rhodo_plot$data%>%
  filter(isTip==T)%>%
  arrange(desc(y)) %>%
  pull(label)

library(viridis)
ANI_heat <- fast_ani %>%
  filter(genome %in% Rhodo_tree$tip.label,
         genome_2 %in% Rhodo_tree$tip.label)%>%
  mutate(genome=factor(genome, levels=order, ordered = TRUE))%>%
  # filter(match(genome, order) < match(genome_2, order)) %>% 
  mutate(label="Average\nNucleotide\nIdentity")%>%
  ggplot(aes(genome, genome_2)) +
  theme_nothing()+
  geom_tile(aes(fill = ANI)) +
  facet_grid(~ label, scales = "free", space = "free") +
  #geom_text(aes(label = round(ANI, 2)),size=2.2) +
  scale_fill_viridis(option="magma", begin = 0.95, end = 0.1)+
  theme(#legend.position = 'bottom',
    axis.title = element_blank(),
    strip.text.x = element_text(size = 8), 
    strip.background = element_rect(fill = "grey90", color = "white"),
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = c(0.8, 0.83),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10))
ANI_heat


tree_met_ANI <-tree_met  %>% aplot::insert_right(ANI_heat, width=0.7)
tree_met_ANI

ggsave("./output/metabolism/Rhodo_tree_met_80_filt.png",tree_met_ANI, height=10, width=20 )



# 
# ANI_heat <- fast_ani %>%
#   filter(genome %in% Rhodo_tree$tip.label,
#          genome_2 %in% Rhodo_tree$tip.label)%>%
#   mutate(genome = factor(genome, levels = order, ordered = TRUE),
#          ANI_category = cut(ANI, breaks = c(75, 80, 85, 90, 95, Inf), 
#                             labels = c("<80", "80-85", "85-90", "90-95", ">95"), 
#                             include.lowest = TRUE)) %>%
#   ggplot(aes(genome, genome_2)) +
#   theme_nothing()+ geom_tile(aes(fill = ANI_category)) +
#   # geom_text(aes(label = round(ANI, 2)), size = 2.2) +
#   scale_fill_manual(values = c("<80" = "white",
#                                "80-85" ="#C2DACA",
#                                "85-90" = "#97CFC0",
#                                "90-95" = "#438B7A",
#                                ">95" = "#22473F"), 
#                     na.value = "#22473F") +
#   theme(legend.position = 'bottom',
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank())
# 
# ANI_heat


################## To change the colors #######################
# methane_colors <- structure(list(`pMMO` = "#b3943c",
#                                  `sMMO` = "darkgreen",
#                                  `calcium (mxa)` = "#5e4fa2",
#                                  `lanthanide (xoxF)` = "#c46ca1",
#                                  putative = "grey75",
#                                  H4MPT = "#d97512",
#                                  GSH = "turquoise4",
#                                  `H4F` = "darkred",
#                                  `Serine cycle` = "darkred",
#                                  `RuMP cycle` = "darkblue" ,
#                                  `CBB cycle` = "#679b60",
#                                  `-` = "#679b60"), class = "character",
#                             row.names = unique_types)
# saveRDS(methane_colors, "palette_methane_colors.rds")
# 

custom_colors <-c("#679b60","#5e4fa2" , "#679b60", "turquoise4", "darkred","#d97512","#c46ca1" , "#b3943c", "grey75","darkblue","darkred","darkgreen" )




####################  Methylocella, not Rhodomicrobium ##################

outgroup_Rhodo <- tax %>%
  filter(Genus=="g__Rhodomicrobium")%>%
  pull(user_genome)


#### Below is to collapse the clades. I will not really use this
out_rho <- getMRCA(tr,outgroup_Rhodo)
Rh_collapsed <- pr %>%  scaleClade(out_rho, .1)%>% collapse(out_rho, clade_name = 'Rhodomicrobium', mode = 'min',
                                                            fill="grey80",
                                                            color= 'black') 

f <- Rh_collapsed %<+% tax + 
  geom_tiplab(aes(label = Species, fill = type),
              geom = "label", label.size = 0.2, hjust = 0, align=TRUE)



# remove outgroup,I use this instead
Rhodo_remove <- ape::drop.tip(tr, tip = outgroup_Rhodo)
f2_Methylocella<-ggtree(Rhodo_remove) %<+% tax + 
  geom_tiplab(aes(label = Species, fill = type),
              geom = "label",  size = 2.5, linesize = 0.1, label.size = 0.0, text.size=6, offset = -0.010, align=TRUE, hjust=0) +
  geom_treescale(x=0.01, y=15)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "darkgreen", "LR-MAG"="#AFEEEE", "SR-MAG"="#5F9EA0"))+
  guides(fill = guide_legend(override.aes = list(label = "\u25A0", 
                                                 size = 0)))

f2_Methylocella



pl <- KEGG_tree %>% 
  filter(Genus=="g__Methylocella")%>%
  filter(Module == 'C1') %>% 
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(., aes(x = gene_label, y = genome))+ 
  geom_tile(data = . %>% filter(presence > 0),
            aes(fill = Type), color="grey90", linewidth=0.1) +  # Set a fixed width for the tiles (e.g., 0.5)
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(values=methane_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_text(size=8.3),
        axis.title = element_blank(),legend.position = "none" # This removes the legend)
  )


##### Filtering the tree to match the heatmap #####

x2<-KEGG_tree %>% 
  filter(Genus=="g__Methylocella")%>%
  filter(Module == 'C1') %>% 
  filter(presence>0)%>%
  select(genome)%>%
  distinct()%>%
  pull(genome)

tree_remove_2 <- ape::keep.tip(Rhodo_remove, tip = x2)


tree2<-ggtree(tree_remove_2)  %<+% tax + 
  geom_tiplab(aes(label = Species, fill = type),
              geom = "label",  size = 2.5, linesize = 0.1, label.size = 0.0, text.size=6, offset = -0.010, align=TRUE, hjust=0) +
  geom_treescale(x=0.01, y=15)+
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "purple1", "LR-MAG"="#AFEEEE", "SR-MAG"="#5F9EA0"))+
  geom_treescale(x=0.01, y=15)

tree2$data



comb <-pl  %>% aplot::insert_left(tree2)
comb
intersect(tree_remove_2$tip.label, pl$data$genome)














########################################################################################################
################################# Getting the legend ###################################################


pl <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, unique(KOs$Metabolic_step))) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, unique(KOs$gene_label))) %>%
  ggplot(., aes(x = gene_label, y = genome))+ 
  geom_tile(
    #data = .  %>% filter(presence > 0),
    aes(fill = Type), color="grey90", linewidth=0.1) +  # Set a fixed width for the tiles (e.g., 0.5)
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(values=methane_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), 
        axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        axis.ticks.y = element_line(linewidth = 0.1), 
        legend.key.size = unit(0.45, "cm"),
        strip.text.x = element_text(size = 8),
        axis.title = element_blank(),legend.position = "bottom" # This removes the legend)
  )+
  guides(fill = guide_legend(nrow = 2, title.theme = element_blank()))



# Getting the legend
legend <- ggpubr::get_legend(pl)
leg<-ggpubr::as_ggplot(legend)
ggsave("./output/metabolism/legend.svg", leg)


