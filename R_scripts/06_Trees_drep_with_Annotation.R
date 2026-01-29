
.libPaths()

#Works as /home/bio.aau.dk/vj52ou/R/rstudio-server/4.4.0

# 
# installed <- installed.packages()[, "Package"]
# print(installed)
# remove.packages(installed)
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("remotes")
# remotes::install_github("YuLab-SMU/ggtree")



#.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.4.0", .libPaths()))
#
#.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.5.0"))
#.libPaths(c("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4"))

#.libPaths("/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")


#.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))

#library(tidyr)

library(ggtree)
library(ggplot2)


#library(cowplot)
#library(vroom)
#library(tidyr)
#library(readxl)

#library(ggtreeExtra)
#library(ape)
#library(ggh4x)
library(treeio)
#library(viridis)
library(dplyr)
library(stringr)
library(ggh4x)
#library(patchwork)
#library(tidyverse, lib.loc ="/home/bio.aau.dk/kl42gg/R/x86_64-pc-linux-gnu-library/4.4")

setwd("~/scripts/MFD/methanotrophs/R_scripts/")

# Tree file

tree_drep <- treeio::read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/with_gtdb_25_03_04/MSA.faa.treefile")
#tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/sp_reps_sylph_mfd_25_03_06/MSA.faa.treefile")



###################################################################################################################
################# OBS! sylph and genome number (label_3 )  ###########
###################################################################################################################




tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

tax<-readRDS("MFD_renamed_tax_25_03_04.rds")%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'Microflora Danica','GTDB'),
         type=if_else(grepl("LIB", user_genome), "Microflora Danica", type))%>%
  mutate(tip.label=as.character(user_genome))%>%
  relocate(tip.label)%>%
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))%>%
  mutate(label_3=gsub("Methylocella sp002890675", "Ca. Methyloaffinis lahnbergensis", label_3))%>%
  mutate(label_3=gsub("Methylocella sp004564215", "Methylocapsa gorgona", label_3))%>%
  mutate(label_3=gsub("Methylocella sp029855125", "Methylocapsa sp. D3K7", label_3))
  
duplicates <- tax %>%
  group_by(Species) %>%
  filter(n() > 1) %>%
  ungroup()%>%
  filter(type=="GTDB")%>%
  pull(tip.label)

tax<-tax%>%filter(!tip.label %in% duplicates)

#tax_old<-readRDS("MFD_renamed_tax.rds")


tr <- phytools::midpoint.root(tree_drep)
pr<-ape::drop.tip(tr, tip = duplicates)
# plot rooted tree
pr_tree <- ggtree::ggtree(pr) 
pr_tree


sylph<-readRDS("output/genome_quantification/sylph_gtdb_25_03_06.rds")%>%
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))


sylph<-sylph%>%
  filter(!is.na(fieldsample_barcode))%>%
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
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))%>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))


levels_hab2<-readRDS("./output/genome_quantification/levels_hab2_sylph_DMS.rds")
palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

hab2_sort<-sylph%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

KEGG<-readRDS("./output/metabolism/KEGG_25_10_10.rds")%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))

methane_colors<-readRDS("palette_methane_colors.rds")%>%
  unlist()
names(methane_colors)[names(methane_colors) == "-"] <- "Formate oxidation"

KOs<-readxl::read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)

KEGG_sort<-KEGG%>%
  arrange(Genus)


hab2_sort_order<-readRDS("./output/hab2_sort_order_ISME.rds")



library(tidyverse)


#################################################################################################
###################### From here, we are working with Rhodomicrobium ############################
#################################################################################################


# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Rhodomicrobium")%>%
  pull(tip.label)
#library(ggtreeExtra)

Rhodo_tree <- ape::keep.tip(pr, tip = outgroup_tips)
Rhodo_plot<-ggtree(Rhodo_tree, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"),guide = "none")+
  theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0))


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
  theme(legend.position = c(0.15, 0.87),
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

Rhodo_plot <- Rhodo_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3,
               inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)


Rhodo_plot 







heat <- sylph %>%
  filter(user_genome %in% Rhodo_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


heat

pl_simple <- KEGG %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(grepl("rbcL - RuBisCo|rbcS - RuBisCo|Rhodomicrobium", Gene))%>%
  mutate(gene_label=if_else(grepl("Rhodomicrobium_umbrella|Rhodomicrobium_like", gene_label), "Rhodomicrobium_mmoX - (1/3)", gene_label))%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Metabolic_step=gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step))%>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(., aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

pl_simple


tree_met <-pl_simple  %>% aplot::insert_left(Rhodo_plot, width=2)
tree_met_heat<-tree_met %>% aplot::insert_right(heat, width=3.5)


#ggsave("./output/metabolism/Rhodo_tree_simpmet_Sylph_drep.png",tree_met_heat, height=6, width=14)



########################### Trying to add to the simple plot #######################
KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "CO and\nH2", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()


pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(!Module=="Oxi. Phos. Alternative")%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value="transparent", values=c(
    "GSH"="#5F9EA0",
    "H4F"="#C7966D",
    "H4MPT"="green",
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
    "Glycolysis"="#1237B7",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    # "Cyt. aa3-600 menaquinol"="darkred",
    "High affin. cyt c bb3-type"="darkblue"))+
  facet_nested(. ~factor(Metabolic_step, levels=c("From.\noxi","Complex\nIV", "CO and\nH2", "N", "S\noxi./ red.", "S\noxi.","PHB")), scales = "free", space = "free") +
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

options("aplot_guides" = "keep")
#pl_simple_2

tree_met <-pl_simple  %>% aplot::insert_left(Rhodo_plot, width=5.8)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=5.5)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=8.6)


ggsave("./output/metabolism/Rhodo_met_condense_Sylph_drep_25_11_10.png",tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 186,
       dpi=300)



ggsave("./output/metabolism/Rhodo_met_condense_Sylph_drep_25_11_10.svg",tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 186,
       dpi=300)



bar_plot <- sylph %>%
  filter(user_genome %in% Rhodo_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_nested(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs, mires\nand fens", "Freshwater\nsediment","Dunes","Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

library(patchwork)
combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))


ggsave("./output/metabolism/Rhodo_Sylph_bar_25_11_10.svg",combined_plot,
       units = c("mm"),
       height = 110,
       width = 86,
       dpi=300)


#################################################################################################
#################################################################################################
######################      From here, we are working with Methylocella      ####################
#################################################################################################
#################################################################################################





# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocella")%>%
  pull(tip.label)

# remove outgroup,I use this instead
Methylocella_tree <- ape::keep.tip(tr, tip = outgroup_tips)


Methylocella_plot<-ggtree(Methylocella_tree) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.1, size = 2.1,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  # geom_polygon(aes(fill = type, x = 0, y = 0)) +
  theme(legend.position = c(0.1, 0.83))+
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


Methylocella_plot <- Methylocella_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=1.2) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.83),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 2)
  ))


ax<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(y)

Methylocella_plot<-Methylocella_plot + geom_segment(aes(x = ax, y = ay, xend = -0.06, yend = ay),
                                                    arrow = arrow(length = unit(0.2, "cm"), type="closed"), lwd=0.2)+
  annotate("text", x = -0.03, y = ay+1.5, label = "Outgroup", size=2.5)+
  theme(legend.position = c(0.15, 0.83)) +
  geom_treescale(x=-0.019, y=22, label="Tree scale", fontsize =2.7)

Methylocella_plot 



pl_simple <- KEGG %>% 
  filter(genome %in% Methylocella_tree$tip.label)%>%
  filter(presence==1)%>%
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo") | Metabolic_step %in% c("Custom\nHMM"))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Metabolic_step=gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step))%>%
  filter(!Gene %in% c("Methyloferula_Methylovirgula", "Beijerinckiaceae_pxmA", "Methylovulum_2", "Homologous_Methylocella"))%>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  ggplot(., aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = Type), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))


# tree_met_simple<-heat %>%aplot::insert_left(pl_simple, width=0.13)
# tree_met_simple_heat<-tree_met_simple %>% aplot::insert_left(Methylocella_plot, width=0.4)
# 
# 
# ggsave("./output/metabolism/Methylocella_tree_simpmet_Sylph_drep_2.png",tree_met_simple_heat, height=6, width=14)
# 


heat <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order, ordered = TRUE)) %>%
 # mutate(mfd_hab2 = factor(mfd_hab2, levels = levels(hab2_sort_order), ordered = TRUE))%>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #facet_nested(~ hab1_label, scales = "free", space = "free") +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))



bar_plot <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_nested(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs, mires\nand fens", "Freshwater\nsediment","Dunes","Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
    scale_y_continuous(expand = c(0,0)) 



combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))
combined_plot



#ggsave("./output/metabolism/Methylocella_with_bar_drep_sylph_DMS.png",combined_plot, height=6, width=8, dpi=800)




################### Making the tree for paper #############################
# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocella")%>%
  filter(!label_3=="Methylocella MAG_10")%>%  ### two members of the same species here, quantified as 1 species in sylph
  pull(tip.label)


tax2<-tax%>%
  mutate(clade = if_else(between(genome_number, 1, 23), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 24, 30), "clade_2", clade))%>%
  mutate(label_3=gsub("Methylocella", "Methylocapsa", label_3))%>%
  mutate(label_3=gsub("Ca. Methyloaffinis", "Ca. M.", label_3))%>%
  mutate(label_3=gsub("Methylocapsa sp003162995", "M. sp003162995", label_3))

tax2_print_phylo<-tax2%>%
  filter(Genus=="g__Methylocella")%>%
  select(tip.label, user_genome, label_3, clade)

#readr::write_csv(tax2_print_phylo, "/home/bio.aau.dk/vj52ou/data/MFD/annotation_combined/pmoA_methylocapsa/pmoA_v2/mags_to_contigs_clades.csv", col_names = FALSE)


# remove outgroup,I use this instead
Methylocella_tree <- ape::keep.tip(tr, tip = outgroup_tips)
Methylocella_plot<-ggtree(Methylocella_tree, linewidth=0.3) %<+% tax2 + 
  geom_tiplab(aes(label = label_3, fill = clade, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("clade_1" = "#B3A2C7", "clade_2"="#8EB4E3", "other"="grey100"),guide = "none")+
  theme(legend.position = "none")+
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


Methylocella_plot <- Methylocella_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  theme(legend.position = c(0.15, 0.87),
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

Methylocella_plot <- Methylocella_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3,
               inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)


Methylocella_plot 


options("aplot_guides" = "keep")

heat <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

bar_plot <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_nested(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs, mires\nand fens", "Freshwater\nsediment","Dunes","Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

pl_simple <- KEGG %>%
  filter(genome %in% Methylocella_tree$tip.label) %>%
  #filter(presence==1) %>%
  mutate(presence = if_else(presence == 0, NA, presence)) %>%
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  #mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step)) %>%
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




KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()


pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Methylocella_tree$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "Denitrification"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative", "Sulphur"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB"))%>%
  mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "Terminal\noxidase")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))




#setdiff(Methylocella_plot$data$label, pl_simple$data$genome)


tree_met <-pl_simple  %>% aplot::insert_left(Methylocella_plot, width=2.12)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=0.77)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.6)


ggsave("./output/metabolism/Methylocapsa_met_condense_Sylph_drep_25_10_30.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 128,
       width = 186,
       dpi=300)




library(patchwork)
combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))

# ggsave("./output/metabolism/Methylocapsa_with_bar_drep_sylph_25_10_25.svg",combined_plot,
#        units = c("mm"),
#        height = 120,
#        width = 93,
#        dpi=300)




#################
################# Making agri comparison
#################



heat_all <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater", "Fields"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
 # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order, ordered = TRUE)) %>%
  # mutate(mfd_hab2 = factor(mfd_hab2, levels = levels(hab2_sort_order), ordered = TRUE))%>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #facet_nested(~ hab1_label, scales = "free", space = "free") +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests", "Fields")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))





tree_heat <-heat_all  %>% aplot::insert_left(Methylocella_plot, width=0.2)
tree_heat










#################################################################################################
###################### From here, we are working with Methylocystis ############################
#################################################################################################


###############

# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocystis")%>%
  filter(tip.label %in% unique(KEGG$genome))%>%
  pull(tip.label)

tax2<-tax%>%
  mutate(clade = if_else(between(genome_number, 1, 25), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 27, 31), "clade_2", clade))%>%
  mutate(clade = if_else(between(genome_number, 37, 44), "clade_3", clade))

# remove outgroup,I use this instead
Methylocystis_tree <- ape::keep.tip(tr, tip = outgroup_tips)


Methylocystis_plot<-ggtree(Methylocystis_tree, linewidth=0.3) %<+% tax2 + 
  geom_tiplab(aes(label = label_3, fill = clade, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(
    values = c("clade_3" = "#C9BFE3", "clade_2" = "#F0D8C0", "clade_1" = "lightblue", "other" = "grey100"),guide = "none")+
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


Methylocystis_plot <- Methylocystis_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 


ax<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(y)


Methylocystis_plot <- Methylocystis_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.3,
    inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)

Methylocystis_plot 




######## Plotting the heatmap ######

options("aplot_guides" = "keep")


heat <- sylph %>%
  filter(user_genome %in% Methylocystis_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))




bar_plot <- sylph %>%
  filter(user_genome %in% Methylocystis_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
   # legend.position = "bottom",
    legend.position = "none",
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.07, "cm"),      # space between rows
    legend.spacing.x = unit(0.07, "cm"),      # space between columns
    legend.title = element_blank(),
    legend.text = element_text(size=5, margin = margin(l = 0.4)),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    panel.spacing = unit(0.03, "cm", data = NULL))+
    scale_y_continuous(expand = c(0,0)) 




#

# heat_raster<-ggrastr::rasterise(heat, layers="Tile", dpi=300)
# bar_raster<-ggrastr::rasterise(bar_plot, layers="Tile", dpi=300)
#library(patchwork)


combined_plot <- heat / bar_plot + patchwork::plot_layout(heights = c(20, 0.5))


# ggsave("./output/metabolism/Methylocystis_with_bar_drep_sylph_25_10_25.svg",combined_plot,
#        units = c("mm"),
#        height = 130,
#        width = 100,
#        dpi=300)



########################## importing pmoA2 line ####################################

otu_merged<-readRDS("output/merged_otu_table_24_06_24.rds")
otu_pmoA2<-otu_merged%>%  select(
  matches("pmoA2"),
  matches("Seq|mfd")
)

otu_pmoA2_long<-pivot_longer(otu_pmoA2, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  mutate(type=if_else(grepl("pxmA", Tax), "pxmA", "pmoA"))%>%
  mutate(type=if_else(grepl("mmoX", Tax), "mmoX", type))%>%
  mutate(type=if_else(grepl("pmoA2", Tax), "pmoA2", type))%>%
  mutate(type=if_else(grepl("Nevskiales_Macondimonas|Binatales|Rhodopila", Tax), "Put. pmoA/mmoX", type))

otu_pmoA2_long<-otu_pmoA2_long %>%
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

hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater")

otu_pmoA2_long<-otu_pmoA2_long%>%
  filter(mfd_hab1 %in% hab_filter)%>%
  filter(!Tax_short %in% c("Methylosinus_pmoA2", "Methylocystis_pmoA2_1"))%>%
  mutate(Tax_short = factor(
    Tax_short,
    levels = c("Methylocystis_pmoA2_3", "Methylocystis_pmoA2_2", "Methylosinus_Methylocystis_pmoA2"),ordered = TRUE))


heat_pmoA2 <- otu_pmoA2_long %>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = Tax_short, x = SeqId, fill = `RPKM`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
       # axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))



combined_plot <- heat / heat_pmoA2 / bar_plot + patchwork::plot_layout(heights = c(20, 1.8, 0.5))
#combined_plot


ggsave("./output/metabolism/Methylocystis_with_bar_pmoA2_drep_sylph_25_12_04.svg",combined_plot,
       units = c("mm"),
       height = 130,
       width = 130,
       dpi=300)


################### Making a dataframe with summarised metabolism ###############

KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()

############################################################################################
############### Methylocystis condense met similar to the Methylocella DMS #################
############################################################################################




genes_keep<-KEGG%>%filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Metabolic_step %in% c("Custom\nHMM"))%>%
  group_by(Gene)%>%
  filter(presence==1)%>%
  select(Gene)%>%
  distinct()%>%
  filter(!Gene%in%c("Homologous_MO", "Rhodomicrobium_umbrella"))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))


KOs<-readxl::read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)


KOs<-KOs%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA1_", Gene), "Methylocystis_pmoA1", Gene))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))

library(tidyverse)

CM <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step)) %>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA1_", Gene), "Methylocystis_pmoA1", Gene))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pxmA", "Beijerinckiaceae_pxmA", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA1_1", "Methylosinus_Methylocystis_pmoA1", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA2_1", "Methylosinus_Methylocystis_pmoA2", Gene))%>%
  filter(Gene %in% genes_keep$Gene | Gene %in% c("mxaF", "mxaI", "xoxF"))%>% #, "rbcL - RuBisCo", "rbcS - RuBisCo"))%>% (no Rubisco, have checked)
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuBis\nCO")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold"),
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

#CM

pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "Denitrification"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative", "Sulphur"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB"))%>%
  mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "Terminal\noxidase")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold"),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

#pl_simple_2


setdiff(pl_simple_2$data$genome, Methylocystis_tree$tip.label)

tree_met <-CM  %>% aplot::insert_left(Methylocystis_plot, width=2.1)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=0.85)

#heat_raster<-ggrastr::rasterise(heat, layers="Tile", dpi=400)

tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.8)
#tree_met_heat_2

#ggsave("./output/metabolism/Methylocystis_met_condense_Sylph_drep_25_06_12.png",tree_met_heat_2, height=7, width=15, dpi=800)



ggsave("./output/metabolism/Methylocystis_met_condense_Sylph_drep_25_10_24.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 190,
       dpi=300)

###





###################################################
################ Looking at HAO ###################


CM <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  #filter(Module == 'C1') %>% 
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo", "hao") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step)) %>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA1_", Gene), "Methylocystis_pmoA1", Gene))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pxmA", "Beijerinckiaceae_pxmA", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA1_1", "Methylosinus_Methylocystis_pmoA1", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA2_1", "Methylosinus_Methylocystis_pmoA2", Gene))%>%
  filter(Gene %in% genes_keep$Gene | Gene %in% c("mxaF", "mxaI", "xoxF", "hao"))%>% #, "rbcL - RuBisCo", "rbcS - RuBisCo"))%>% (no Rubisco, have checked)
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "Nitrification")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors, "Nitrification"="#EBBC3E")) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold"),
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


















########### Checking where the F the recovered MAGs can be found ##########

heat <- sylph %>%
  filter(user_genome %in% Methylocystis_tree$tip.label)%>%
  #filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  #filter(!mfd_hab2=="Mire")%>%
  #filter(!mfd_sampletype=="Water")%>%
 # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
#  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))



pl_simple <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  #filter(Metabolic_step %in% c("Methanol\noxidation")| gene_label %in% c("pmoA - (1/3)", "pxmA - (1/3)", "mmoX - (1/3)"))%>%
  filter(Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation"))%>%
  #filter(type %in% c('LR-MAG', 'GTDB'))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(., aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        #axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))


genes_keep<-KEGG%>%filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(gene_label)%>%
  filter(presence==1)%>%
  select(gene_label)%>%
  distinct()


pl_simple_v2 <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Metabolic_step=="Custom\nHMM")%>%
  filter(gene_label %in% genes_keep$gene_label)%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(Gene = factor(Gene, levels=unique(KOs$Gene))) %>%
  ggplot(., aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.1 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        #axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))



tree_met_simple <-pl_simple_v2  %>% aplot::insert_left(Methylocystis_plot, width=3)
tree_met_simple_heat<-tree_met_simple %>% aplot::insert_right(heat, width=5)
tree_met_simple_heat














#############################################################
###############  Gammaproteobacteria   ######################
#############################################################


# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Class=="c__Gammaproteobacteria")%>%
  pull(tip.label)

Gamma_tree <- ape::keep.tip(pr, tip = outgroup_tips)
Gamma_plot<-ggtree(Gamma_tree, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"),guide = "none")+
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



# remove outgroup,I use this instead


# put on bootstrap values
df1 <- Gamma_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


Gamma_plot <- Gamma_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=1) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.83),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 2)
  ))


ax<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(y)

Gamma_plot<-Gamma_plot + geom_segment(aes(x = ax, y = ay, xend = -0.05, yend = ay),
                                      arrow = arrow(length = unit(0.2, "cm"), type="closed"), lwd=0.2)+
  annotate("text", x = -0.025, y = ay+1.5, label = "Outgroup", size=2.5)+
  theme(legend.position = c(0.15, 0.83)) +
  geom_treescale(x=0.035, y=9, label="Tree scale", fontsize =2.5)

Gamma_plot 




heat <- sylph %>%
  filter(user_genome %in% Gamma_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
 # mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


tree_sylph<-heat%>% aplot::insert_left(Gamma_plot, width=0.7)

ggsave("./output/genome_quantification/gammas_sylph_gtdb.png", tree_sylph, height=20, width=15)


############################# Okay, clearly, I need to only select some of the interesting genomes and potentially leave the remeaning in supplementary ######



Gamma_plot <- ggtree(Gamma_tree) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align = TRUE, offset = 0.1, size = 2.1, label.size = 0.0,
              hjust = 1, linesize = 0.2, alpha = 1) + 
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0")) + 
  theme(legend.position = c(0.15, 0.83)) + 
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(aes(label = node), hjust = -0.5, size = 2, color = "black")  # Add node IDs

Gamma_plot

tree_sylph<-heat%>% aplot::insert_left(Gamma_plot, width=0.7)



#ggsave("./output/genome_quantification/gammas_sylph_gtdb_temp.png", tree_sylph, height=20, width=15)


# Remove node 168,. 220
# Single tips: 13, 19
pruned_tips_f <- reduce(offspring(Gamma_tree, c("168", "220", "228", "216", "196", "191", "151", "177")), c)
Single_tips<-c("13", "19", "31", "33", "30", "35", "38", "39", "54", "51", "52", "53", "4")
pruned_tips<-Gamma_plot$data%>%filter(isTip==TRUE)%>%filter(node %in% union(pruned_tips_f, Single_tips))%>%pull(label)

extra <- c("GCA_026399025.1", "GCF_002934365.1",
           "GCF_030680395.1", "GCA_028712225.1","GCF_002005105.1",
           "GCA_024636875.1", "GCA_016715325.1", "GCA_030679035.1", "GCA_028707085.1")

pruned_tips <- unique(c(pruned_tips, extra))


# remove outgroup,I use this instead
Gamma_tree_pruned <- ape::drop.tip(Gamma_tree, tip = pruned_tips)



# Recreate the plot with the pruned tree
Gamma_plot_pruned <- ggtree(Gamma_tree_pruned) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align = TRUE, offset = 0.1, size = 2.1, label.size = 0.0,
              hjust = 1, linesize = 0.2, alpha = 1) + 
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0")) + 
  theme(legend.position = c(0.15, 0.83)) + 
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(aes(label = node), hjust = -0.5, size = 2, color = "black")  # Add node IDs

# Display the plot
Gamma_plot_pruned


setdiff(Gamma_plot$data$label, Gamma_plot_pruned$data$label)


tree_sylph<-heat%>% aplot::insert_left(Gamma_plot_pruned, width=0.7)


Gamma_plot_pruned <- ggtree(Gamma_tree_pruned, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align=TRUE, offset = 0.56, size = 5/.pt, hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0"),guide = "none") + 
  scale_x_discrete(expand = c(0, 0)) 



# put on bootstrap values
df1 <- Gamma_tree_pruned %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Gamma_plot_pruned <- Gamma_plot_pruned %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 

ax<-Gamma_plot_pruned$data%>%filter(node==min(parent))%>%pull(x)
ay<-Gamma_plot_pruned$data%>%filter(node==min(parent))%>%pull(y)

Gamma_plot_pruned<-Gamma_plot_pruned+
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
   aes(x = x, y = y, xend = xend, yend = yend),
   linewidth = 0.3,
   inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)

Gamma_plot_pruned 

# 
# tree_sylph<-heat%>% aplot::insert_left(Gamma_plot_pruned, width=0.7)
# 
# ggsave("./output/genome_quantification/gammas_sylph_gtdb_temp2.png", tree_sylph, height=10, width=15)
# 
# 

heat <- sylph %>%
  filter(user_genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Bogs\nand fens", "Freshwater\nsediment", "Grassland\nformations", "Forests", "Green-\nspaces")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))





################################################################################
########################### Nice, then it is metabolism time ###################
################################################################################



genes_keep<-KEGG%>%filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(Metabolic_step %in% c("Custom\nHMM"))%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  group_by(Gene)%>%
  filter(presence==1)%>%
  select(Gene)%>%
  distinct()%>%
  filter(!Gene%in%c("Homologous_MO", "Rhodomicrobium_umbrella"))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))


KOs<-readxl::read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)


KOs<-KOs%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))
  



CM <- KEGG %>% 
  filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Type %in% c("RuMP cycle") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  mutate(Gene=if_else(Type=="RuMP cycle", Gene_collapsed, Gene))%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  filter(Gene %in% genes_keep$Gene | Gene %in% c("mxaF", "mxaI", "xoxF", "hxlB/hps-phi", "hxlA/fae-hps/hps-phi"))%>% 
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuMP", Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuMP")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

CM



KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "DNRA / denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step),
         Metabolic_step=gsub("DNRA / denitrification", "N", Metabolic_step))%>%
  ungroup()



pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "aprAB", "DsrABL", "DsrC", "dsrMKJOP", "dsrEFH"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB", "Complex\nIV"))%>%
 # mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("N\nfix", "N", Metabolic_step))%>%
  mutate(Metabolic_step=if_else(grepl("S\noxi.", Metabolic_step), "S oxi.\n/red.", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "N", "Terminal\noxidase", "S oxi.\n/red.")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "Assimilatory sulphate reduction"="#0DB7A8",
    "sat"="#6C961E",
    "SOX"="#92C76C"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

pl_simple_2


setdiff(Gamma_plot_pruned$data$label, pl_simple_2$data$genome)

tree_met <-CM  %>% aplot::insert_left(Gamma_plot_pruned, width=3.1)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=1)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.8)
#tree_met_heat_2


ggsave("./output/metabolism/Gammas_met_condense_Sylph_drep_25_11_10.png",
       tree_met_heat_2,
       units = c("mm"),
       height = 160,
       width = 190,
       dpi=300)


ggsave("./output/metabolism/Gammas_met_condense_Sylph_drep_25_11_10.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 160,
       width = 190,
       dpi=300)




#


bar_plot <- sylph %>%
  filter(user_genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ factor(hab1_label, levels=c("Bogs\nand fens", "Freshwater\nsediment", "Grassland\nformations", "Forests", "Green-\nspaces")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))


ggsave("./output/metabolism/Gamma_Sylph_bar_25_11_10.svg",combined_plot,
       units = c("mm"),
       height = 110,
       width = 86,
       dpi=300)



#################### I want to look at the full alternative metabolism ##############




























#############################################################
#######  Methylomirabilaceae and Verrucomicrobia ############
#############################################################


############ Something here is off haha


# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus %in% c("g__Methylomirabilis", "g__Methylacidimicrobium", "g__Methylacidiphilum", "g__Methylacidithermus"))%>%
  pull(tip.label)

# remove outgroup,I use this instead
mira_acidi_tree <- ape::keep.tip(pr, tip = outgroup_tips)
mira_acidi_plot<-ggtree(mira_acidi_tree) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.1, size = 2.1,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"))+
  # geom_polygon(aes(fill = type, x = 0, y = 0)) +
  theme(legend.position = c(0.15, 0.83))+
  scale_x_discrete(expand = c(0,0))

mira_acidi_plot

# put on bootstrap values
df1 <- mira_acidi_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


mira_acidi_plot <- mira_acidi_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=1) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.83),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 2)
  ))


ax<-mira_acidi_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-mira_acidi_plot$data%>%filter(node==min(parent))%>%pull(y)

mira_acidi_plot<-mira_acidi_plot + geom_segment(aes(x = ax, y = ay, xend = -0.05, yend = ay),
                                      arrow = arrow(length = unit(0.2, "cm"), type="closed"), lwd=0.2)+
  annotate("text", x = -0.025, y = ay+1.5, label = "Outgroup", size=2.5)+
  theme(legend.position = c(0.15, 0.83)) +
  geom_treescale(x=0.035, y=9, label="Tree scale", fontsize =2.5)

mira_acidi_plot 




heat <- sylph %>%
  filter(user_genome %in% mira_acidi_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


tree_sylph<-heat%>% aplot::insert_left(mira_acidi_plot, width=0.7)


### Not aligning at all :)

#ggsave("./output/genome_quantification/gammas_sylph_gtdb.png", tree_sylph, height=20, width=15)


























#############################################################
#################### To get the legend ######################
#############################################################


# Define the list of habitats

sylph %>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% 
  select(hab1_label, mfd_hab2)%>%
  arrange(hab1_label)%>%
  distinct()

# 
ISME_colours<-c( 
  `Calcareous fens` = "#cdad00",
  `Sphagnum acid bogs` = "#32849f",
  `Wet thicket` = "#7301A8FF",
  `Sea dunes`="#cdad00",
  `Inland dunes`="#32849f",
  `Deciduous trees` = "#0D0887FF" ,
  `Forests no MFDO2` = "#32849f",
  `Non-native trees` = "#cd2626",
  `Coniferous forest` = "#cc33ff",
  Beech = "#8CB88F",
  `Alluvial woodland` = "#33ccff",
  Oak = "#ff9999",
  Willow= "#7301A8FF",
  `Bog woodland`  = "#cdad00",
  Spruce = "#000000",
  `Urban enclosed water`="#cd2626",
  `Rainwater basin`= "#8CB88F" ,
  `Running freshwater`= "#ff9999",
  `Standing freshwater, lake`="#7301A8FF",
  `Standing freshwater, other` = "#32849f",
  `Semi-natural dry grasslands` = "#32849f",
  `Semi-natural humid meadows` = "#cdad00",
  `Natural grasslands` = "#7301A8FF",          
  Parks = "#cdad00",
  Other = "#32849f",
  `Temperate heath` = "#cdad00",
  `Sclerophyllous scrub` = "#cd2626")


saveRDS(ISME_colours, "./palette_mfd_hab2_ISME.rds")


habitats <- c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater")

# Loop through each habitat
for (habitat in habitats) {
  # Filter the data for the current habitat
  plot_data <- sylph %>%
    filter(bin %in% Methylocella_tree$tip.label) %>%
    filter(mfd_hab1 == habitat) %>%
    mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
    mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2)))%>%
    filter(!mfd_hab2=="Mire")%>%
    filter(!mfd_sampletype=="Water")
  
  # Create the plot
  bar_plot <- ggplot(plot_data, aes(x = SeqId, fill = mfd_hab2)) +
    geom_tile(aes(y = 1)) +
    facet_nested(. ~ hab1_label, scales = "free", space = "free") +
    scale_fill_manual(values = ISME_colours) +
    labs(fill = habitat) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0.1, "lines"),
      legend.text = element_text(size=11),
      legend.title = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.5, "cm")
    ) +
    scale_y_continuous(expand = c(0, 0))+
    guides(fill = guide_legend(nrow = 3))
  
  legend <- ggpubr::get_legend(bar_plot)
  leg<-ggpubr::as_ggplot(legend)
  assign(paste0("leg_", gsub(", ", "_", gsub(" ", "_", habitat))), ggpubr::as_ggplot(legend))
  
  # Save the plot as SVG
  #ggsave(filename = paste0("./output/legend_habitat/bar_plot_", gsub(", ", "_", gsub(" ", "_", habitat)), ".svg"), plot = leg)
}

leg_comb<-`leg_Bogs,_mires_and_fens` + leg_Dunes + leg_Forests + leg_Freshwater + leg_Grassland_formations + leg_Temperate_heath_and_scrub + plot_layout()

ggsave("./output/legend_habitat/combined_DMS.svg", height=5, width=20)




############################################################################################################
######################################## Getting the sylph levels ##########################################
############################################################################################################
hab2_sort_order<-c("Temperate heath", "Sclerophyllous scrub", "Sphagnum acid bogs","Wet thicket", "Calcareous fens", "Standing freshwater, lake",
                   "Standing freshwater, other", "Running freshwater", "Rainwater basin", "Urban enclosed water",
                   "Sea dunes", "Inland dunes",
                   "Semi-natural humid meadows", "Semi-natural dry grasslands", "Natural grasslands",
                   "Beech", "Deciduous trees", "Non-native trees", "Oak", "Coniferous forest", "Forests no MFDO2","Bog woodland", "Alluvial woodland", "Parks", "Other")

OTU_wide<-sylph%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% 
  filter(mfd_hab2 %in% hab2_sort_order)%>%
  pivot_wider(names_from = bin, values_from = `Taxonomic_abundance`, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)


hab2_sort_order<-factor(hab2_sort_order)


OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
groups<-sylph%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  arrange(mfd_hab2)%>%
  select(mfd_hab2)%>%
  distinct()%>%
  pull(mfd_hab2)


for (i in seq_along(groups)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- sylph %>%
    filter(mfd_hab2 %in% groups[i]) %>%
    pivot_wider(names_from = bin, values_from = `Taxonomic_abundance`, id_cols = c(SeqId)) %>% ## Change to tax curated
    as.data.frame(.) %>%
    filter(rowSums(select(., -1)) > 0)
  
  # Assign the result to the dynamically named variable
  assign(gr_name, current_group)
  
  # Continue with the rest of your code using the dynamic variable name
  OTU_sums <- get(gr_name)[,-1]
  row.names(OTU_sums) <- get(gr_name)$SeqId
  
  bray_curtis_dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
  hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")
  ward_dendrogram <- as.dendrogram(hclust_ward)
  ward_order <- order.dendrogram(ward_dendrogram)
  assign(paste0("levels_", gr_name), hclust_ward$labels[order.dendrogram(ward_dendrogram)])
  
  
}

levels_hab2 <- c(
  levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22, levels_gr_23, levels_gr_24, levels_gr_25)

saveRDS(levels_hab2, "./output/genome_quantification/levels_hab2_sylph_DMS.rds")











####################################################################################

##################### Trying to play with the bar - below not working #############

####################################################################################

# 
# x<-ggplot_build(bar_plot)
# xx<-x$data[[1]]%>%select(fill, xmax)%>%distinct()
# yy<-unique(x$plot$data$SeqId)
# 
# 
# # Convert y to a named vector for easy access in ggplot
# color_mapping <- setNames(xx$fill, yy)
# 
# hepl<-unlist(color_mapping)
# df <- tibble(mfd_hab2 = names(palette_mfd_hab2), color = unlist(palette_mfd_hab2))
# 
# col<-merge(genom_abun, df, by="mfd_hab2")%>%
#   select(mfd_hab2, color, SeqId)%>%
#   distinct()
# 
# 
# x<-genom_abun %>%
#   filter(bin %in% Methylocella_tree$tip.label)%>%
#   filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
#   filter(!mfd_hab2=="Mire")%>%
#   filter(!mfd_sampletype=="Water")%>%
#   pull(SeqId)
# 
# 
# col_v <- col %>%
#   filter(SeqId %in% x)%>%
#   mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
#   arrange(SeqId)
# 
# # Extract the 'color' column and set its names to the corresponding 'SeqId' values
# color_vector <- col_v %>%
#   pull(color)
# names(color_vector) <- col_v %>%
#   pull(SeqId)
# 
# # Display the named vector
# color_vector



# df <- tibble(mfd_hab2 = names(palette_mfd_hab2), color = unlist(palette_mfd_hab2))
# 
# col<-merge(genom_abun, df, by="mfd_hab2")%>%
#   select(mfd_hab2, color, SeqId)%>%
#   distinct()
# 
# col_v <- col %>%
#   mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
#   arrange(SeqId)
# 
# # Extract the 'color' column and set its names to the corresponding 'SeqId' values
# color_vector <- col_v %>%
#   pull(color)
# names(color_vector) <- col_v %>%
#   pull(SeqId)
# 
# # Display the named vector
# color_vector


# install.packages("ggtext")
# library(ggtext)
# 
# heat <- genom_abun %>%
#   filter(bin %in% Methylocella_tree$tip.label)%>%
#   filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
#   filter(!mfd_hab2=="Mire")%>%
#   filter(!mfd_sampletype=="Water")%>%
#   mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
#   ggplot(aes(y = bin, x = SeqId, fill = `Relative Abundance (%)`)) +
#   geom_tile() +
#   scale_fill_viridis_c(name = "Relative\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 0.4), na.value = "#F8E622")+
#   facet_nested(~ hab1_label, scales = "free", space = "free") +
#   #facet_nested(~ hab1_label+mfd_hab2, scales = "free", space = "free") +
#   #theme_minimal() +
#   theme(#axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_blank(),
#         #axis.text.y = element_text(size = 8.2),
#         axis.text.x = element_markdown(colour = color_mapping, size=0.01, angle = 90),
#         legend.position = "right",
#         strip.text.x = element_text(size = 8),  # Size for x-axis facet labels
#         #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
#         strip.background = element_rect(fill = "grey90", color = "white"),
#         text = element_text(family = "Arial"),
#         plot.background = element_rect(fill = "transparent"),
#         panel.spacing = unit(0.03, "cm", data = NULL))+
#   scale_y_discrete(expand = c(0,0))+
#   scale_x_discrete(expand = c(0,0))

ggsave("./output/metabolism/test.png",heat, height=10, width=26)

#heat %>% aplot::insert_bottom(bar_plot, height = 0.02)

ggarrange(heat, bar_plot,
          ncol = 1, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = F)





heat1 <- genom_abun %>%
  filter(bin %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens"))%>%
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = bin, x = SeqId, fill = `Relative Abundance (%)`)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 0.4), na.value = "#F8E622")+
  # facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        strip.text.x = element_text(size = 8),  # Size for x-axis facet labels
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        strip.background = element_rect(fill = "grey90", color = "white"),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


heat1
#unique(heat$data$SeqId)

bar_plot1 <- genom_abun %>%
  filter(bin %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens"))%>% filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  # facet_nested(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.spacing = unit(0.1, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) 


t1<-heat1  %>% aplot::insert_bottom(bar_plot1, height=0.02)
t1
#





heat2 <- genom_abun %>%
  filter(bin %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Forests"))%>%
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = bin, x = SeqId, fill = `Relative Abundance (%)`)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 0.4), na.value = "#F8E622")+
  # facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        strip.text.x = element_text(size = 8),  # Size for x-axis facet labels
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        strip.background = element_rect(fill = "grey90", color = "white"),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


heat2
#unique(heat$data$SeqId)

bar_plot2 <- genom_abun %>%
  filter(bin %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Forests"))%>% filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  # facet_nested(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.spacing = unit(0.1, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) 


t2<-heat2  %>% aplot::insert_bottom(bar_plot2, height=0.02)
t2

t3<-t1  %>% aplot::insert_left(t2)

combined_plot <- ggarrange(t2, t1)  + plot_layout(widths = c(20, 20))

t4<-bar_plot2%>%aplot::insert_left(bar_plot1)
t4






####################################################################################

##################### If I want to add some ANI #############################

####################################################################################

# Look here when I want to label the clades of the tree https://4va.github.io/biodatasci/r-ggtree.html 
fast_ani<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/Rhodomicrobium_24_07_09/ANI_Rhodo_Methylocella/fastANI.tsv", 
                col_names = c("genome", "genome_2", "ANI","aligned_contigs","total_contigs"))%>%
  mutate(genome = gsub("_genomic.fna", "", genome),
         genome = gsub(".*\\/", "", genome),
         genome = gsub(".fa","", genome),
         genome_2 = gsub("_genomic.fna", "", genome_2),
         genome_2 = gsub(".*\\/", "", genome_2),
         genome_2 = gsub(".fa","", genome_2))


ANI_heat <- fast_ani %>%
  filter(genome %in% tree_remove_2$tip.label,
         genome_2 %in% tree_remove_2$tip.label)%>%
  mutate(genome=factor(genome, levels=order, ordered = TRUE))%>%
  ggplot(aes(genome, genome_2)) +
  theme_nothing()+
  geom_tile(aes(fill = ANI)) +
  #geom_text(aes(label = round(ANI, 2)),size=2.2) +
  scale_fill_gradientn(colors = c( "#deede8", "#7dbaa4","#31504f", na.value = "#8c493d")) +
  theme(legend.position = 'bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())
ANI_heat

ANI_heat <- fast_ani %>%
  filter(genome %in% tree_remove_2$tip.label,
         genome_2 %in% tree_remove_2$tip.label) %>%
  mutate(genome = factor(genome, levels = order, ordered = TRUE),
         ANI_category = cut(ANI, breaks = c(75, 80, 85, 90, 95, Inf), 
                            labels = c("75-80", "80-85", "85-90", "90-95", ">95"), 
                            include.lowest = TRUE)) %>%
  ggplot(aes(genome, genome_2)) +
  theme_nothing()+
  geom_tile(aes(fill = ANI_category)) +
  # geom_text(aes(label = round(ANI, 2)), size = 2.2) +
  scale_fill_manual(values = c("75-80" = "#deede8",
                               "80-85" ="#91C5B2",
                               "85-90" = "#59a78b",
                               "90-95" = "#31504f",
                               ">95" = "#8c493d"), 
                    na.value = "#8c493d") +
  theme(legend.position = 'bottom',
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())

ANI_heat

comb_3 <-comb_2  %>% aplot::insert_right(ANI_heat)
comb_3


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



###################################################################################
################################# New Sylph #######################################
###################################################################################



.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggtree)
library(ape)
library(ggh4x)
library(treeio)

options("aplot_guides" = "keep")

setwd("~/scripts/MFD/methanotrophs/R_scripts/")


# MFD_gtdb<-vroom("/home/bio.aau.dk/wz65bi/mfd_sylph_quant/analysis/MFD_drep_gtdb_quant/MFD_drep_gtdb_tax_relative_abundance.tsv", delim = "\t")%>%
#   filter(grepl("g__Methylocella|g__Methylocystis|g__Methylosinus|g__Methylobacter", clade_name))%>%
#   filter(grepl("s__", clade_name))%>%
#   pivot_longer(cols=starts_with("/projects"), names_to="SeqId", values_to = "Taxonomic_abundance", names_prefix="/projects/microflora_danica/data/sequencing/data_flat/trimmed/")
# 
# MFD_gtdb<-MFD_gtdb%>%
#   mutate(SeqId=gsub("_R1.fastq.gz", "", SeqId))%>%
#   separate(clade_name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|")
# 
# MFD_gtdb<-MFD_gtdb%>%distinct()
# 
# saveRDS(MFD_gtdb, "./dataframes/sylph_selection_25_02_12.rds")
MFD_gtdb<-readRDS("./dataframes/sylph_selection_25_02_12.rds")


linkage <- vroom("2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  relocate(SeqId)%>%
  filter(after_total_reads>500000)%>%
  select(SeqId, fieldsample_barcode)

metadata.sub <- readxl::read_excel('2024-02-13_mfd_db.xlsx') %>%
  left_join(linkage)%>%
  filter(SeqId %in% linkage$SeqId)%>%
  relocate(SeqId)%>%
  select(SeqId:mfd_hab3, cell.10km, cell.1km) %>%
  mutate(complex = str_c(mfd_sampletype, mfd_areatype, mfd_hab1, sep = ", "))

MFD_gtdb<-MFD_gtdb%>%left_join(metadata.sub)

MFD_gtdb<-MFD_gtdb%>%
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))
  
  




sylph_heat<-MFD_gtdb%>%
  filter(!is.na(fieldsample_barcode))%>%
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
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))%>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))


levels_hab2<-readRDS("./output/genome_quantification/levels_hab2_sylph_DMS.rds")
palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")


heat <- sylph_heat %>%
  filter(Genus %in% c("g__Methylocella", "g__Methylocystis"))%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes","Freshwater", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(y = Species, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_nested(Genus~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 8, face="bold"),
        strip.text.y = element_text(angle=0),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#heat


bar_plot <- sylph_heat %>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes","Freshwater", "Forests"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_nested(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour="black", fill=NA) ,
    panel.spacing = unit(0.03, "cm", data = NULL),
    panel.background = element_blank()
  )+   scale_y_continuous(expand = c(0,0)) 


#

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.4))
#combined_plot


ggsave("./output/genome_quantification/Sylph_MFD_GTDB_25_02_13.png",combined_plot, height=8, width=15, dpi=300)

