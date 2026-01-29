#!/usr/bin/env Rscript

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
#BiocManager::install("ggtreeExtra")
library(ggnewscale)
library(ggtreeExtra)
library(ggbreak)
library(ggpattern)


setwd("~/scripts/MFD/methanotrophs/R_scripts/")

### Making the table of genomes recovered ####
tax_filter<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/taxonomy_filter_r220.txt", delim = "\t", col_names = "tax")%>%
  mutate(tax=gsub("o__Methylococcales;f__UBA1147;g__UBA1147", "g__UBA1147", tax))%>%pull(tax)

linkage<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/MFD_LR_SR_r220_sp_cluster.tsv", delim="\t", col_names = c("user_genome","classification"))%>%
  mutate(user_genome=gsub(".fa","", user_genome))%>%
  separate(classification, into = c("Domain", "Phylum", "Class","Order", "Family", "Genus", "Species_cluster"), sep = ';', remove=T)%>%
  select(user_genome, Species_cluster)

#quality<-vroom("/projects/microflora_danica/deep_metagenomes/mags_v4/MFD_mags_V4.tsv")%>%
quality<-vroom("/projects/microflora_danica/MFD-LR/mags_v4/MFD_mags_V4.tsv")%>%
  select(bin, MAG_status)%>%rename(user_genome=bin)

quality2<-vroom("/projects/microflora_danica/MFD-SR/results/mags_shallow_all.tsv")%>%
  select(bin, MAG_status)%>%rename(user_genome=bin)
qual<-rbind(quality, quality2)

tax<-readRDS("tax_MAGs_gtdb.rds")%>%
  left_join(qual)%>%
  filter(str_detect(classification, str_c(tax_filter, collapse = "|")))%>%
  left_join(linkage)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'LR-MAG','GTDB'),
         type=if_else(grepl("LIB", user_genome), "SR-MAG", type))%>%
  mutate(novel = if_else(type %in% c("LR-MAG", "SR-MAG") & Species == "s__", "yes", "no"))


df<-tax%>%mutate(Genus=if_else(Genus=="g__", paste0(Family,";",Genus), Genus))%>%
  group_by(Genus)%>% summarise(
    Genomes_recovered = sum(type %in% c("LR-MAG", "SR-MAG")),  # Count of LR-MAG & SR-MAG
    HQ_genomes = sum(type %in% c("LR-MAG", "SR-MAG") & MAG_status=="HQ"),
    Novel_species_recovered = sum(type %in% c("LR-MAG", "SR-MAG") & !is.na(Species_cluster) & novel == "yes"),
    HQ_genomes_novel = sum(type %in% c("LR-MAG", "SR-MAG") & !is.na(Species_cluster) & novel == "yes" & MAG_status=="HQ")# Count with conditions
  )%>%filter(Genomes_recovered>0)



#readr::write_csv(df, "/home/bio.aau.dk/vj52ou/data/MFD/genome_phylogeny/methanotrophs/genomes_recovered_25_02_18.csv", col_names = T)




# Tree file
tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/with_gtdb_v2/MSA.faa.treefile")

tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

tax<-readRDS("tax_MAGs_gtdb.rds")%>%filter(user_genome %in% tree_drep$tip.label)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'LR-MAG','GTDB'),
         type=if_else(grepl("LIB", user_genome), "SR-MAG", type))%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", user_genome ))%>%
  mutate(Family_label=gsub("f__","", Family))


tr <- phytools::midpoint.root(tree_drep)
# plot rooted tree
pr <- ggtree(tr) 
pr

### Now, I need to remove the ones that are duplicates 

tax_unique <- tax %>%
  filter(case_when(
    type == "GTDB" ~ TRUE,
    type %in% c("LR-MAG", "SR-MAG") & Species == "s__" ~ TRUE,
    TRUE ~ FALSE
  ))%>%
  mutate(hit=if_else(type=="GTDB", NA, "yes"))


### Then, I will remove those from the tree
unique_tips <- tax_unique %>%
  pull(user_genome)

# remove outgroup,I use this instead
Unique_tree <- ape::keep.tip(tr, tip = unique_tips)
Unique_plot<-ggtree(Unique_tree,size=0.35, layout="fan") %<+% tax + 
  geom_treescale(x=0.3, y=Unique_tree$Nnode/2, label="Tree scale", fontsize =3.5, width = 0.3, offset.label = 1.5, offset=-4) 



#rotate_tree(Unique_plot, 180)


#### Now I need to make a color scale for the Family ######
unique(tax_unique$Order)
Order_colors = c("#734c78",'#bf4e4e',"#ACABFC90","#826765" ,'#9c760560',"#AFEEEE",'#3559e6', '#b2df8a','#e31a1c50','#418cb2','gold4','#34978F90', 'orange3' )
Order_colors <- rev(Order_colors)

your_family_colors <- c(
  "Methylomonadaceae" = "#9c760560",
  "Methylococcaceae" = "#b2df8a",
  "Beijerinckiaceae" = "#34978F90",
  "Xanthobacteraceae" = "#734c78",
  "Rhodomicrobiaceae"= "#ACABFC90",
  "Methyloligellaceae"="#3559e6",
  "Azospirillaceae"="orange3",
  "Methylomirabilaceae"="#AFEEEE",
  "Methylacidiphilaceae"="#e31a1c50",
  "JACCXJ01"="#418cb2",
  "Methylothermaceae"="#826765",
  "DRLZ01"="gold4",
  "UBA1147"="#bf4e4e"
)


tip_data <- tax_unique %>% select(user_genome, hit) 
p_save <- Unique_plot %<+% tip_data 

# combine
p_savephyl <- p_save +
  geom_tippoint(mapping=aes(color=hit),
                size=2, stroke=0)+
  scale_color_manual(values = '#dd0305',
                     na.value = NA, guide="none")


p_savephyl_fam <- p_savephyl + new_scale_fill() + geom_fruit(
  geom = geom_tile,
  mapping = aes(y=tip_label, group=label,
                fill=Genus),
  #alpha = 0.8,
  width=0.25,
  offset= -0.02 #before -0.03
) #+
scale_fill_manual('Family',
                  values = your_family_colors)

# 
# x<- p_savephyl  +
#   geom_text(aes(x=branch, label=node), size=3,
#             vjust=-.3, color="firebrick")
# 
# Node_plot <- x + new_scale_fill() + geom_fruit(
#   geom = geom_tile,
#   mapping = aes(y=tip_label, group=label,
#                 fill=Genus),
#   #alpha = 0.8,
#   width=0.25,
#   offset= -0.02 #before -0.03
# ) 
# 
# 
# Node_plot<-ggtree(Unique_tree,size=0.35) %<+% tax + 
#   geom_treescale(x=0.3, y=Unique_tree$Nnode/2, label="Tree scale", fontsize =3.5, width = 0.3, offset.label = 1.5, offset=-4) +
#   geom_tippoint(mapping=aes(color=Genus),size=2, stroke=0)+
#   geom_fruit(
#     geom = geom_tile,
#     mapping = aes(y=tip_label, group=label,
#                   fill=Genus),
#     #alpha = 0.8,
#     width=0.25,
#     offset= -0.02 #before -0.03
#   ) 
# 
# ggsave("node.png", Node_plot, width = 40, height=40)
# 



p_savephyl <- p_save +
  geom_tippoint(mapping=aes(color=hit),
                size=2, stroke=0)+
  scale_color_manual(values = '#dd0305',
                     na.value = NA, guide="none")


p_savephyl_fam <- p_savephyl + new_scale_fill() + geom_fruit(
  geom = geom_tile,
  mapping = aes(y=tip_label, group=label,
                fill=Family_label),
  #alpha = 0.8,
  width=0.25,
  offset= -0.02 #before -0.03
) +
  scale_fill_manual('Family',
                    values = your_family_colors)







p_savephyl_fam


dat <- data.frame(
  node = c(969, 899,1019,589,624,888, 694, 720, 750, 862) ,#, 673),
  name = c("Methylocella", "Methylocystis","Rhodomicrobium","Methylobacter_A", "Methylobacter_C","USCg-Taylor", "UBA10906", "Methyloglobulus", "Methylomonas", "Methylumidiphilus"))#, "Methylovulum"))

dat

p_savephyl_fam_annot<-p_savephyl_fam+
  geom_cladelab(data=dat, mapping = aes(node=node,label=name),
                fontsize = 3.3,
                align = TRUE,
                show.legend = FALSE,
                angle = "auto",
                #hjust=0.5,
                hjust = "center",
                barsize = 0.3,
                #  vjust=0.5,
                horizontal =F,offset= 0.12, offset.text=0.08
  )

p_savephyl_fam_annot
#ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_10_29.svg', height=8, width = 12)
#ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_10_29.png', height=8, width = 12)



# Count the number of added genomes per genus
genus_summary <- tax_unique %>%
  group_by(Genus) %>%
  summarize(
    added_genomes = sum(hit == "yes", na.rm = TRUE),   # Number of genomes being added
    gtdb_genomes = sum(type == "GTDB", na.rm = TRUE),  # Number of GTDB genomes per genus
    total_genomes = n()                                # Total number of genomes per genus
  ) %>%
  mutate(
    percent_increase = round((added_genomes / (gtdb_genomes + added_genomes)) * 100, 0) # Percent increase
  ) %>%
  ungroup()

# Merge the genus_summary back to the original tax_unique dataframe
tax_final <- tax_unique %>%
  left_join(genus_summary, by = "Genus")%>%
  select(!c("user_genome", "classification", "Species", "tree_label", "type", "hit"))%>%
  distinct()%>%
  filter(percent_increase>10)

# View the final dataframe
print(tax_final)



################################################################
####### Adding bar-plots with No of genomes recovered ##########
################################################################

tax<-readRDS("tax_MAGs_gtdb.rds")%>%
  left_join(qual)%>%
  filter(str_detect(classification, str_c(tax_filter, collapse = "|")))%>%
  left_join(linkage)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'LR-MAG','GTDB'),
         type=if_else(grepl("LIB", user_genome), "SR-MAG", type))%>%
  mutate(novel = if_else(type %in% c("LR-MAG", "SR-MAG") & Species == "s__", "yes", "no"))


df<-tax%>%mutate(Genus=if_else(Genus=="g__", paste0(Family,";",Genus), Genus))%>%
  group_by(Family,Genus)%>% summarise(
    'Total genomes' = sum(type %in% c("LR-MAG", "SR-MAG")),  # Count of LR-MAG & SR-MAG
    'HQ genomes' = sum(type %in% c("LR-MAG", "SR-MAG") & MAG_status=="HQ"),
    'New Species' = sum(type %in% c("LR-MAG", "SR-MAG") & !is.na(Species_cluster) & novel == "yes"),
    'HQ new Species' = sum(type %in% c("LR-MAG", "SR-MAG") & !is.na(Species_cluster) & novel == "yes" & MAG_status=="HQ")# Count with conditions
  )%>%filter('Total genomes' >0)%>%
  pivot_longer(cols=c("Total genomes", "HQ genomes", "New Species", "HQ new Species"), names_to = "genome_stat", values_to = "count")



p1<-df%>%
  group_by(Genus) %>%
  filter(!Family=="f__Methylomonadaceae")%>%
  filter(any(genome_stat == "New Species" & count > 0)) %>%
  ungroup()%>%
  mutate(Family=gsub("f__", "", Family))%>%
  mutate(Genus=gsub("^g__", "", Genus))%>%
  mutate(genome_stat = fct_relevel(genome_stat,"Total genomes", "New Species","HQ genomes",  "HQ new Species"))%>%
  mutate(Family = factor(Family, levels = c("Methylococcaceae", "JACCXJ01", "Beijerinckiaceae", "Rhodomicrobiaceae", "Methylomirabilaceae")))%>%
  mutate(Genus=gsub("Methylocella", "Methylocapsa", Genus))%>%
  ggplot(., aes(x = Genus, y = count, fill = Family, pattern=genome_stat)) +
  geom_bar_pattern(
    stat = "identity",
    color = "black",
    position = position_dodge(),
    pattern_fill="white",
    width = 0.8,
    pattern_spacing=0.11, pattern_frequency=0.5, pattern_units='cm'
  ) +
  scale_pattern_manual(values = c(
    "Total genomes" = "none", 
    "New Species" = "stripe",              
    "HQ genomes" = "circle",         
    "HQ new Species" = "crosshatch"                  
  )) +
  scale_pattern_type_manual(values=c(NA, NA, NA, 'plain')) +
  geom_text(aes(label = ifelse(count > 0, count, "")), position = position_dodge(width = 0.8), vjust = -0.35, hjust= 0.4, size = 5 / .pt)+
  #geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
  theme_minimal() +
  labs(y = "No. of MAGs") +
  scale_fill_manual(values = your_family_colors) +
  facet_grid(~ Family, scales = "free_x", space="free_x", switch = "x") +  
 # facet_nested(~ Family, scales = "free_x", space="free_x", switch = "x", nest_line = element_line(color="grey20", linewidth = 0.4), resect=unit(2, "pt"),strip = strip_nested(clip = "off")) +
  theme(
    panel.spacing = unit(2, "pt"),
    strip.clip = "off",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, margin = margin(t = -1)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    strip.text.x = element_text(size = 5, margin = margin(t= -0.5)),
    strip.background = element_blank(),
    panel.border = element_blank()
  )+
  scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 5, 10, 20, 50, 90, 130))#+
#  scale_y_break(c(90, 130))
# scale_y_cut(breaks=c(5, 50), which=c(1,2,3), scales=c(1,1,3), space=.1)



p2<-df%>%
  group_by(Genus) %>%
  filter(Family=="f__Methylomonadaceae")%>%
  filter(any(genome_stat == "New Species" & count > 0)) %>%
  ungroup()%>%
  mutate(Family=gsub("f__", "", Family))%>%
  mutate(Genus=gsub("f__", "", Genus))%>%
  mutate(Genus=gsub("^g__", "", Genus))%>%
  mutate(genome_stat = fct_relevel(genome_stat,"Total genomes", "New Species","HQ genomes",  "HQ new Species"))%>%
  ggplot(., aes(x = Genus, y = count, fill = Family, pattern=genome_stat)) +
  geom_bar_pattern(
    stat = "identity",
    color = "black",
    position = position_dodge(),
    pattern_fill="white",
    width = 0.8,
    pattern_spacing=0.11, pattern_frequency=0.5, pattern_units='cm'
  ) +
  scale_pattern_manual(values = c(
    "Total genomes" = "none", 
    "New Species" = "stripe",              
    "HQ genomes" = "circle",         
    "HQ new Species" = "crosshatch"                  
  )) +
  #geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
  labs(y = "No. of MAGs") +
  geom_text(aes(label = ifelse(count > 0, count, "")), position = position_dodge(width = 0.8), vjust = -0.35, hjust=0.4, size = 5 / .pt)+
  theme_minimal() +
  guides(fill = "none", pattern_fill = "none") + 
  scale_fill_manual(values = your_family_colors) +
  facet_grid(~ Family, scales = "free_x", space="free_x", switch = "x") +  
  theme(
    panel.spacing = unit(2, "pt"),
    strip.clip = "off",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5, margin = margin(t = -1)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    strip.text.x = element_text(size = 5, margin = margin(t= 1)),
    strip.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(0, 2, 0, 0)
  )+
  scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 5, 15))#+
#scale_y_break(c(50, 100))
# scale_y_cut(breaks=c(5, 50), which=c(1,2,3), scales=c(1,1,3), space=.1)

# First, recreate the plot but with legend visible
p2_legend <- p2 +
  geom_bar_pattern(
    stat = "identity",
    color = "black",
    position = position_dodge(),
    pattern_fill="white",
    width = 0.2,
    pattern_spacing=0.11, pattern_frequency=0.5, pattern_units='cm', pattern_key_scale_factor=0.9
  ) +
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill = "white",           # 👈 white background behind pattern
        pattern_fill = "white",   # 👈 ensures pattern fill is white
        colour = "black"          # optional: border color
      )
    ),
    fill = "none",               # 👈 hide fill legend
  ) +
  theme(legend.position = "right",
        legend.key.size = unit(0.25, "cm"),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1))


# Extract only the legend
legend_only <- get_legend(p2_legend)
p_leg<-plot_grid(legend_only)




# Top row: p1 and legend side by side
top_row <- p1 + p_leg + plot_layout(widths = c(6, 1))

# Combine with bottom row (p2), using relative heights
combined_plot <- top_row / p2 +
  plot_layout(heights = c(3.5, 1.3))+ #, widths=c(1,1)) +  # p2 is half the height of top_row
  plot_annotation(theme = theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ))



#
ggsave("output/barplot_genomes_25_11_07.png", combined_plot, 
       #width=15, 
       #height=7, 
       units = c("mm"),
       height = 90,
       width = 140,
       dpi=300)


ggsave("output/barplot_genomes_25_11_07.svg", combined_plot, 
       #width=15, 
       #height=7, 
       units = c("mm"),
       height = 90,
       width = 140,
       dpi=300)


#





############################### Old tree ################################

# Tree file
tree_drep <- read.tree("~/data/MFD/genome_phylogeny/methanotrophs/Methanotrophs_genome_mapping/with_gtdb/MSA.faa.treefile")

tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

tax<-readRDS("tax_MAGs_gtdb.rds")%>%filter(user_genome %in% tree_drep$tip.label)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'LR-MAG','GTDB'),
         type=if_else(grepl("LIB", user_genome), "SR-MAG", type))%>%
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", user_genome ))%>%
  mutate(Family_label=gsub("f__","", Family))


tr <- phytools::midpoint.root(tree_drep)
# plot rooted tree
pr <- ggtree(tr) 
pr

### Now, I need to remove the ones that are duplicates 

tax_unique <- tax %>%
  filter(case_when(
    type == "GTDB" ~ TRUE,
    type %in% c("LR-MAG", "SR-MAG") & Species == "s__" ~ TRUE,
    TRUE ~ FALSE
  ))%>%
  mutate(hit=if_else(type=="GTDB", NA, "yes"))


### Then, I will remove those from the tree
unique_tips <- tax_unique %>%
  pull(user_genome)

# remove outgroup,I use this instead
Unique_tree <- ape::keep.tip(tr, tip = unique_tips)
Unique_plot<-ggtree(Unique_tree,size=0.35, layout="fan") %<+% tax + 
  geom_treescale(x=0.3, y=Unique_tree$Nnode/2, label="Tree scale", fontsize =3.5, width = 0.3, offset.label = 1.5, offset=-4) 



#rotate_tree(Unique_plot, 180)


#### Now I need to make a color scale for the Family ######
unique(tax_unique$Order)
Order_colors = c('#F5F5F5','#bf4e4e',"#ACABFC90",'#9c760560',"#AFEEEE",'#3559e6', '#b2df8a','#e31a1c50','#418cb2','gold4','#34978F90', 'orange3' )
Order_colors <- rev(Order_colors)


tip_data <- tax_unique %>% select(user_genome, hit) 
p_save <- Unique_plot %<+% tip_data 

# combine
p_savephyl <- p_save +
  geom_tippoint(mapping=aes(color=hit),
                size=2, stroke=0)+
  scale_color_manual(values = '#dd0305',
                     na.value = NA, guide="none")


p_savephyl_fam <- p_savephyl + new_scale_fill() + geom_fruit(
  geom = geom_tile,
  mapping = aes(y=tip_label, group=label,
                fill=Family_label),
  #alpha = 0.8,
  width=0.25,
  offset= -0.02 #before -0.03
) +
  scale_fill_manual('Family',
                    values = Order_colors)



p_savephyl_fam

#### If I want to add annotations #####
# 
# filt<-tax%>%filter(Genus%in%c("g__Methylobacter","g__USCg-Taylor", "g__Methylobacter_A","g__Methylobacter_C"))
# p2<-ggtree(Unique_tree) %<+% filt+
#   geom_text(aes(x=branch, label=node), size=3,
#             vjust=-.3, color="firebrick")
# 
# p2+new_scale_fill() + geom_fruit(
#   geom = geom_tile,
#   mapping = aes(y=tip_label, group=label,
#                 fill=Genus),
#   alpha = 0.6,
#   width=0.2,
#   offset= -0.2 #before -0.03
# )

dat <- data.frame(
  node = c(951, 881,1001,690,578,870),
  name = c("Methylocella", "Methylocystis","Rhodomicrobium","Methylobacter_A", "Methylobacter_C","USCg-Taylor")
)

dat
# p_savephyl_zooglan+
#   geom_cladelab(data=dat, mapping = aes(node=node,label=name),
#                 fontsize = 3.3,
#                 align = TRUE,
#                 show.legend = FALSE,
#                 angle = c(30,0,51,29,-17,-53),
#                 hjust=c(0.5,0.5,0.5,0.5,0.5,0.5),
#                 vjust=c(0.5,0.5,0.5,0.5,0.5,0.5),
#                 horizontal =F,offset= 0.12, offset.text=0.08
#   ) 

p_savephyl_fam_annot<-p_savephyl_fam+
  geom_cladelab(data=dat, mapping = aes(node=node,label=name),
                fontsize = 3.3,
                align = TRUE,
                show.legend = FALSE,
                angle = "auto",
                #hjust=0.5,
                hjust = "center",
                barsize = 0.3,
                #  vjust=0.5,
                horizontal =F,offset= 0.12, offset.text=0.08
  )

p_savephyl_fam_annot
ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_07_29.svg', height=8, width = 12)
ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_07_29.png', height=8, width = 12)



# Count the number of added genomes per genus
genus_summary <- tax_unique %>%
  group_by(Genus) %>%
  summarize(
    added_genomes = sum(hit == "yes", na.rm = TRUE),   # Number of genomes being added
    gtdb_genomes = sum(type == "GTDB", na.rm = TRUE),  # Number of GTDB genomes per genus
    total_genomes = n()                                # Total number of genomes per genus
  ) %>%
  mutate(
    percent_increase = round((added_genomes / (gtdb_genomes + added_genomes)) * 100, 0) # Percent increase
  ) %>%
  ungroup()

# Merge the genus_summary back to the original tax_unique dataframe
tax_final <- tax_unique %>%
  left_join(genus_summary, by = "Genus")%>%
  select(!c("user_genome", "classification", "Species", "tree_label", "type", "hit"))%>%
  distinct()%>%
  filter(percent_increase>10)

# View the final dataframe
print(tax_final)







##########################################
########### Naming #######################


linkage<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/MFD_LR_SR_r220_sp_cluster.tsv", delim="\t", col_names = c("user_genome","classification"))%>%
  mutate(user_genome=gsub(".fa","", user_genome))%>%
  separate(classification, into = c("Domain", "Phylum", "Class","Order", "Family", "Genus", "Species_cluster"), sep = ';', remove=T)%>%
  select(user_genome, Species_cluster)

#quality<-vroom("/projects/microflora_danica/deep_metagenomes/mags_v4/MFD_mags_V4.tsv")%>%
quality<-vroom("/projects/microflora_danica/MFD-LR/mags_v4/MFD_mags_V4.tsv")%>%
  select(bin, MAG_status, gtdb_classification, Completeness_CheckM1, Contamination_CheckM1, Completeness_CheckM2, Contamination_CheckM2)%>%rename(user_genome=bin)

quality2<-vroom("/projects/microflora_danica/MFD-SR/results/mags_shallow_all.tsv")%>%
  select(bin, MAG_status, gtdb_classification, Completeness_CheckM1, Contamination_CheckM1, Completeness_CheckM2, Contamination_CheckM2)%>%rename(user_genome=bin)
qual<-rbind(quality, quality2)

tax_og<-readRDS("tax_MAGs_gtdb.rds")%>%
  left_join(qual)%>%
 # filter(str_detect(classification, str_c(tax_filter, collapse = "|")))%>%
  left_join(linkage)%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'LR-MAG','GTDB'),
         type=if_else(grepl("LIB", user_genome), "SR-MAG", type))%>%
  mutate(novel = if_else(type %in% c("LR-MAG", "SR-MAG") & Species == "s__", "yes", "no"))

name<-readxl::read_excel("dataframes/naming_methanotrophs.xlsx")%>%
  rename(sp_label=Species)%>%
  rename(g_label=Genus)



name_add_TUSC<-tax_og%>%filter(str_detect(classification, str_c(name$g_label, collapse = "|")))%>%
  filter(grepl("o__JACPRU01", Order))


### None are HQ - I will not persure this


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

tax<-tax%>%
  mutate(label_3=gsub("Methylocella", "Methylocapsa", label_3))%>%
  mutate(label_3=gsub("Ca. Methyloaffinis", "Ca. M.", label_3))%>%
  mutate(label_3=gsub("Methylocapsa sp003162995", "M. sp003162995", label_3))
  


name_filt<-name%>%left_join(tax, by = c("sp_label" = "label_3"))%>%
  left_join(qual)%>%
  #mutate(MAG_id=tip.label)%>%
  mutate(Qual=MAG_status)%>%
  select(!c("gtdb_classification", "MAG_id", "MAG_status", "tip.label", "type", "Domain", "MAG_Flag", "genome_number", "label_2", "Phylum", "Class"))
  

### getting sp_clust for the new genomes ####
sp_clust<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/methanotrophs/Troms_enrich_oxi_deep_seq_25_10_10/drep_output/data_tables/Cdb.csv")%>%
  mutate(genome=gsub(".fa", "", genome))

tax_new<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/methanotrophs/Troms_enrich_oxi_deep_seq_25_10_10/mmlong2_bins.tsv")%>%
  filter(bin %in% sp_clust$genome)%>%
  select(bin, completeness_checkm1, completeness_checkm2, contamination_checkm1, contamination_checkm2, r_abund, bin_status, gtdb_tax)%>%
  separate(gtdb_tax, into = c("Domain", "Phylum", "Class","Order", "Family", "Genus", "Species"), sep = ';', remove=T)%>%
  left_join(sp_clust, by=c("bin"="genome"))


sp_shared<-name_filt%>%left_join(sp_clust, by=c("user_genome"="genome"))%>%
  filter(secondary_cluster %in% tax_new$secondary_cluster)

sp_new<-tax_new%>%
  filter(!secondary_cluster %in% sp_shared$secondary_cluster)%>%
  filter(!bin=="MFD04830-1.bin.c.10")%>%
  select(!c("Domain","Phylum", "Class", "threshold", "secondary_cluster", "cluster_method", "comparison_algorithm", "primary_cluster"))


by_cols <- c("Order", "Family", "Genus", "Species",
             "Completeness_CheckM1" = "completeness_checkm1",
             "Contamination_CheckM1" = "contamination_checkm1",
             "Completeness_CheckM2" = "completeness_checkm2",
             "Contamination_CheckM2" = "contamination_checkm2",
             "Qual"="bin_status",
              "user_genome"="bin")

combined_df <- full_join(name_filt, sp_new, by = by_cols)%>%
  mutate(g_label=if_else(is.na(g_label), paste0(gsub("g__", "", Genus)), g_label))%>%filter(!g_label=="JACPRU01")


writexl::write_xlsx(combined_df, path = "dataframes/naming_methanotrophs_filt_R.xlsx")


                            