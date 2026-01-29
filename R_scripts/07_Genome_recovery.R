
#!/usr/bin/env Rscript
## ============================================================================
## 07_Genome_recovery.R
## ============================================================================
setwd("path/to/your/repo/MFD_methanotrophs_DK/")

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
library(ggnewscale)
library(ggtreeExtra)
library(ggbreak)
library(ggpattern)


setwd("path/to/your/repo/MFD_methanotrophs_DK/")


# Tree file
tree_drep <- read.tree("data/phylogenomic_new_spe_recovery.treefile")

tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

tax<-readRDS("data/tax_MAGs_gtdb.rds")%>%filter(user_genome %in% tree_drep$tip.label)%>%
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
ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_10_29.svg', height=8, width = 12)
ggsave(p_savephyl_fam_annot, filename = './output/Genome_recovery_24_10_29.png', height=8, width = 12)



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

tax<-readRDS("data/tax_MAGs_gtdb.rds")%>%
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
        fill = "white",         
        pattern_fill = "white",   
        colour = "black"          
      )
    ),
    fill = "none",              
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

