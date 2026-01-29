#!/usr/bin/env Rscript

library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(patchwork)


setwd("path/to/your/repo/MFD_methanotrophs_DK/")




p_combine<-readRDS("../palette_mfd_hab2_ISME.rds") 
OTU_filtered_long<-readRDS("OTU_filtered_long_25_09_01.rds")


# summary<-OTU_filtered_long%>%
#   filter(mfd_hab2=="Asterids")%>%
#   select(SeqId, mfd_hab3)%>%
#   distinct()%>%
#   count(mfd_hab3)



#Ordering the tax:

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
         Tax_curated=gsub("Methylosinus_Methylocystis_pmoA2; ", "", Tax_curated),
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


#### Simple prevance stuff #####
Methylococcaceae_presence <- OTU_filtered_long %>%
  filter(Tax_short == "Methylococcaceae") %>%
  group_by(SeqId, mfd_hab1) %>%
  summarise(
    present = any(RPKM > 0),
    RPKM = RPKM,
    .groups = "drop"
  )

prevalence_by_habitat <- Methylococcaceae_presence %>%
  group_by(mfd_hab1) %>%
  summarise(
    n_samples = n(),
    n_positive = sum(present),
    prevalence = n_positive / n_samples,
    median_RPKM = median(RPKM),
    mean_RPKM=mean(RPKM)
  )



OTU_wide<-OTU_filtered_long%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
#saveRDS(file="OTU_wide.rds", OTU_wide) ## To use in PCA
groups<-unique(OTU_filtered_long$complex_long)


#length(groups2)
#Trying to plot with the clustering of MFD_hab02
#Cluster individually

groups2<-unique(OTU_filtered_long$mfd_hab2)


for (i in seq_along(groups2)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
  
  
  # OTU_filtered_long<-OTU_filtered_long%>%
  #   mutate(mfd_hab2=factor(mfd_hab2, levels = groups2[i]), ordered = TRUE)
  
}

levels_hab2 <- c(
  levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22, levels_gr_23, levels_gr_24, levels_gr_25,
  levels_gr_26, levels_gr_27, levels_gr_28, levels_gr_29, levels_gr_30,
  levels_gr_31, levels_gr_32, levels_gr_33, levels_gr_34, levels_gr_35,
  levels_gr_36, levels_gr_37, levels_gr_38, levels_gr_39, levels_gr_40,
  levels_gr_41, levels_gr_42, levels_gr_43, levels_gr_44, levels_gr_45,
  levels_gr_46, levels_gr_47, levels_gr_48, levels_gr_49, levels_gr_50,
  levels_gr_51, levels_gr_52, levels_gr_53, levels_gr_54, levels_gr_55,
  levels_gr_56, levels_gr_57)




#### Making a massive plot just so I can see which groups to exclude #####



heat <- OTU_filtered_long %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # filter(Tax=="Root; o_Rhizobiales_pmoA")%>%
  #mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
  # filter(complex_long %in% groups[i]) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  geom_tile(width=2.0, linewidth=0.0) +
  scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
  labs(x = "", y = "") +
  facet_grid(type ~ mfd_areatype, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.5),
        legend.position = "right",
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "white"),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.05, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))


heat



ggsave("heatmap_all.png",
       heat,
       height = 10,
       width = 35)


######################################################################################################################
###################### Now I will zoom in to the same habitats as for the genomes - for ISME presentation ############
######################################################################################################################


Tax_filter<-c("TUSC", "Likely_mmoX", "AVCC01_Methylotenera_clade", "Methylomonadaceae; Methyloprofundus_WTBX01_clade",  #"USCg_Methyloglobulus",
              "Methylococcales_unknown", "o_Rhizobiales_pmoA", "Methylococcaceae; Methylogaea", "Methylomagnum", 
              "Rhizobiales; Beijerinckiaceae; Methyloferula_Methylovirgula", "Methylomonadaceae; SXIZ01", "Methylomonadaceae; Methyloprofundus", "pxmA_umbrella", "Methylohalobius", "Methylococcales_mmoX_umbrella", "Beijerinckiaceae_pmoA1", "Methylomonadaceae; Methylomarinum")

hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater")

OTU_filtered_long_ISME<-OTU_filtered_long %>%
  filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))




################# The re-do this #####################

OTU_wide<-OTU_filtered_long_ISME%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
#saveRDS(file="OTU_wide_SIME.rds", OTU_wide) ## To use in PCA
OTU_wide<-readRDS("OTU_wide_SIME.rds")
groups<-unique(OTU_filtered_long_ISME$complex_long)


#length(groups2)
#Trying to plot with the clustering of MFD_hab02
#Cluster individually

groups2<-unique(OTU_filtered_long_ISME$mfd_hab2)


for (i in seq_along(groups2)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long_ISME %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
  
  
  # OTU_filtered_long<-OTU_filtered_long%>%
  #   mutate(mfd_hab2=factor(mfd_hab2, levels = groups2[i]), ordered = TRUE)
  
}
######


levels_hab2 <- c(
  levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22, levels_gr_23, levels_gr_24, levels_gr_25,
  levels_gr_26, levels_gr_27, levels_gr_28, levels_gr_29)


hab2_sort<-OTU_filtered_long_ISME%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)
#saveRDS(hab2_sort_order, "hab2_sort_order_ISME.rds")

palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

OTU_filtered_long_ISME <- OTU_filtered_long_ISME %>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))

OTU_filtered_long_ISME<-OTU_filtered_long_ISME%>%
  mutate(type=gsub("mmoX", "*mmoX*", type),
         type=gsub("pxmA", "*pxmA*", type),
         type=gsub("pmoA", "*pmoA*", type),
         type = gsub("\\*pmoA\\*2", "*pmoA2*", type))          

levels(OTU_filtered_long_ISME$type)<-c("*mmoX*",
                                  "*pmoA*",
                                  "*pmoA2*", "*pxmA*")


heat <- OTU_filtered_long_ISME %>%
  mutate(Tax_curated=gsub("Methylocella", "Methylocapsa", Tax_curated))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  # filter(!type=="Put. pmoA/mmoX")%>%
  # filter(!Tax_curated %in% Tax_filter)%>%
  # filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  ggrastr::geom_tile_rast() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5, margin = margin(r = -2)),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.position = "none",
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 5, face="bold", margin=margin(l=0.8, r=0.8)),
        strip.text.y.left = ggtext::element_markdown(),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#heat

bar_plot <- OTU_filtered_long_ISME %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.08, "lines", data = NULL),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,-1,0,0), "cm"))



# Step 3: Combine the heatmap plot and the bar plot

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot



ggsave("ISME_gene_abundance_red.png",
       combined_plot,
       units = c("mm"),
       height = 95,
       width = 170,
       dpi=300)

ggsave("ISME_gene_abundance_red.svg",
       combined_plot,
       units = c("mm"),
       height = 95,
       width = 170,
       dpi=300)

########################################


#########  Getting the legend ######


#######################################


ISME_colours<-readRDS("./palette_mfd_hab2_ISME.rds")

habitats <- c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Greenspaces", "Freshwater")

# Loop through each habitat
for (habitat in habitats) {
  # Filter the data for the current habitat
  plot_data <- OTU_filtered_long_ISME %>%
    filter(mfd_hab1 == habitat) %>%
    mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
    mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
    filter(!mfd_hab2=="Mire")%>%
    filter(!mfd_sampletype=="Water")
  
  # Create the plot
  bar_plot <- ggplot(plot_data, aes(x = SeqId, fill = mfd_hab2)) +
    geom_tile(aes(y = 1)) +
    facet_nested(. ~ hab1_label, scales = "free", space = "free") +
    scale_fill_manual(values = ISME_colours) +
    labs(fill = habitat) +
    theme(
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom",
      panel.spacing = unit(0.1, "lines"),
      legend.text = element_text(size=11),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0))+
    guides(fill = guide_legend(nrow = 4))
  
  legend <- ggpubr::get_legend(bar_plot)
  leg<-ggpubr::as_ggplot(legend)
  assign(paste0("leg_", gsub(", ", "_", gsub(" ", "_", habitat))), ggpubr::as_ggplot(legend))
  
  # Save the plot as SVG
  #ggsave(filename = paste0("./output/legend_habitat/bar_plot_", gsub(", ", "_", gsub(" ", "_", habitat)), ".svg"), plot = leg)
}

leg_comb<-`leg_Bogs,_mires_and_fens` + leg_Dunes + leg_Forests + leg_Freshwater + leg_Grassland_formations + leg_Greenspaces + leg_Temperate_heath_and_scrub + plot_layout(nrow=1)

ggsave("./legend_habitat/ISME_gene_legend.svg", height=3, width=30)








###############################################################################
###### Checking the abundance of the TUSC and Binatales groups #####
###############################################################################









Tax_filter<-c("Likely_mmoX", "AVCC01_Methylotenera_clade", "Methylomonadaceae; Methyloprofundus_WTBX01_clade",  #"USCg_Methyloglobulus",
              "Methylococcales_unknown", "o_Rhizobiales_pmoA", "Methylococcaceae; Methylogaea", "Methylomagnum", 
              "Rhizobiales; Beijerinckiaceae; Methyloferula_Methylovirgula", "Methylomonadaceae; SXIZ01", "Methylomonadaceae; Methyloprofundus", "pxmA_umbrella", "Methylohalobius", "Methylococcales_mmoX_umbrella", "Beijerinckiaceae_pmoA1", "Methylomonadaceae; Methylomarinum")

hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater", "Fields")

OTU_filtered_long_put<-OTU_filtered_long %>%
  #filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater", "Fields"))%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))




################# The re-do this #####################

OTU_wide<-OTU_filtered_long_put%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
#saveRDS(file="OTU_wide_SIME.rds", OTU_wide) ## To use in PCA
#OTU_wide<-readRDS("OTU_wide_SIME.rds")
groups<-unique(OTU_filtered_long_put$complex_long)


#length(groups2)
#Trying to plot with the clustering of MFD_hab02
#Cluster individually

groups2<-unique(OTU_filtered_long_put$mfd_hab2)


for (i in seq_along(groups2)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long_put %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
  
  
  # OTU_filtered_long<-OTU_filtered_long%>%
  #   mutate(mfd_hab2=factor(mfd_hab2, levels = groups2[i]), ordered = TRUE)
  
}
######



levels_hab2 <- unlist(mget(paste0("levels_gr_", 1:length(groups2))))


hab2_sort<-OTU_filtered_long_put%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)
#saveRDS(hab2_sort_order, "hab2_sort_order_ISME.rds")

palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

OTU_filtered_long_put <- OTU_filtered_long_put %>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))


heat <- OTU_filtered_long_put %>%
  mutate(type=if_else(Tax_curated=="TUSC", paste0("Put."), type))%>%
  mutate(type=if_else(type=="Put. pmoA/mmoX", paste0("Put."), type))%>%
  filter(type=="Put.")%>%
  filter(!grepl("Burkhold", Tax_curated))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  ggrastr::geom_tile_rast() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5, margin = margin(r = -1)),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.position = "none",
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin=margin(b=0.8)),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 5, face="bold", margin=margin(l=1.6, r=0.8)),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0.04,0), "cm"))




colors<-unlist(readRDS("../../../nitrifiers/R_scripts/palette_mfd_hab2.rds"))
palette_mfd_hab2<-unlist(readRDS("./palette_mfd_hab2_ISME.rds"))

cols_combi<-c(colors, palette_mfd_hab2)

bar_plot <- OTU_filtered_long_put %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = cols_combi) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    panel.spacing = unit(0.08, "lines", data = NULL),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,-1,0,0), "cm"))



# Step 3: Combine the heatmap plot and the bar plot


combined_plot <- heat / bar_plot + patchwork::plot_layout(heights = c(5, 0.8)) + theme(panel.background = element_blank(), plot.background = element_blank())






ggsave("Gene_abundance_put_only_25_12_11.png",
       combined_plot,
       units = c("mm"),
       height = 23,
       width = 175,
       dpi=300)

ggsave("Gene_abundance_put_only_25_2_11.svg",
       combined_plot,
       units = c("mm"),
       height = 23,
       width = 175,
       dpi=300)






###############################################################################
######################### Investigating habitats that could be interesting ####
###############################################################################

#####  Habitats for deep-dive
# Dunes: Rhodomicrobium + Methylocystis
# Grasslands: Rhodomicrobium + Gammas + cystis
# Temperate heath: as above + pmoA2
# Agriculture: some cereals with USCg and USCα hotspots. Could this be location/reflected in other things




Tax_filter<-c("TUSC", "Likely_mmoX", "AVCC01_Methylotenera_clade", "Methylomonadaceae; Methyloprofundus_WTBX01_clade", #"USCg_Methyloglobulus",
               "Methylococcales_unknown", "o_Rhizobiales_pmoA", "Methylococcaceae; Methylogaea", "Methylomagnum", 
              "Rhizobiales; Beijerinckiaceae; Methyloferula_Methylovirgula", "Methylomonadaceae; SXIZ01", "Methylomonadaceae; Methyloprofundus", "pxmA_umbrella", "Methylohalobius", "Methylococcales_mmoX_umbrella", "Beijerinckiaceae_pmoA1", "Methylomonadaceae; Methylomarinum")

hab_filter<-c("Grassland formations", "Temperate heath and scrub", "Dunes")

OTU_filtered_long_hab3<-OTU_filtered_long %>%
  filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_hab1 %in% hab_filter)%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))%>%
  filter(!is.na(mfd_hab3))
 # mutate(mfd_hab3=if_else(is.na(mfd_hab3), paste0(mfd_hab2, " - NA"), mfd_hab3))%>%
#  filter(!mfd_hab3 %in% c("Inland dunes - NA", "Sea dunes - NA", "Sclerophyllous scrub - NA"))



## clustering at level hab3
OTU_wide<-OTU_filtered_long_hab3%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId
groups<-unique(OTU_filtered_long_hab3$complex_long)


groups3<-unique(OTU_filtered_long_hab3$mfd_hab3)


for (i in seq_along(groups3)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long_hab3 %>%
    filter(mfd_hab3 %in% groups3[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
######


levels_hab3 <- c(
  levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11, levels_gr_12, levels_gr_13, levels_gr_14, levels_gr_15,
  levels_gr_16, levels_gr_17, levels_gr_18, levels_gr_19, levels_gr_20,
  levels_gr_21, levels_gr_22)


hab2_sort<-OTU_filtered_long_hab3%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)



OTU_filtered_long_hab3<-OTU_filtered_long_hab3%>%
  mutate(type=gsub("mmoX", "*mmoX*", type),
         type=gsub("pxmA", "*pxmA*", type),
         type=gsub("pmoA", "*pmoA*", type),
         type = gsub("\\*pmoA\\*2", "*pmoA2*", type))          

levels(OTU_filtered_long_hab3$type)<-c("*mmoX*",
                                       "*pmoA*",
                                       "*pmoA2*", "*pxmA*")



heat <- OTU_filtered_long_hab3 %>%
  mutate(Tax_curated=gsub("Methylocella", "Methylocapsa", Tax_curated))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  ggrastr::geom_tile_rast() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5, margin = margin(r = -1)),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.position = "none",
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin=margin(b=0.8)),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 5, face="bold", margin=margin(l=1.6, r=0.8)),
        strip.text.y.left = ggtext::element_markdown(),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#heat

#palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

bar_plot <- OTU_filtered_long_hab3 %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab3)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = c("#cdad00",  "orange","#000000", "#5a946d", "#32849f", "#a81818", "#1f6f87", "#0D0887", "#663300","#e6c200", "#7301A8",
                               "#9855d4", "#99cc00", "#e066ff", "#cd2626", "#ff4d4d", "#ff9999", "#33ccff", "#0099cc", 
                               "#005580",  "#8CB88F", "yellow2", "#ff6600", "#444444")) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
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



#library(patchwork)

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot

ggsave("./dune_grass_heath_hab3_25_11_14.svg",
       combined_plot,
       units = c("mm"),
       height = 110,
       width = 180,
       dpi=300)


#





############## Now to fields


OTU_filtered_long_cereal<-OTU_filtered_long %>%
  filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_hab2 =="Poales, Cereal")%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))%>%
  filter(!is.na(mfd_hab3))
# mutate(mfd_hab3=if_else(is.na(mfd_hab3), paste0(mfd_hab2, " - NA"), mfd_hab3))%>%
#  filter(!mfd_hab3 %in% c("Inland dunes - NA", "Sea dunes - NA", "Sclerophyllous scrub - NA"))


not_include<-OTU_filtered_long_cereal%>%select(SeqId, mfd_hab3)%>%distinct()%>%group_by(mfd_hab3)%>%summarise(count=n())%>%filter(count<10)%>%pull(mfd_hab3)
OTU_filtered_long_cereal<-OTU_filtered_long_cereal%>%filter(!mfd_hab3 %in% not_include)

## clustering at level hab3
OTU_wide<-OTU_filtered_long_cereal%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

groups3<-unique(OTU_filtered_long_cereal$mfd_hab3)


for (i in seq_along(groups3)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long_cereal %>%
    filter(mfd_hab3 %in% groups3[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
######


levels_hab3 <- c(
  levels_gr_1, levels_gr_2, levels_gr_3, levels_gr_4, levels_gr_5,
  levels_gr_6, levels_gr_7, levels_gr_8, levels_gr_9, levels_gr_10,
  levels_gr_11)


hab2_sort<-OTU_filtered_long_cereal%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)

heat <- OTU_filtered_long_cereal %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=9),
        legend.position = "right",
        strip.clip = "off",
        strip.text.x = element_text(size = 9, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 9, angle = 0, face="bold"),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
heat

#palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")





bar_plot <- OTU_filtered_long_cereal %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab3, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab3)) +
  geom_tile(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = c("#cdad00",  "orange","#000000", "#5a946d", "#32849f", "#a81818", "#1f6f87", "#0D0887", "#663300","#e6c200", "#7301A8",
                               "#9855d4", "#99cc00", "#e066ff", "#cd2626", "#ff4d4d", "#ff9999", "#33ccff", "#0099cc", 
                               "#005580",  "#8CB88F", "yellow2", "#ff6600", "#444444")) + # Specify your colors
  
  
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8),
    legend.position = "bottom",
    panel.spacing = unit(0.08, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))




combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())

combined_plot

ggsave("./cereal_hab3.png",
       combined_plot,
       height = 10,
       width = 15,
       dpi=400)

#


########################### And all agriculture ###############


OTU_filtered_long_agri<-OTU_filtered_long %>%
  filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_areatype =="Agriculture")%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))%>%
  filter(!is.na(mfd_hab3))
# mutate(mfd_hab3=if_else(is.na(mfd_hab3), paste0(mfd_hab2, " - NA"), mfd_hab3))%>%
#  filter(!mfd_hab3 %in% c("Inland dunes - NA", "Sea dunes - NA", "Sclerophyllous scrub - NA"))


not_include<-OTU_filtered_long_agri%>%select(SeqId, mfd_hab2)%>%distinct()%>%group_by(mfd_hab2)%>%summarise(count=n())%>%filter(count<9)%>%pull(mfd_hab2)
OTU_filtered_long_agri<-OTU_filtered_long_agri%>%filter(!mfd_hab2 %in% not_include)

## clustering at level hab2
OTU_wide<-OTU_filtered_long_agri%>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId))%>% ##Change to tax curated once possible
  as.data.frame(.)%>%
  filter(rowSums(select(., -1)) > 0)
OTU_sums<-OTU_wide[,-1]
row.names(OTU_sums)<-OTU_wide$SeqId

groups2<-unique(OTU_filtered_long_agri$mfd_hab2)


for (i in seq_along(groups2)) {
  # Use paste0 to dynamically create variable names gr_1, gr_2, etc.
  gr_name <- paste0("gr_", i)
  
  # Use the filter condition to select the current group
  current_group <- OTU_filtered_long_agri %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = c(SeqId)) %>% ## Change to tax curated
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
######



levels_hab2 <- unlist(mget(paste0("levels_gr_", 1:length(groups2))))

hab2_sort<-OTU_filtered_long_agri%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)





OTU_filtered_long_agri<-OTU_filtered_long_agri%>%
  mutate(type=gsub("mmoX", "*mmoX*", type),
         type=gsub("pxmA", "*pxmA*", type),
         type=gsub("pmoA", "*pmoA*", type),
         type = gsub("\\*pmoA\\*2", "*pmoA2*", type))          

levels(OTU_filtered_long_agri$type)<-c("*mmoX*",
                                       "*pmoA*",
                                       "*pmoA2*", "*pxmA*")



heat <- OTU_filtered_long_agri %>%
  mutate(Tax_curated=gsub("Methylocella", "Methylocapsa", Tax_curated))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  ggrastr::geom_tile_rast() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free", switch = "y") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5, margin = margin(r = -1)),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.position = "none",
        strip.clip = "off",
        strip.text.x = element_text(size = 5, face="bold", margin=margin(b=0.8)),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 5, face="bold", margin=margin(l=1.6, r=0.8)),
        strip.text.y.left = ggtext::element_markdown(),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#heat

#palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")
colors<-unlist(readRDS("../../../nitrifiers/R_scripts/palette_mfd_hab2.rds"))



# ✅ Update Fabids to grey
colors["Fabids"] <- "grey30"  # grey

# ✅ Swap Fallow and Mixed crops
temp <- colors["Fallow"]
colors["Fallow"] <- colors["Mixed crops"]
colors["Mixed crops"] <- temp


bar_plot <- OTU_filtered_long_agri %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = colors)+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
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



#library(patchwork)

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())
#combined_plot

ggsave("./fields_hab2_25_11_14.svg",
       combined_plot,
       units = c("mm"),
       height = 95,
       width = 180,
       dpi=300)




heat <- OTU_filtered_long_agri %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  geom_tile() +
  #scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 5), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #labs(x = "", y = "") +
  facet_grid(type ~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=9),
        legend.position = "right",
        strip.clip = "off",
        strip.text.x = element_text(size = 9, face="bold"),  # Size for x-axis facet labels
        strip.text.y = element_text(size = 9, angle = 0, face="bold"),
        plot.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        strip.background = element_blank(),
        text = element_text(family = "Arial"),
        panel.spacing = unit(0.08, "lines", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
heat

#palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

p_combine<-unlist(readRDS("../../../nitrifiers/R_scripts/palette_mfd_hab2.rds"))




bar_plot <- OTU_filtered_long_agri %>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  geom_tile(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = p_combine)+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.key.size = unit(0.45, "cm"),
    legend.title = element_text(size=8),
    legend.text = element_text(size=8),
    legend.position = "bottom",
    panel.spacing = unit(0.08, "lines", data = NULL)
  )+   scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))




combined_plot <- heat / bar_plot + patchwork::plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())

combined_plot

ggsave("./fields_hab2.png",
       combined_plot,
       height = 10,
       width = 15,
       dpi=400)




###################################################################################################################
###################################################################################################################
                                        #"pmoA2 investigation"
###################################################################################################################
###################################################################################################################




pmoA2<-OTU_filtered_long %>%
  filter(Tax_curated==" Methylocystis_pmoA2")

not_include<-pmoA2%>%select(SeqId, mfd_hab3)%>%distinct()%>%group_by(mfd_hab3)%>%summarise(count=n())%>%filter(count<10)%>%pull(mfd_hab3)

pmoA2<-pmoA2%>%arrange(mfd_hab2)%>%arrange(mfd_hab1)

plot_test2 <- pmoA2%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Greenspaces", "Freshwater"))%>%
  filter(!mfd_hab3 %in% not_include)%>%
  filter(!is.na(mfd_hab3))%>%
  mutate(mfd_hab3 = factor(mfd_hab3, levels = unique(pmoA2$mfd_hab3), ordered = TRUE)) %>%
  ggplot(., aes(x = mfd_hab3, y = RPKM)) + 
  geom_jitter(aes(fill=mfd_hab1), size = 1.3, alpha = 0.8, color = "black", pch = 21,
              position = position_jitter(width = 0.2, height = 0.03) # Adjust jitter spread if needed
  ) +
  geom_violin(aes(fill=mfd_hab1),  # Add median quantile line
              color = "black",     # Median line color
              size = 0.6,
              alpha=0.4
              #fill="transparent"# Median line thickness
  ) +
  #  scale_fill_manual(values = p_combine)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.8, size = 0.6, color = "red3") +
  # stat_compare_means(label = "..p.format..", method = "wilcox.test",
  #                    ref.group = ".all.", method.args = list(alternative = "greater"))+ ###Adjust here for one sided "less" I think
  
  ggpubr::geom_pwc(method = "wilcox_test",
                   label = 'p.adj.signif',
                   ref.group = ".all.",
                   hide.ns = T,
                   #angle = 90, vjust = 0.75,
                   p.adjust.method = "fdr",
                   group.by = "x.var",
                   p.adjust.by = "panel",
                   method.args = list(alternative = "less"),
                   # method.args = list(alternative = "greater"),
                   remove.bracket = T)+
  #geom_hline(yintercept = mean(non_fields%>%pull(drep_abundance)), linetype = 2)+ # Add horizontal line at base mean
  scale_y_sqrt() +
  #  facet_grid(.~mfd_hab2, scales="free_x", space = "free")+
  theme_minimal()+
  theme(legend.position="bottom", panel.grid=element_blank(), 
        axis.text.x = element_text(angle=45, hjust = 1))+
  xlab("MFD ontology level 3")+
  ylab("RPKM of pmoA2 (Methylocystis)")

plot_test2

ggsave("./pmoA2_hab3.png",
       plot_test2,
       height = 10,
       width = 15,
       dpi=400)


plot_test2$layers


stats<-pmoA2%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Greenspaces", "Freshwater"))%>%
  filter(!mfd_hab3 %in% not_include)%>%
  filter(!is.na(mfd_hab3))%>%
  ungroup()%>%
  #group_by(mfd_hab3)%>%
  rstatix::wilcox_test(RPKM~mfd_hab3,alternative = "less", p.adjust.method = "fdr", ref.group = "all", detailed=T)%>%print(n = 24)




##################### Okay so I can try to do correlation within these enviornments? #############


lm_cladeB_P <- lm(`Nitrospira_clade_B`~`Palsa-1315*`, data = Nitrospira_BOG)
summary(lm_cladeB_P)

num_samples <- n_distinct(Nitrospira_BOG$SeqId)

habs<-stats%>%filter(!p.adj.signif=="ns")%>%pull(group2)

df<-OTU_filtered_long%>%filter(mfd_hab3 %in% habs)%>%
  filter(Tax_curated %in% c(" Methylocystis_pmoA2", "Beijerinckiaceae; Methylosinus_Methylocystis_pmoA1"))%>%
  select(SeqId, Tax_curated, RPKM, mfd_hab1, mfd_hab2, mfd_hab3)%>%
  pivot_wider(names_from = Tax_curated, values_from = RPKM)
  

plot <-df %>%
  ggplot(., aes(
    y = `Beijerinckiaceae; Methylosinus_Methylocystis_pmoA1`,
    x= ` Methylocystis_pmoA2`)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  #scale_fill_manual(values=p_combine, name = paste0("Bogs, Mires and Fens,\nn = ",num_samples))+
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  #theme_classic() +
  facet_wrap(.~mfd_hab3)+
  coord_cartesian(xlim = c(0, max(df$` Methylocystis_pmoA2`)), ylim = c(0,max(df$`Beijerinckiaceae; Methylosinus_Methylocystis_pmoA1`)))+
  theme(legend.position = c(0.88, 0.2),
        legend.background = element_rect(fill = "transparent"))+
  guides(fill = guide_legend(ncol = 2))

plot


ggsave("./pmoA2_hab3_linear_pmoA1.png",
       plot,
       height = 10,
       width = 15,
       dpi=400)







####################################################################################################################
####### Filtering samples for CoverM ######
CoverM_samples_1<-OTU_filtered_long%>%
  filter(!mfd_hab1=="Fields")%>%
  filter(!mfd_hab1=="Bogs, mires and fens")%>%
  filter(!mfd_hab1=="Saltwater")%>%
  filter(!mfd_hab1 %in% c("Wastewater", "Biogas", "Other"))%>%
  select(SeqId)%>%
  distinct()

readr::write_csv(CoverM_samples_1, "/home/bio.aau.dk/vj52ou/data/MFD/genome_mapping/CoverM_samples_no_fields.txt", col_names = FALSE)

CoverM_samples_fields<-OTU_filtered_long%>%
  filter(mfd_hab1=="Fields")%>%
  select(SeqId)%>%
  distinct()

readr::write_csv(CoverM_samples_fields, "/home/bio.aau.dk/vj52ou/data/MFD/genome_mapping/CoverM_samples_only_fields.txt", col_names = FALSE)

######################################## Investigating the subterranian part ##########################
colors <- c( "gold3","#33a08c","coral","#1a0060", "darkred", "#F8E622","darkgreen","#6a3d9a","#1f78b4", "firebrick3")

hab2_sort<-OTU_filtered_long%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
  arrange(SeqId)



sub_t <- OTU_filtered_long %>%
  filter(mfd_areatype=="Subterranean")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  # filter(Tax=="Root; o_Rhizobiales_pmoA")%>%
  #mutate(Tax_curated=factor(Tax_curated, levels = rev(tax2$Tax_order), ordered=TRUE))%>%
  # filter(complex_long %in% groups[i]) %>%
  ggplot(aes(y = Tax, x = SeqId, fill = RPKM)) +
  geom_tile(width=10.0, linewidth=0.0) +
  scale_fill_viridis_c(name = "RPKM", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  # scale_fill_gradient(low = "white", high = "darkred",  na.value = "darkred", name = "RPKM", trans = "sqrt", limits = c(0.00, 3), guide = "none") +
  labs(x = "", y = "") +
  facet_grid(type ~ mfd_areatype, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        #  axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.5),
        legend.position = "right",
        strip.text = element_text(size = 9.3, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "white"),
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent") )+
  scale_y_discrete(expand = c(0,0))

tile <- OTU_filtered_long %>%
  filter(mfd_areatype=="Subterranean")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = unique(hab2_sort$mfd_hab2))) %>%
  ggplot(., aes(x = SeqId, fill = mfd_hab2, y = 1)) +
  geom_tile()+
  theme(legend.position = "bottom",
        legend.text = element_text(size=9), legend.key.size = unit(0.5, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial")) +
  scale_fill_manual(values = colors)+
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(nrow = 4, title.theme = element_blank()))

assign(paste0("HM_plot_Subterranean"),
       sub_t %>% aplot::insert_bottom(tile, height = 0.02))


HM_plot_Subterranean










# 
# 
# 
# 
p_combine <- structure(list(`Poales, Cereal` = "#cdad00",
                            `Mixed crops` = "#32849f",
                            `Poales, grass` = "#7301A8FF",
                            Malvids = "#33ccff",
                            Superasterids = "#F8E622",
                            Fabids = "#ff9999",
                            Asparagales = "#CC4678FF",
                            `Running freshwater` = "#cdad00",
                            `Standing freshwater` = "#32849f",
                            `Deciduous trees` = "#0D0887FF" ,
                            `Forests - no MFDO2` = "#32849f",
                            `Non-native trees (exotic)` = "#cd2626",
                            `Coniferous forest` = "#cc33ff",
                            Beech = "#33ccff",
                            `Alluvial woodland` = "#F8E622",
                            Oak = "#ff9999",
                            `Bog woodland` = "#7301A8FF",
                            Willow = "#cdad00",
                            Spruce = "#000000",
                            `Mixed crops` = "#32849f",
                            Fallow = "#0D0887FF",
                            Asterids = "#cd2626",
                            `Running freshwater` = "#cdad00",
                            `Standing freshwater` = "#32849f",
                            `Rainwater basin, City` = "#cdad00",
                            `Rainwater basin, Dried` = "#32849f",
                            `Enclosed water` = "#7301A8FF",
                            `Rainwater basin, Roadside` = "#cd2626",
                            `Semi-natural dry grasslands` = "#32849f",
                            `Semi-natural tall-herb humid meadows` = "#cdad00",
                            `Natural grasslands` = "#7301A8FF",
                            Parks = "#cdad00",
                            Other = "#32849f",
                            `Calcareous fens` = "#cdad00",
                            `Sphagnum acid bogs` = "#32849f",
                            `Wet thicket` = "#7301A8FF",
                            `Activated sludge` = "#cdad00",
                            Influent = "#32849f"), class = "character",
                       row.names = c("Poales, Cereal",
                                     "Mixed crops",
                                     "Poales, grass",
                                     "Malvids",
                                     "Superasterids",
                                     "Fabids",
                                     "Asparagales",
                                     "Running freshwater",
                                     "Standing freshwater",
                                     "Deciduous trees",
                                     "Forests - no MFDO2",
                                     "Non-native trees (exotic)",
                                     "Coniferous forest",
                                     "Beech",
                                     "Alluvial woodland",
                                     "Oak",
                                     "Bog woodland",
                                     "Willow",
                                     "Spruce",
                                     "Mixed crops",
                                     "Fallow",
                                     "Asterids",
                                     "Running freshwater",
                                     "Standing freshwater",
                                     "Rainwater basin, City",
                                     "Rainwater basin, Dried",
                                     "Enclosed water",
                                     "Rainwater basin, Roadside",
                                     "Semi-natural dry grasslands",
                                     "Semi-natural tall-herb humid meadows",
                                     "Natural grasslands",
                                     "Parks",
                                     "Other",
                                     "Calcareous fens",
                                     "Sphagnum acid bogs",
                                     "Wet thicket",
                                     "Activated sludge",
                                     "Influent"))

