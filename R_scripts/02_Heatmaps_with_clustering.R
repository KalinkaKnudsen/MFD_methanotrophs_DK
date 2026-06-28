#!/usr/bin/env Rscript

# 02_Heatmaps_with_clustering.R
# Purpose: generate heatmaps of gene abundance (RPKM) across samples,
#          cluster samples within habitat groups, and export publication-ready
#          figures and legends used in the ISME presentation.
# Usage: set working directory to repository root (or edit the `setwd()` below)
# Dependencies: ggplot2, cowplot, vroom, tidyverse, readxl, patchwork, vegan, ggrastr, ggpubr, rstatix

library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
library(patchwork)


## NOTE: adjust this path to your local clone if needed
setwd("path/to/your/repo/MFD_methanotrophs_DK/")




## Load color palettes and pre-processed OTU table (long format)
### OBS, color palette is located in the accessory files, and must be unzipped for use
p_combine<-readRDS("data/palette_mfd_hab2_ISME.rds") 
OTU_filtered_long<-readRDS("output/OTU_filtered_long_25_09_01.rds") ## Generated in 01_scale_OTU.R, contains RPKM values and curated taxonomy for plotting



## ----- Taxonomy curation and ordering -----
# The following block standardises and simplifies long `Tax` strings into
# `Tax_curated` for clearer plotting labels and grouping.

#Ordering the tax:

tax<-OTU_filtered_long%>%
  select(Tax)%>%
  distinct()

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



# Merge curated taxonomy back onto the main OTU dataframe
## This keeps original rows and attaches the simplified `Tax_curated` used later in plots
OTU_filtered_long<-merge(OTU_filtered_long, tax_curated, by= "Tax")

#Cleaning up the habitat names a bit
## Clean and normalise habitat labels, create composite labels used for faceting
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


#### Simple prevalence & summary statistics #####
# Compute presence/absence and simple summary stats (median, mean) for Methylococcaceae
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



## Create a wide OTU matrix (rows = SeqId, cols = Tax) and drop rows with all zeros

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


groups2 <- unique(OTU_filtered_long$mfd_hab2)

for (i in seq_along(groups2)) {
  
  gr_name <- paste0("gr_", i)
  
  current_group <- OTU_filtered_long %>%
    filter(mfd_hab2 %in% groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = SeqId) %>%
    as.data.frame(.)
  
  nonzero <- current_group %>%
    filter(rowSums(select(., -SeqId)) > 0)
  
  zero_ids <- setdiff(current_group$SeqId, nonzero$SeqId)
  
  # --- CASE 1: 2 or more non-zero SeqIds ---
  if (nrow(nonzero) >= 2) {
    
    OTU_sums <- nonzero[,-1]
    row.names(OTU_sums) <- nonzero$SeqId
    
    bray_curtis_dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
    hclust_ward <- hclust(bray_curtis_dist, method = "ward.D2")
    ward_order <- hclust_ward$labels[order.dendrogram(as.dendrogram(hclust_ward))]
    
    levels_gr_i <- c(ward_order, zero_ids)
    
    # --- CASE 2: exactly 1 non-zero SeqId ---
  } else if (nrow(nonzero) == 1) {
    
    levels_gr_i <- c(nonzero$SeqId, zero_ids)
    
    # --- CASE 3: zero non-zero SeqIds ---
  } else {
    
    levels_gr_i <- current_group$SeqId
  }
  
  assign(paste0("levels_", gr_name), levels_gr_i)
}

levels_hab2 <- unique(unlist(mget(paste0("levels_gr_", seq_along(groups2)))))





### Filtering to selected genes for display

Tax_filter<-c("TUSC", "Likely_mmoX", "AVCC01_Methylotenera_clade", "Methylomonadaceae; Methyloprofundus_WTBX01_clade",  #"USCg_Methyloglobulus",
              "Methylococcales_unknown", "o_Rhizobiales_pmoA", "Methylococcaceae; Methylogaea", "Methylomagnum", 
              "Rhizobiales; Beijerinckiaceae; Methyloferula_Methylovirgula", "Methylomonadaceae; SXIZ01", "Methylomonadaceae; Methyloprofundus", "pxmA_umbrella", "Methylohalobius", "Methylococcales_mmoX_umbrella", "Beijerinckiaceae_pmoA1", "Methylomonadaceae; Methylomarinum")

hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater")

OTU_filtered_long<-OTU_filtered_long %>%
  filter(!type=="Put. pmoA/mmoX")%>%
  filter(!Tax_curated %in% Tax_filter)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>%
  mutate(Tax_curated=gsub("_mmoX_umbrella","", Tax_curated))%>%
  mutate(Tax_curated=gsub("_pmoA1", "", Tax_curated))




################# The re-do this #####################


# ---------------------------------------------------------
# 0. Build OTU_wide
# ---------------------------------------------------------
OTU_wide <- OTU_filtered_long %>%
  pivot_wider(names_from = Tax, values_from = RPKM, id_cols = SeqId) %>%
  as.data.frame() %>%
  filter(rowSums(select(., -SeqId)) > 0)

OTU_sums <- OTU_wide[,-1]
row.names(OTU_sums) <- OTU_wide$SeqId

# ---------------------------------------------------------
# 1. Extract ORIGINAL habitat order (your correct ordering)
# ---------------------------------------------------------
hab2_levels <- OTU_filtered_long %>%
  arrange(hab1_label) %>%
  distinct(mfd_hab2) %>%
  pull(mfd_hab2)

# ---------------------------------------------------------
# 2. Cluster SeqIds WITHIN each habitat (including zero-SeqIds)
# ---------------------------------------------------------
groups2 <- hab2_levels
levels_list <- vector("list", length(groups2))

for (i in seq_along(groups2)) {
  
  current_group <- OTU_filtered_long %>%
    filter(mfd_hab2 == groups2[i]) %>%
    pivot_wider(names_from = Tax, values_from = RPKM, id_cols = SeqId) %>%
    as.data.frame()
  
  nonzero <- current_group %>% filter(rowSums(select(., -SeqId)) > 0)
  zero_ids <- setdiff(current_group$SeqId, nonzero$SeqId)
  
  if (nrow(nonzero) >= 2) {
    OTU_sums <- nonzero[,-1]
    row.names(OTU_sums) <- nonzero$SeqId
    
    dist <- vegan::vegdist(vegan::decostand(OTU_sums, method = "hellinger"))
    hc <- hclust(dist, method = "ward.D2")
    ord <- hc$labels[order.dendrogram(as.dendrogram(hc))]
    
    levels_list[[i]] <- c(ord, zero_ids)
    
  } else if (nrow(nonzero) == 1) {
    levels_list[[i]] <- c(nonzero$SeqId, zero_ids)
    
  } else {
    levels_list[[i]] <- current_group$SeqId
  }
}

# ---------------------------------------------------------
# 3. Final SeqId ordering (preserves habitat block order)
# ---------------------------------------------------------
levels_hab2 <- unique(unlist(levels_list, use.names = FALSE))

# ---------------------------------------------------------
# 4. Apply ordering + clean type labels
# ---------------------------------------------------------
OTU_filtered_long <- OTU_filtered_long %>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow")) %>%
  mutate(
    Tax_curated = gsub("Methylocella", "Methylocapsa", Tax_curated),
    SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE),
    mfd_hab2 = factor(mfd_hab2, levels = hab2_levels),
    type = gsub("mmoX", "*mmoX*", type),
    type = gsub("pxmA", "*pxmA*", type),
    type = gsub("pmoA", "*pmoA*", type),
    type = gsub("\\*pmoA\\*2", "*pmoA2*", type)
  )

levels(OTU_filtered_long$type) <- c("*mmoX*", "*pmoA*", "*pmoA2*", "*pxmA*")

# ---------------------------------------------------------
# 5. HEATMAP
# ---------------------------------------------------------
heat <- OTU_filtered_long %>%
  ggplot(aes(y = Tax_curated, x = SeqId, fill = RPKM)) +
  ggrastr::geom_tile_rast() +
  scale_fill_gradientn(
    name = "RPKM",
    colors = c("white", "#82A3CD", "darkorange", "darkred", "black"),
    trans = "sqrt",
    limits = c(0, 5),
    na.value = "black"
  ) +
  facet_grid(type ~ hab1_label, scales = "free", space = "free", switch = "y") +
  theme(
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 5, margin = margin(r = -2)),
    legend.position = "none",
    strip.clip = "off",
    strip.text.x = element_text(size = 5, face = "bold"),
    strip.text.y = element_text(size = 5, face = "bold", margin = margin(l = 0.8, r = 0.8)),
    strip.text.y.left = ggtext::element_markdown(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2),
    strip.background = element_blank(),
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    panel.spacing = unit(0.08, "lines")
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# ---------------------------------------------------------
# 6. BAR PLOT (habitat color bar)
# ---------------------------------------------------------
bar_plot <- OTU_filtered_long %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ hab1_label, scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.08, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.2)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0, -1, 0, 0), "cm"))



combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5)) + theme(panel.background = element_blank(), plot.background = element_blank())
combined_plot



ggsave("output/gene_abundance_red.png",
       combined_plot,
       units = c("mm"),
       height = 95,
       width = 170,
       dpi=300)

ggsave("output/gene_abundance_red.svg",
       combined_plot,
       units = c("mm"),
       height = 95,
       width = 170,
       dpi=300)

########################################



saveRDS(levels_hab2, "26_12_06_seqid_order_hab2.rds")
saveRDS(hab2_levels, "26_12_06_habitat_order_hab2.rds")


#########  Extract per-habitat legends ######
## Creates separate small bar plots for each habitat and extracts their legends
## These are later combined into a single legend image for use in figures.

#######################################


ISME_colours<-readRDS("./palette_mfd_hab2_ISME.rds")

habitats <- c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Greenspaces", "Freshwater")

# Loop through each habitat
for (habitat in habitats) {
  # Filter the data for the current habitat
  plot_data <- OTU_filtered_long %>%
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

ggsave("output/ISME_gene_legend.svg", height=3, width=30)




###############################################################################
###### Investigate abundance of specific groups (TUSC, Binatales, putative groups) #####
###############################################################################

### Filtering to selected genes for display
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
  
  
}
######



levels_hab2 <- unlist(mget(paste0("levels_gr_", 1:length(groups2))))


hab2_sort<-OTU_filtered_long_put%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

hab2_sort_order<-hab2_sort%>%
  select(mfd_hab2)%>%distinct()%>%
  pull(mfd_hab2)
#saveRDS(hab2_sort_order, "hab2_sort_order.rds")

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






ggsave("output/Gene_abundance_put_only_25_12_11.png",
       combined_plot,
       units = c("mm"),
       height = 23,
       width = 175,
       dpi=300)

ggsave("output/Gene_abundance_put_only_25_2_11.svg",
       combined_plot,
       units = c("mm"),
       height = 23,
       width = 175,
       dpi=300)






###############################################################################
########## Investigating sub-habitats that could be interesting ####
###############################################################################

#####  Habitats for deep-dive
# Dunes: Rhodomicrobium + Methylocystis
# Grasslands: Rhodomicrobium + Gammas + cystis
# Temperate heath: as above + pmoA2
# Agriculture and USCa

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

ggsave("output/dune_grass_heath_hab3_25_11_14.svg",
       combined_plot,
       units = c("mm"),
       height = 110,
       width = 180,
       dpi=300)


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
scale_fill_gradientn(
    name = "RPKM",  
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  
    limits = c(0.00, 5), 
    na.value = "black"  
  ) +
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
        strip.text.x = element_text(size = 5, face="bold", margin=margin(b=0.8)),  
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

ggsave("output/fields_hab2_25_11_14.svg",
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

ggsave("output/fields_hab2.png",
       combined_plot,
       height = 10,
       width = 15,
       dpi=400)





############## Constructing a custom palette for the  heatmap ##############

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

