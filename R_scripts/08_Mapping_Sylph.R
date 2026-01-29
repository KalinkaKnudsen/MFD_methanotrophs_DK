# Script for mapping genome (sylph) abundances on DK-map
# Kalinka Sand Knudsen

setwd("path/to/your/repo/MFD_methanotrophs_DK/")

### Setup env
library(tidyverse)
library(ggrepel)
library(ggpp)
library(ggspatial)
library(sf)
library(grid)
library(rnaturalearth)
library(rnaturalearthdata)
library(vroom)
library(patchwork)
library(ggtext)

## load data 

tax<-readRDS("data/MFD_renamed_tax_25_03_04.rds")%>%
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

# duplicates <- tax %>%
#   group_by(Species) %>%
#   filter(n() > 1) %>%
#   ungroup()%>%
#   filter(type=="GTDB")%>%
#   pull(tip.label)
# 
# tax<-tax%>%filter(!tip.label %in% duplicates)



sylph<-readRDS("output/sylph_gtdb_25_03_06.rds")%>%
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


levels_hab2<-readRDS("output/levels_hab2_sylph.rds")
palette_mfd_hab2<-readRDS("data/palette_mfd_hab2_ISME.rds")

hab2_sort<-sylph%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)


seq_meta <- vroom("data/2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  select(SeqId, fieldsample_barcode)

meta<-readxl::read_excel("data/2025-02-19_mfd_db.xlsx")%>%left_join(seq_meta)%>%select(SeqId, fieldsample_barcode, longitude, latitude)

### Create map of Denmark
points <- meta %>%
  select(fieldsample_barcode, SeqId, longitude, latitude)

world <- st_read(paste0("data/CNTR_RG_01M_2024_4326.shp")) %>% #### Needs to be downloaded from Eurogeographics 
  st_transform(crs = 4326)%>%
  filter(NAME_ENGL%in%c("Denmark", "Sweden", "Germany"))
# 
# ## DK map with points
# map0 <- ggplot(data = world) + 
#   geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
#   geom_point(data = points, aes(x = longitude, y = latitude), size = 4, alpha = 1, shape = 21) +
#   theme_bw(base_size = 16) + 
#   #scale_fill_manual(values = color_vector.9) +
#   annotation_scale(location = "bl", width_hint = 0.5) + 
#   annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
#   coord_sf(xlim = c(7.5, 13), ylim = c(54.5, 58), expand = FALSE) + 
#   labs(fill = "mfd_hab1") +
#   xlab("Longitude") + 
#   ylab("Latitude") + 
#   theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
#         panel.background = element_rect(fill = "aliceblue"),
#         axis.title = element_blank(),
#         legend.position = "right")
# 
# map0



##################################################################
############### Getting methylocystis_mmoX cluster ###############
##################################################################


m_silviterra_clade<-sylph%>%filter(Genus=="g__Methylocystis")%>%filter(Species%in%c("s__Methylocystis silviterrae", "s__Methylocystis hirsuta"))%>%left_join(meta)


m_silviterra_clade_sediment<-m_silviterra_clade%>%
  filter(hab1_label=="Freshwater\nsediment")%>%
  group_by(fieldsample_barcode) %>%
  summarise(
    abundance_sum = sum(Taxonomic_abundance, na.rm = TRUE),          # or use mean(RPKM) if you prefer
    longitude = first(longitude),                # keep one value per group
    latitude  = first(latitude),
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3),
    SeqId     = first(SeqId)
  ) %>%
  ungroup%>%
  arrange(abundance_sum)%>%
  filter(!is.na(longitude))





################## loop start ################

n_samples <- length(unique(m_silviterra_clade_sediment$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.5,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.3,
  text = paste0("*Methylocystis silviterrae*<br>Freshwater sediment (n = ", n_samples, ")")
)


map_m_silviterra <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = m_silviterra_clade_sediment, aes(x = longitude, y = latitude, fill = abundance_sum),               
              size = 3.2, shape = 21, stroke = 0.25,  position = position_jitter(width = 0.04, height = 0.04, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "Taxonomic",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_m_silviterra



unique(m_silviterra_clade_sediment$mfd_hab2)


##########################################
# Only for Urban enclosed water and Rainwater basin


m_silviterra_clade_sediment_urban<-m_silviterra_clade_sediment%>%filter(mfd_hab2%in%c("Urban enclosed water", "Rainwater basin"))%>%arrange(Taxonomic_abundance)


n_samples <- length(unique(m_silviterra_clade_sediment_urban$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.5,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis silviterrae*<br>Urban enclosed water and<br>Rainwater basin (n = ", n_samples, ")")
)


map_m_silviterra_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = m_silviterra_clade_sediment_urban, aes(x = longitude, y = latitude, fill = abundance_sum),               
              size = 3.2, shape = 21, stroke = 0.25,  position = position_jitter(width = 0.04, height = 0.04, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_m_silviterra_urban




# MAG 14 MFD02809.bin.2.311
#Methylocystis rosea MAG2 MFD02809.bin.2.26


m_rosea_clade<-sylph%>%filter(Genus=="g__Methylocystis")%>%filter(user_genome%in%c("MFD02809.bin.2.26"))%>%left_join(meta)%>%filter(mfd_hab2%in%c("Urban enclosed water", "Rainwater basin"))%>%arrange(Taxonomic_abundance)
mag_14_clade<-sylph%>%filter(Genus=="g__Methylocystis")%>%filter(user_genome%in%c("MFD02809.bin.2.311"))%>%left_join(meta)%>%filter(mfd_hab2%in%c("Urban enclosed water", "Rainwater basin"))%>%arrange(Taxonomic_abundance)



n_samples <- length(unique(m_rosea_clade$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.5,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis rosea*<br>Urban enclosed water and<br>Rainwater basin (n = ", n_samples, ")")
)


map_m_rosea_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = m_rosea_clade, aes(x = longitude, y = latitude, fill = Taxonomic_abundance),               
              size = 3.2, shape = 21, stroke = 0.25,  position = position_jitter(width = 0.04, height = 0.04, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_m_rosea_urban



n_samples <- length(unique(mag_14_clade$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.5,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis MAG_14*<br>Urban enclosed water and<br>Rainwater basin (n = ", n_samples, ")")
)


map_m_MAG14_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = mag_14_clade, aes(x = longitude, y = latitude, fill = Taxonomic_abundance),               
              size = 3.2, shape = 21, stroke = 0.25,  position = position_jitter(width = 0.04, height = 0.04, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_m_MAG14_urban




############################################################################
########################## Only Urban enclosed water #######################



m_silviterra_clade_sediment_urban<-m_silviterra_clade_sediment%>%filter(mfd_hab2%in%c("Urban enclosed water"))%>%arrange(abundance_sum)


n_samples <- length(unique(m_silviterra_clade_sediment_urban$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.2,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis silviterrae*<br>Urban enclosed water<br>(n = ", n_samples, ")")
)


map_m_silviterra_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = m_silviterra_clade_sediment_urban, aes(x = longitude, y = latitude, fill = abundance_sum),               
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.02, height = 0.02, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "Tax\nabund.\n[%]",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5, height = unit(0.08, "cm"), line_width = 0.3, text_cex = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "cm"), pad_y = unit(0.5, "cm"), style = north_arrow_fancy_orienteering(text_size = 5, line_width = 0.4),
  #                        height = unit(0.5, "cm"), width = unit(0.5, "cm")) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        plot.margin = margin(0.0, 1, 0, 0),
        axis.text = element_text(size=5),
        legend.title = element_text(size=5, margin = margin(b = 2)),  # small bottom margin
        legend.text = element_text(size=5, margin = margin(l = 1.6, unit = "pt")),
        legend.position = c(0.93, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.ticks.length = unit(2, "pt"), 
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1))

#map_m_silviterra_urban

# 
# ggsave("output/map_m_silviterra_urban_test.png", map_m_silviterra_urban, 
#        #width=15, 
#        #height=7, 
#        units = c("mm"),
#        height = 90,
#        width = 90,
#        dpi=300)




m_rosea_clade<-sylph%>%filter(Genus=="g__Methylocystis")%>%filter(user_genome%in%c("MFD02809.bin.2.26"))%>%left_join(meta)%>%filter(mfd_hab2%in%c("Urban enclosed water"))%>%arrange(Taxonomic_abundance)
mag_14_clade<-sylph%>%filter(Genus=="g__Methylocystis")%>%filter(user_genome%in%c("MFD02809.bin.2.311"))%>%left_join(meta)%>%filter(mfd_hab2%in%c("Urban enclosed water"))%>%arrange(Taxonomic_abundance)



n_samples <- length(unique(m_rosea_clade$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.2,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis rosea*<br>Urban enclosed water<br>(n = ", n_samples, ")")
)


map_m_rosea_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = m_rosea_clade, aes(x = longitude, y = latitude, fill = Taxonomic_abundance),               
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.02, height = 0.02, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "Tax\nabund.\n[%]",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5, height = unit(0.08, "cm"), line_width = 0.3, text_cex = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "cm"), pad_y = unit(0.4, "cm"), style = north_arrow_fancy_orienteering(text_size = 6, line_width = 0.5),
  #                        height = unit(0.8, "cm"), width = unit(0.8, "cm")) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        plot.margin = margin(0.0, 1, 0, 0),
        axis.text = element_text(size=5),
        axis.text.y=element_blank(),
        legend.title = element_text(size=5, margin = margin(b = 2)),  # small bottom margin
        legend.text = element_text(size=5, margin = margin(l = 1.6, unit = "pt")),
        legend.position = c(0.93, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
        axis.ticks.y=element_blank())

#map_m_rosea_urban



n_samples <- length(unique(mag_14_clade$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(m_silviterra_clade_sediment$longitude) - 0.2,
  latitude  = max(m_silviterra_clade_sediment$latitude) + 0.25,
  text = paste0("*Methylocystis MAG_14*<br>Urban enclosed water<br>(n = ", n_samples, ")")
)


map_m_MAG14_urban <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) +    scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+   ggnewscale::new_scale_fill() +
  geom_jitter(data = mag_14_clade, aes(x = longitude, y = latitude, fill = Taxonomic_abundance),               
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.02, height = 0.02, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "Tax\nabund.\n[%]",  # The label
    colors = c(
      "#FFFFFF4C",  # white, alpha = 0.3
      "#82A3CD99",  # #82A3CD, alpha = 0.6
      "#FF8C00CC",  # darkorange, alpha = 0.6
      "#8B0000E6",  # darkred, alpha = 0.9
      "#000000E6"   # black, alpha = 0.9
    ),  
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 3),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5, height = unit(0.08, "cm"), line_width = 0.3, text_cex = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "cm"), pad_y = unit(0.4, "cm"), style = north_arrow_fancy_orienteering(text_size = 6, line_width = 0.5),
  #                        height = unit(0.8, "cm"), width = unit(0.8, "cm")) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = F) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.text = element_text(size=5),
        axis.text.y=element_blank(),
        legend.title = element_text(size=5, margin = margin(b = 2)),  # small bottom margin
        legend.text = element_text(size=5, margin = margin(l = 1.6, unit = "pt")),
        legend.position = c(0.93, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
        axis.ticks.y=element_blank())
#map_m_MAG14_urban




p_combine<- (map_m_silviterra_urban+map_m_rosea_urban+map_m_MAG14_urban) 


ggsave("output/Methylocystis_MAGs_Urban_enclosed_water_25_10_16.png",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 183,
       dpi=300)


ggsave("output/Methylocystis_MAGs_Urban_enclosed_water_25_10_16.svg",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 183,
       dpi=300)



