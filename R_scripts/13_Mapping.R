# Script for mapping gene abundances on DK-map
# Kalinka Sand Knudsen
.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
setwd("~/scripts/MFD/methanotrophs/R_scripts/output/")


### Setup env
library(tidyverse)
library(ggplot2)
library(ggrepel)
#library(ggpp)
library(ggspatial)
library(sf)
#library(grid)
#library(rnaturalearth)
#library(rnaturalearthdata)
library(vroom)
library(patchwork)
library(ggtext)

## load data 
OTU_filtered_long<-readRDS("OTU_filtered_long_24_06_24.rds")
seq_meta <- vroom("../2023-10-11_samples_minimal_metadata_collapsed.csv", delim = ",") %>%
  mutate(flat_name=gsub(".fastq.gz","", flat_name))%>%
  rename(SeqId = flat_name) %>%
  select(SeqId, fieldsample_barcode)

meta<-readxl::read_excel("../2025-02-19_mfd_db.xlsx")%>%left_join(seq_meta)



### Create map of Denmark
points <- meta %>%
  select(fieldsample_barcode, SeqId, longitude, latitude, mfd_hab1)

## World data
#world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf", country=c("Denmark", "Sweden", "Germany"))


world <- st_read(paste0("../dataframes/CNTR_RG_01M_2024_4326.shp/CNTR_RG_01M_2024_4326.shp")) %>%
  st_transform(crs = 4326)%>%
  filter(NAME_ENGL%in%c("Denmark", "Sweden", "Germany"))

## DK map with points
map0 <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_point(data = points, aes(x = longitude, y = latitude, fill = mfd_hab1), size = 4, alpha = 1, shape = 21) +
  theme_bw(base_size = 16) + 
#  scale_fill_manual(values = color_vector.9) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  labs(fill = "mfd_hab1") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = "right")

map0



long_point<-7.95168
lat_point<-57.96323




points_pmoA2<-OTU_filtered_long%>%filter(type=="pmoA2")%>%filter(Tax_short=="Methylocystis_pmoA2")%>%left_join(meta)%>%arrange(RPKM)

points_pmoA2_heath<-points_pmoA2%>%filter(mfd_hab3=="Wet heath")

### Methylocystis_pmoA2 (only relevant one)
# Map in the habs: Wet heath, Molina meadows, Decalcified Empetrum dunes, humid dune slacks


### USCg in: calcareous grasslands and Xeric sand calcereous grasslands






################## loop start ################

n_samples <- length(unique(points_pmoA2_heath$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(points_pmoA2_heath$longitude) - 0.2,
  latitude  = max(points_pmoA2_heath$latitude) + 0.4,
  text = paste0("*Methylocystis pmoA2*<br>Wet heath<br>(n = ", n_samples, ")")
)


map_heath <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) + 
  scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+
  ggnewscale::new_scale_fill() +
  geom_jitter(data = points_pmoA2_heath, aes(x = longitude, y = latitude, fill = RPKM),               
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.05, height = 0.05, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        plot.margin = margin(0, 2, 0, 0),
        axis.text = element_text(size=5),
        #axis.text.y=element_blank(),
        legend.title = element_text(size=5, margin = margin(b = 2)),  # small bottom margin
        legend.text = element_text(size=5, margin = margin(l = 1.6, unit = "pt")),
        legend.position = c(0.9, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
       # axis.ticks.y=element_blank()
       )
#map_heath



points_pmoA2_dune_slack<-points_pmoA2%>%filter(mfd_hab3=="Humid dune slacks")
n_samples <- length(unique(points_pmoA2_dune_slack$fieldsample_barcode))


label_df <- data.frame(
  longitude = min(points_pmoA2_heath$longitude) - 0.2,
  latitude  = max(points_pmoA2_heath$latitude) + 0.4,
  text = paste0("*Methylocystis pmoA2*<br>Humid dune slacks<br>(n = ", n_samples, ")")
)





map_dune_slack <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) + 
  scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+
  ggnewscale::new_scale_fill() +
  geom_jitter(data = points_pmoA2_dune_slack, aes(x = longitude, y = latitude, fill = RPKM), 
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.05, height = 0.05, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
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
        legend.position = c(0.9, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
        axis.ticks.y=element_blank())

#map_dune_slack
# 
# ggsave("maps/pmoA2_wet_heath.png",
#        map_heath,
#        height = 10,
#        width = 10)


#


p_combine<- (map_heath+map_dune_slack)


ggsave("maps/Methylocystis_pmoA2_heath_dune_25_10_16.png",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 100,
       dpi=300)


#
ggsave("maps/Methylocystis_pmoA2_heath_dune_25_10_16.svg",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 120,
       dpi=300)




points_pmoA2_Molinia_meadows<-points_pmoA2%>%filter(mfd_hab3=="Molinia meadows")


n_samples <- length(unique(points_pmoA2_Molinia_meadows$fieldsample_barcode))


label_df <- data.frame(
  longitude = min(points_pmoA2_heath$longitude) - 0.5,
  latitude  = max(points_pmoA2_heath$latitude) + 0.3,
  text = paste0("*Methylocystis pmoA2*<br>Molinia meadows (n = ", n_samples, ")")
)

map_meadow <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_jitter(data = points_pmoA2_Molinia_meadows, aes(x = longitude, y = latitude, fill = RPKM),
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_meadow




# ggsave("maps/pmoA2_Molinia_meadows.png",
#        map_meadow,
#        height = 10,
#        width = 10)








points_pmoA2_Decalcified_Empetrum_dunes<-points_pmoA2%>%filter(mfd_hab3=="Decalcified Empetrum dunes")
n_samples <- length(unique(points_pmoA2_Decalcified_Empetrum_dunes$fieldsample_barcode))

label_df <- data.frame(
  longitude = min(points_pmoA2_heath$longitude) - 0.5,
  latitude  = max(points_pmoA2_heath$latitude) + 0.3,
  text = paste0("*Methylocystis pmoA2*<br>Decalcified Empetrum dunes (n = ", n_samples, ")")
)


map_Decalcified_Empetrum_dunes <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_jitter(data = points_pmoA2_Decalcified_Empetrum_dunes, aes(x = longitude, y = latitude, fill = RPKM),              
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())


map_Decalcified_Empetrum_dunes



# ggsave("maps/pmoA2_Decalcified_Empetrum_dunes.png",
#        map_Decalcified_Empetrum_dunes,
#        height = 10,
#        width = 10)



points_pmoA2_dune_slack<-points_pmoA2%>%filter(mfd_hab3=="Humid dune slacks")
n_samples <- length(unique(points_pmoA2_dune_slack$fieldsample_barcode))


label_df <- data.frame(
  longitude = min(points_pmoA2_heath$longitude) - 0.5,
  latitude  = max(points_pmoA2_heath$latitude) + 0.3,
  text = paste0("*Methylocystis pmoA2*<br>Humid dune slacks (n = ", n_samples, ")")
)




map_dune_slack <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_jitter(data = points_pmoA2_dune_slack, aes(x = longitude, y = latitude, fill = RPKM), 
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.86),
        legend.title=element_text(size=12),
        legend.background = element_blank())

map_dune_slack



# ggsave("maps/pmoA2_dune_slack.png",
#        map_dune_slack,
#        height = 10,
#        width = 10)



#




p_combine<- (map_heath+map_meadow) / (map_Decalcified_Empetrum_dunes+map_dune_slack)


ggsave("maps/pmoA2_combined.png",
       p_combine,
       height = 16,
       width = 16)












#######################################################################################
############################ USCg #####################################################
#######################################################################################


points_USCg<-OTU_filtered_long%>%filter(type=="pmoA")%>%filter(grepl("Root; USCg", Tax))%>%left_join(meta)%>%arrange(RPKM)


points_summary <- points_USCg %>%
  group_by(fieldsample_barcode) %>%
  summarise(
    RPKM_sum = sum(RPKM, na.rm = TRUE),          # or use mean(RPKM) if you prefer
    longitude = first(longitude),                # keep one value per group
    latitude  = first(latitude),
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3)
  ) %>%
  ungroup%>%
  arrange(RPKM_sum)



points_USCg_calcareous<-points_summary%>%
  filter(mfd_hab3%in%c("Calcareous grassland")) #"Xeric sand calcareous grasslands"

n_samples <- length(unique(points_USCg_calcareous$fieldsample_barcode))


label_df <- data.frame(
  longitude = long_point,
  latitude  = lat_point,
  text = paste0("*USCg*<br>Calcareous grassland<br>(n = ", n_samples, ")")
)




map_calcareous <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) + 
  scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+
  ggnewscale::new_scale_fill() +
  geom_jitter(data = points_USCg_calcareous, aes(x = longitude, y = latitude, fill = RPKM_sum), 
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.05, height = 0.05, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        plot.margin = margin(0, 2, 0, 0),
        axis.text = element_text(size=5),
   #     axis.text.y=element_blank(),
        legend.title = element_text(size=5, margin = margin(b = 2)),  # small bottom margin
        legend.text = element_text(size=5, margin = margin(l = 1.6, unit = "pt")),
        legend.position = c(0.9, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
    #    axis.ticks.y=element_blank()
   )

map_calcareous



# ggsave("maps/pmoA2_dune_slack.png",
#        map_dune_slack,
#        height = 10,
#        width = 10)




points_USCg_xeric<-points_summary%>%left_join(meta)%>%
  filter(mfd_hab3%in%c("Xeric sand calcareous grasslands")) #""

n_samples <- length(unique(points_USCg_xeric$fieldsample_barcode))


label_df <- data.frame(
  longitude = long_point,
  latitude  = lat_point,
  text = paste0("*USCg*<br>Xeric sand calcareous grasslands<br>(n = ", n_samples, ")")
)



map_xeric <- ggplot(data = world) + 
  geom_sf(aes(fill=NAME_ENGL), linewidth=0.2, show.legend = FALSE) + 
  scale_fill_manual(values = c("Denmark" = "antiquewhite", "Sweden" = "grey90","Germany" = "grey90"))+
  ggnewscale::new_scale_fill() +
  geom_jitter(data = points_USCg_xeric, aes(x = longitude, y = latitude, fill = RPKM_sum), 
              size = 1.2, shape = 21, stroke = 0.1,  position = position_jitter(width = 0.05, height = 0.05, seed = 5)) +
  geom_richtext(
    data = label_df, aes(x = longitude, y = latitude, label = text), 
    hjust = 0, vjust = 1, size = 5/.pt, fill=NA,  label.color = NA  )+
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
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  # annotation_scale(location = "bl", width_hint = 0.5) + 
  # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.8, 13.15), ylim = c(54.5, 58), expand = FALSE) + 
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
        legend.position = c(0.9, 0.86),
        legend.background = element_blank(), 
        legend.key.size = unit(0.2, "cm"),
        axis.line = element_line(color = "black", linewidth = 0.1),
        axis.ticks = element_line(color = "black", linewidth = 0.1),
        axis.ticks.length = unit(2, "pt"), 
        axis.ticks.y=element_blank())






p_combine<- map_calcareous + map_xeric


ggsave("maps/USCg_combined_25_10_16.png",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 100,
       dpi=300)


ggsave("maps/USCg_combined_25_10_16.svg",
       p_combine,
       units = c("mm"),
       height = 90,
       width = 100,
       dpi=300)











### Create DK map with repelled points
map1 <- ggplot(data = world) + 
  geom_sf(fill = "antiquewhite") + 
  geom_point(data = points_pmoA2_heath, aes(x = longitude, y = latitude), size = 1, color = "black") +
  geom_point_s(data = points_pmoA2_heath, aes(x = longitude, y = latitude, fill = RPKM), 
               size = 3, alpha = 1, shape = 21,
               position = position_nudge_center(x = 0.1,
                                                y = 0.1,
                                                direction = "radial")) +
  theme_bw(base_size = 16) + 
  scale_fill_gradientn(
    name = "RPKM",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(7.5, 13.5), ylim = c(54.5, 58), expand = FALSE) + 
  labs(fill = "Ecospace") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        legend.position = "right")

map1
