# Script for analyzing genome (sylph) correlations
# Kalinka Sand Knudsen

setwd("path/to/your/repo/MFD_methanotrophs_DK/")


### Setup env
library(tidyverse)
library(vroom)
library(patchwork)
library(ggtext)

## load data 
OTU_filtered_long<-readRDS("output/OTU_filtered_long_25_09_01.rds")

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



OTU_filtered_long<-merge(OTU_filtered_long, tax_curated, by= "Tax")

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




hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater")

OTU_filtered_long<-OTU_filtered_long%>%
  filter(mfd_hab1 %in% hab_filter)
  

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


sylph<-readRDS("data/sylph_gtdb_25_03_06.rds")%>%
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



##################################
######## Adding Fisher's exact test 
##################################

jaccard_with_fisher <- function(list1, list2) {
  # Convert to presence/absence
  list1_pa <- ifelse(list1 > 0, 1, 0)
  list2_pa <- ifelse(list2 > 0, 1, 0)
  N <- length(list1_pa)
  
  # Compute contingency table counts
  a <- sum(list1_pa == 1 & list2_pa == 1)  # both present
  b <- sum(list1_pa == 1 & list2_pa == 0)  # species1 only
  c <- sum(list1_pa == 0 & list2_pa == 1)  # species2 only
  d <- sum(list1_pa == 0 & list2_pa == 0)  # neither present
  
  # Jaccard metrics
  intersection <- a
  union <- a + b + c
  similarity <- ifelse(union > 0, intersection / union, NA)
  distance <- 1 - similarity
  
  # Expected overlap under independence
  expected_overlap <- ((a + b) * (a + c)) / N
  
  # Fisher tests
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fisher_two_sided <- fisher.test(contingency, alternative = "two.sided")$p.value
  fisher_greater   <- fisher.test(contingency, alternative = "greater")$p.value
  fisher_less      <- fisher.test(contingency, alternative = "less")$p.value
  
  # Return tidy tibble with all counts as separate columns
  tibble(
    similarity = similarity,
    distance = distance,
    intersection = intersection,
    union = union,
    a_both_present = a,
    b_species1_only = b,
    c_species2_only = c,
    d_neither = d,
    expected_overlap = expected_overlap,
    fisher_two_sided_p = fisher_two_sided,
    fisher_greater_p = fisher_greater,
    fisher_less_p = fisher_less
  )
}





###############################################################################
clades_of_alpha<-sylph%>%filter(Genus=="g__Methylocella")%>%
  mutate(clade = if_else(between(genome_number, 1, 23), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 24, 30), "clade_2", clade))



clades_of_alpha_sum<-clades_of_alpha%>%
  group_by(fieldsample_barcode, clade) %>%
  summarise(
    abundance_sum = sum(Taxonomic_abundance, na.rm = TRUE),          # or use mean(RPKM) if you prefer
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3),
    SeqId     = first(SeqId)
  ) %>%
  ungroup%>%
  pivot_wider(names_from = clade, values_from = abundance_sum)





p_combine<-readRDS("data/palette_mfd_hab2_ISME.rds") 

P <-clades_of_alpha_sum %>%
  filter(!(clade_1 == 0 & clade_2 == 0))  %>%
  ggplot(., aes(
    y = clade_1,
    x= clade_2)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine)+ #name = paste0("Bogs, Mires and Fens,\nn = ",num_samples)
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
#  theme_classic() +
 # coord_cartesian(xlim = c(0, max(Methylocella_Forest$`Methylocella_pxmA`)), ylim = c(0,max(Methylocella_Forest$Methylocella_USCa)))+
  theme(#legend.position = c(0.84,0.15),
        legend.background = element_rect(fill = "transparent"))+
  facet_wrap(.~mfd_hab1)+
  guides(fill = guide_legend(ncol = 2)) + 
  ggpubr::stat_cor(method=c("spearman"), r.digits = 3,  cor.coef.name = c("rho"))

P



P <-clades_of_alpha_sum %>%
  filter(mfd_hab1%in%c("Grassland formations", "Forests", "Temperate heath and scrub", "Dunes"))%>%
  ggplot(., aes(
    y = clade_1,
    x= clade_2)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine)+ #name = paste0("Bogs, Mires and Fens,\nn = ",num_samples)
 # geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  theme_bw() +
  # coord_cartesian(xlim = c(0, max(Methylocella_Forest$`Methylocella_pxmA`)), ylim = c(0,max(Methylocella_Forest$Methylocella_USCa)))+
  theme(#legend.position = c(0.84,0.15),
    legend.background = element_rect(fill = "transparent"))+
  facet_wrap(.~mfd_hab1)+
  guides(fill = guide_legend(ncol = 2)) + 
  ggpubr::stat_cor(method=c("spearman"), r.digits = 3,  cor.coef.name = c("rho"))

P




############# making Jaccard #############
# Create Function for Calculating Jaccard Similarity and Distance
jaccard <- function(list1, list2) {
  # Convert to presence/absence (0/1)
  list1_pa <- ifelse(list1 > 0, 1, 0)
  list2_pa <- ifelse(list2 > 0, 1, 0)
  
  # Check equal length
  if (length(list1_pa) != length(list2_pa)) {
    stop("Both lists must have the same length.")
  }
  
  # Compute intersection and union as counts
  intersection <- sum(list1_pa & list2_pa)
  union <- sum(list1_pa | list2_pa)
  
  similarity <- intersection / union
  distance <- 1 - similarity
  
  
  # Combine vectors into a dataframe for inspection
  comparison_df <- data.frame(
    list1 = list1,
    list1_pa = list1_pa,
    list2 = list2,
    list2_pa = list2_pa,
    intersection_row = as.integer(list1_pa & list2_pa),
    union_row = as.integer(list1_pa | list2_pa)
  )
  
  return(list(
    similarity = similarity,
    distance = distance,
    intersection = intersection,
    union = union,
    comparison = comparison_df
  ))
}


#### Testing if it changes 
clades_of_alpha_sum_non_zero <-clades_of_alpha_sum%>%
  filter(!(clade_1 == 0 & clade_2 == 0))

res1<-jaccard(clades_of_alpha_sum_non_zero$clade_1, clades_of_alpha_sum_non_zero$clade_2)
res2<-jaccard(clades_of_alpha_sum$clade_1, clades_of_alpha_sum$clade_2)

res1$similarity
res1$distance

##################################
######## Adding Fisher's exact test 
##################################

jaccard_with_fisher <- function(list1, list2) {
  # Convert to presence/absence
  list1_pa <- ifelse(list1 > 0, 1, 0)
  list2_pa <- ifelse(list2 > 0, 1, 0)
  N <- length(list1_pa)
  
  # Compute contingency table counts
  a <- sum(list1_pa == 1 & list2_pa == 1)  # both present
  b <- sum(list1_pa == 1 & list2_pa == 0)  # species1 only
  c <- sum(list1_pa == 0 & list2_pa == 1)  # species2 only
  d <- sum(list1_pa == 0 & list2_pa == 0)  # neither present
  
  # Jaccard metrics
  intersection <- a
  union <- a + b + c
  similarity <- ifelse(union > 0, intersection / union, NA)
  distance <- 1 - similarity
  
  # Expected overlap under independence
  expected_overlap <- ((a + b) * (a + c)) / N
  
  # Fisher tests
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fisher_two_sided <- fisher.test(contingency, alternative = "two.sided")$p.value
  fisher_greater   <- fisher.test(contingency, alternative = "greater")$p.value
  fisher_less      <- fisher.test(contingency, alternative = "less")$p.value
  
  # Return tidy tibble with all counts as separate columns
  tibble(
    similarity = similarity,
    distance = distance,
    intersection = intersection,
    union = union,
    a_both_present = a,
    b_species1_only = b,
    c_species2_only = c,
    d_neither = d,
    expected_overlap = expected_overlap,
    fisher_two_sided_p = fisher_two_sided,
    fisher_greater_p = fisher_greater,
    fisher_less_p = fisher_less
  )
}


df<-clades_of_alpha_sum
species_cols <- grep("clade|other", names(df), value = TRUE)
habitats <- unique(df$mfd_hab1)

results_list <- list()

for (hab in habitats) {
  subset_df <- df %>% filter(mfd_hab1 == hab)
  mat <- subset_df[, species_cols]
  
  # Unique unordered pairs
  combs <- combn(species_cols, 2, simplify = FALSE)
  
  for (pair in combs) {
    res <- jaccard_with_fisher(mat[[pair[1]]], mat[[pair[2]]])
    
    # Add metadata
    res$species1 <- pair[1]
    res$species2 <- pair[2]
    res$mfd_hab1 <- hab
    
    results_list <- append(results_list, list(res))
  }
}

# Combine all into one tidy dataframe
results <- bind_rows(results_list)%>%
  mutate(
    fisher_two_sided_p_BH = p.adjust(fisher_two_sided_p, method = "BH"),
    fisher_greater_p_BH   = p.adjust(fisher_greater_p, method = "BH"),
    fisher_less_p_BH      = p.adjust(fisher_less_p, method = "BH")
  )%>%relocate(mfd_hab1, fisher_less_p, fisher_less_p_BH)

results%>%filter(!grepl("other", species2))%>%arrange(similarity)




############################################################################
############################################################################
############ Methylocystis clades comparison ###############################

clades_of_cystis<-sylph%>%filter(Genus=="g__Methylocystis")%>%
  mutate(clade = if_else(between(genome_number, 1, 19), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 27, 33), "clade_2", clade))



clades_of_cystis_sum<-clades_of_cystis%>%
  group_by(fieldsample_barcode, clade) %>%
  summarise(
    abundance_sum = sum(Taxonomic_abundance, na.rm = TRUE),          
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3),
    SeqId     = first(SeqId)
  ) %>%
  ungroup%>%
  pivot_wider(names_from = clade, values_from = abundance_sum)




df<-clades_of_cystis_sum%>%filter(mfd_hab1=="Bogs, mires and fens" |
                                  mfd_hab2 %in% c("Semi-natural humid meadows",
                                                  "Alluvial woodland",
                                                  "Bog woodland"))
species_cols <- grep("clade|other", names(df), value = TRUE)
habitats <- unique(df$mfd_hab1)
habitats2 <- unique(df$mfd_hab2)

results_list <- list()

for (hab in habitats) {
  subset_df <- df %>% filter(mfd_hab1 == hab)
  mat <- subset_df[, species_cols]
  
  # Unique unordered pairs
  combs <- combn(species_cols, 2, simplify = FALSE)
  
  for (pair in combs) {
    res <- jaccard_with_fisher(mat[[pair[1]]], mat[[pair[2]]])
    
    # Add metadata
    res$species1 <- pair[1]
    res$species2 <- pair[2]
    res$mfd_hab1 <- hab
    
    results_list <- append(results_list, list(res))
  }
}

habitats3 <- unique(df$mfd_hab3)
results_list <- list()
for (hab in habitats3) {
  subset_df <- df %>% filter(mfd_hab3 == hab)
  mat <- subset_df[, species_cols]
  
  # Unique unordered pairs
  combs <- combn(species_cols, 2, simplify = FALSE)
  
  for (pair in combs) {
    res <- jaccard_with_fisher(mat[[pair[1]]], mat[[pair[2]]])
    
    # Add metadata
    res$species1 <- pair[1]
    res$species2 <- pair[2]
    res$mfd_hab3 <- hab
    
    results_list <- append(results_list, list(res))
  }
}

# Combine all into one tidy dataframe
results <- bind_rows(results_list)%>%
  mutate(
    fisher_two_sided_p_BH = p.adjust(fisher_two_sided_p, method = "BH"),
    fisher_greater_p_BH   = p.adjust(fisher_greater_p, method = "BH"),
    fisher_less_p_BH      = p.adjust(fisher_less_p, method = "BH")
  )%>%relocate(mfd_hab3, fisher_less_p, fisher_less_p_BH)

results%>%filter(!grepl("other", species2))%>%arrange(similarity)

