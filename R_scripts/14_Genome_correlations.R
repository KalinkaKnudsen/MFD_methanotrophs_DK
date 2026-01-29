# Script for mapping gene abundances on DK-map
# Kalinka Sand Knudsen
.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
setwd("~/scripts/MFD/methanotrophs/R_scripts/output/")


### Setup env
library(tidyverse)
library(vroom)
library(patchwork)
library(ggtext)

## load data 
OTU_filtered_long<-readRDS("OTU_filtered_long_25_09_01.rds")

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
  

tax<-readRDS("../MFD_renamed_tax_25_03_04.rds")%>%
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


sylph<-readRDS("genome_quantification/sylph_gtdb_25_03_06.rds")%>%
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


levels_hab2<-readRDS("genome_quantification/levels_hab2_sylph_DMS.rds")
palette_mfd_hab2<-readRDS("./palette_mfd_hab2_ISME.rds")

hab2_sort<-sylph%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)


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





p_combine<-readRDS("../palette_mfd_hab2_ISME.rds") 

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



## A little test below to find the number of presence/absence in forest types
# 
# clades_of_alpha_sum %>%
#   filter(mfd_hab1 == "Forests") %>%
#   mutate(
#     clade_1_pa = ifelse(clade_1 > 0, 1, 0),
#     clade_2_pa = ifelse(clade_2 > 0, 1, 0)
#   ) %>%
#   group_by(mfd_hab2) %>%
#   summarise(
#     clade_1_present = sum(clade_1_pa, na.rm = TRUE),
#     clade_2_present = sum(clade_2_pa, na.rm = TRUE),
#     total_rows = n()
#   ) %>%
#   ungroup()
#   





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


########## pmoA2 and the second clade ######


#Interestingly, we observe that the distribution of this second clade correlates
#with the pmoA2 of Methylocystis which was found almost exclusively in acidic bogs 
#and only in lower abundance in few dune, bog woodland, and grassland and heath samples 



points_pmoA2<-OTU_filtered_long%>%filter(type=="pmoA2")%>%
  filter(Tax_short=="Methylocystis_pmoA2")%>%
    mutate(pmoA2_pa = ifelse(RPKM > 0, 1, 0))%>%
  select(pmoA2_pa, SeqId)
  



df<-clades_of_cystis_sum%>%left_join(points_pmoA2)%>%filter(!is.na(pmoA2_pa))%>%
                          filter(mfd_hab2=="Sphagnum acid bogs" |
                                 mfd_hab3 %in% c("Humid dune slacks", # "Alluvial woodland",
                                                  "Bog woodland",
                                                  "Wet heath", "Molinia meadows"))

species_cols <- grep("clade|other|pmo", names(df), value = TRUE)
habitats <- unique(df$mfd_hab2)

results_list <- list()

for (hab in habitats) {
  subset_df <- df %>% filter(mfd_hab2 == hab)
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
  )%>%relocate(mfd_hab1, fisher_greater_p, fisher_greater_p_BH)





results%>%filter(grepl("pmo", species2) & grepl("clade_2", species1))%>%arrange(similarity)


### I think this might work better with correlations tho????


points_pmoA2<-OTU_filtered_long%>%filter(type=="pmoA2")%>%
  filter(Tax_short=="Methylocystis_pmoA2")%>%
 # mutate(pmoA2_pa = ifelse(RPKM > 0, 1, 0))%>%
  select(RPKM, SeqId)




df<-clades_of_cystis_sum%>%left_join(points_pmoA2)%>%filter(!is.na(RPKM))%>%
  filter(mfd_hab2=="Sphagnum acid bogs" |
           mfd_hab3 %in% c("Humid dune slacks","Decalcified Empetrum dunes", # "Alluvial woodland",
                           "Bog woodland",
                           "Wet heath", "Molinia meadows"))





# Get unique mfd_hab1 groups
habitats <- df%>%filter(!mfd_hab1=="Freshwater")%>%
  pull(mfd_hab1)%>%unique()

# Initialize empty list to store plots
plot_list <- list()

# Loop over habitats
for (h in habitats) {
  p <- df %>%
    #filter(!(Rhodomicrobium_pmoA == 0 & Rhodomicrobium_mmoX == 0))  %>% Removing the double-zeroes here influences the correlation like crazy!
    filter(!is.na(mfd_hab2), mfd_hab1 == h) %>%
    ggplot(aes(
      y = clade_2,
      x = RPKM)) +
    geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6,
               color = "black", pch = 21, stroke=.8) +
  #  scale_fill_manual(values = p_combine) +
    geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
    labs(
      y = "Methylocystis_clade_2 [Rel. Abund]",
      x = "Methylocystis_pmoA2 [RPKM]",
      title = h
    ) +
    coord_cartesian(
      xlim = c(0, max(df$clade_2, na.rm=TRUE)),
      ylim = c(0, max(df$clade_2, na.rm=TRUE))
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.justification = c(0, 1), 
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    ggpubr::stat_cor(label.x = 0, label.y=12, method="spearman", r.digits = 3, p.digits = 0.00001, cor.coef.name = "rho")
  
  plot_list[[h]] <- p
}

# Combine with patchwork
P_all <- patchwork::wrap_plots(plot_list)
P_all


ggsave("correlations/Methylocystis_clade_2_pmoA2.png",
       P_all,
       height = 9,
       width = 12,
       dpi=400)

#












p_combine


########################## Now I want to plot it in one "big" plot with one line per habitat ########################

df<-clades_of_cystis_sum%>%left_join(points_pmoA2)%>%filter(!is.na(RPKM))%>%
  filter(mfd_hab2=="Sphagnum acid bogs" |
           mfd_hab3 %in% c("Humid dune slacks","Decalcified Empetrum dunes", # "Alluvial woodland",
                           "Bog woodland",
                           "Wet heath", "Molinia meadows"))


min_double <- .Machine$double.xmin  # R's smallest double-precision value

format_p <- function(p) {
  ifelse(p < min_double,
         paste0("< ", formatC(min_double, format = "e", digits = 2)),
         formatC(p, format = "e", digits = 2))
}

habitats<-unique(df$mfd_hab2)
results_list <- list()
for (hab in habitats) {
  subset_df <- df %>% filter(mfd_hab2 == hab)
  mat <- subset_df[, species_cols]
  
  # Unique unordered pairs
  combs <- combn(species_cols, 2, simplify = FALSE)
  
  for (pair in combs) {
    res <- jaccard_with_fisher(mat[[pair[1]]], mat[[pair[2]]])
    
    # Add metadata
    res$species1 <- pair[1]
    res$species2 <- pair[2]
    res$mfd_hab2 <- hab
    
    results_list <- append(results_list, list(res))
  }
}


# Combine all into one tidy dataframe
results <- bind_rows(results_list)%>%
  mutate(
    fisher_two_sided_p_BH = p.adjust(fisher_two_sided_p, method = "BH"),
    fisher_greater_p_BH   = p.adjust(fisher_greater_p, method = "BH"),
    fisher_less_p_BH      = p.adjust(fisher_less_p, method = "BH")
  )%>%relocate(mfd_hab2, fisher_greater_p, fisher_greater_p_BH)
results_filt<-results%>%filter(species1=="clade_2",species2=="RPKM")%>%select(mfd_hab2, fisher_greater_p_BH, similarity)

cor_stats <- df %>%
  left_join(results_filt,by = join_by(mfd_hab2))%>%
  group_by(mfd_hab2) %>%
  summarise(n_samples = n(),
            rho_val = cor(RPKM, clade_2, method = "spearman", use = "complete.obs"),
            rho_p   = cor.test(RPKM, clade_2, method = "spearman")$p.value,
            lm_fit  = list(lm(clade_2 ~ RPKM)),
            coef_val = coef(lm(clade_2 ~ RPKM))[2],
            coef_p   = summary(lm(clade_2 ~ RPKM))$coefficients[2, 4],
            r2       = summary(lm(clade_2 ~ RPKM))$r.squared,
            J        = first(similarity),  
            fisher_greater_p_BH = first(fisher_greater_p_BH)) %>%
  mutate(rho_p_fmt   = format_p(rho_p),
         coef_p_fmt  = format_p(coef_p),
         J_fmt      = round(J, 2),
         fisher_fmt = if_else(fisher_greater_p_BH==1, paste0("1*"), format_p(fisher_greater_p_BH )))%>%
  mutate(label = paste0(mfd_hab2, " (n = ", n_samples, ")\n",
                        "\u03C1: ", round(rho_val, 2), " (p=", rho_p_fmt, ")\n",
                        "R²: ", round(r2, 2), " (p=", coef_p_fmt, ")\n",
                        "J = ", J_fmt, " (p=", fisher_fmt, ")"
                        # "Reg. coeff.: ", round(coef_val, 3), " (p=", coef_p_fmt, ") ",
  )
  )


# cor_stats <- df %>%
#   group_by(mfd_hab2) %>%
#   summarise(n_samples = n(),
#     rho_val = cor(RPKM, clade_2, method = "spearman", use = "complete.obs"),
#     rho_p   = cor.test(RPKM, clade_2, method = "spearman")$p.value,
#     lm_fit  = list(lm(clade_2 ~ RPKM)),
#     coef_val = coef(lm(clade_2 ~ RPKM))[2],
#     coef_p   = summary(lm(clade_2 ~ RPKM))$coefficients[2, 4],
#     r2       = summary(lm(clade_2 ~ RPKM))$r.squared) %>%
#   mutate(rho_p_fmt   = format_p(rho_p),
#     coef_p_fmt  = format_p(coef_p),
#     label = paste0(mfd_hab2, " (n = ", n_samples, ")\n",
#      # "Reg. coeff.: ", round(coef_val, 3), " (p=", coef_p_fmt, ") ",
#       "R²: ", round(r2, 2), " (p=", coef_p_fmt, ")\n",
#       "\u03C1: ", round(rho_val, 2), " (p=", rho_p_fmt, ")"
#     )
#   )
# 
# 

legend_labels <- setNames(cor_stats$label, cor_stats$mfd_hab2)


# Fade line colors for transparency
hab_colors_faded <- scales::alpha(hab_colors, 0.8)

line_data <- df %>%
  group_by(mfd_hab2) %>%
  do({
    model <- lm(clade_2 ~ RPKM, data = .)
    data.frame(
      RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200),
      clade_2 = predict(model, newdata = data.frame(RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200))),
      mfd_hab2 = unique(.$mfd_hab2)
    )
  })


# Plot
p_all <- ggplot(df, aes(x = RPKM, y = clade_2)) +
  # Transparent regression lines
  geom_line(data = line_data, aes(x = RPKM, y = clade_2, color = mfd_hab2, group = mfd_hab2),
            linetype = 1, linewidth = 0.4, show.legend = F)+
  
  # Points with outline-only shapes and bold colors
  geom_point(aes(fill = mfd_hab2), shape = 21,
             size = 0.9, alpha = 0.7, stroke = 0.2) +
  
  # Manual scales
 # scale_linetype_manual(values = c("dotted" = "33")) +
  scale_color_manual(values = hab_colors, labels = legend_labels) +
  scale_fill_manual(values = hab_colors, labels = legend_labels) +
#  scale_shape_manual(values = c(0, 1, 2, 5, 8, 6, 4), labels = legend_labels) +
  
  # Labels and layout
  labs(
    x = "Methylocystis_pmoA2 [RPKM]",
    y = "Methylocystis_clade_2 [Rel. Abund]",
  #  title = "Habitat-wise correlation between pmoA2 and clade_2"
  ) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
  theme_classic() +
  theme(
    legend.spacing.y = unit(1, "cm"),
    legend.spacing.x = unit(0.07, "cm"),
    legend.position = c(0.55, 0.7),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5,margin = margin(t = 2, b = 2, l=2)),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.key.size = unit(0.1, "cm"),
    plot.title = element_text(size = 5, face = "bold", vjust = 1)
  ) +
  guides(
      shape = guide_legend(override.aes = list(size = 0.8, stroke = 0.5, alpha=0.9))
  )





# ggsave("correlations/Methylocystis_clade_2_pmoA2_25_10_16.png",
#        p_all,
#        units = c("mm"),
#        height = 70,
#        width = 80,
#        dpi=300)
# 
# ggsave("correlations/Methylocystis_clade_2_pmoA2_25_10_16.svg",
#        p_all,
#        units = c("mm"),
#        height = 70,
#        width = 80,
#        dpi=300)
# 





# Optional: define custom colors for habitats
hab_colors <- c(
  "Sphagnum acid bogs" = "#32849f",
  "Sea dunes" = "#33ccff",
  "Bog woodland" = "#cdad00",
  "Temperate heath" = "red3",
  "Semi-natural humid meadows" = "#FF7F00"
)

# Plot


library(ggtext)


p_all <- ggplot(df, aes(x = RPKM, y = clade_2)) +
  # Transparent regression lines
  geom_line(data = line_data, aes(x = RPKM, y = clade_2, color = mfd_hab2, group = mfd_hab2),
            linetype = 1, linewidth = 0.4, show.legend = F)+
  geom_point(aes(fill = mfd_hab2), shape = 21,
             size = 1.1, alpha = 0.66, stroke = 0.3) +
  scale_color_manual(values = hab_colors, labels = legend_labels) +
  scale_fill_manual(values = hab_colors, labels = legend_labels) +
  labs(
    x = "Methylocystis pmoA2 [RPKM]",
    y = "Methylocystis clade 2 [Rel. Tax. Abund]",
    #  title = "Habitat-wise correlation between pmoA2 and clade_2"
  ) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
  theme_classic() +
  theme(
    legend.spacing.y = unit(1, "cm"),
    legend.spacing.x = unit(0.07, "cm"),
    legend.position = c(0.63, 0.7),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5,margin = margin(t = 2, b = 2, l=2)),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.key.size = unit(0.1, "cm"),
    plot.title = element_text(size = 5, face = "bold", vjust = 1)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 1.5, stroke = 0.4, alpha = 0.9))
  )




ggsave("correlations/Methylocystis_clade_2_pmoA2_wide_25_11_27.png",
       p_all,
       units = c("mm"),
       height = 70,
       width = 100,
       dpi=300)

ggsave("correlations/Methylocystis_clade_2_pmoA2_wide_25_11_27.svg",
       p_all,
       units = c("mm"),
       height = 70,
       width = 100,
       dpi=300)




####################################################################################
######################## Below not in use ##########################################
####################################################################################








####################################################################################
########### Making comparisons pmoA2 Methylocystis 2025-11-27
####################################################################################


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



points_pmoA2<-OTU_filtered_long%>%filter(type=="pmoA2")%>%
  filter(Tax_short=="Methylocystis_pmoA2")%>%
  # mutate(pmoA2_pa = ifelse(RPKM > 0, 1, 0))%>%
  select(RPKM, SeqId, mfd_hab3, mfd_hab2)


points_pmoA2%>%arrange(desc(RPKM))

mfd_hab3_filter<-c("Decalcified Empetrum dunes", "Humid dune slacks", "Hydrophilous tall-herb swamp", "Molinia meadows", "Wet heath", "Bog woodland", "Degraded raised bog", "Quaking bogs", "Active raised bogs")
#mfd_hab2_filter<-c("Standing freshwater, lake", "Urban enclosed water")

df <- clades_of_cystis_sum %>%
  left_join(points_pmoA2) %>%
  filter(!is.na(RPKM)) %>%
  filter(
  #  mfd_hab2 %in% mfd_hab2_filter |
      mfd_hab3 %in% mfd_hab3_filter
  )


habitats<-unique(df$mfd_hab3)
# Initialize empty list to store plots





species_cols <- grep("clade|RPKM", names(df), value = TRUE)
results_list <- list()

for (hab in habitats) {
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
  )%>%relocate(mfd_hab3, fisher_greater_p, fisher_greater_p_BH)




plot_list <- list()

# Loop over habitats
for (h in habitats) {
  p <- df %>%
    #filter(!(Rhodomicrobium_pmoA == 0 & Rhodomicrobium_mmoX == 0))  %>% Removing the double-zeroes here influences the correlation like crazy!
    filter(mfd_hab3 == h) %>%
    ggplot(aes(
      y = clade_2,
      x = RPKM)) +
    geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6,
               color = "black", pch = 21, stroke=.8) +
    #  scale_fill_manual(values = p_combine) +
    geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
    labs(
      y = "Methylocystis_clade_2 [Rel. Abund]",
      x = "Methylocystis_pmoA2 [RPKM]",
      title = h
    ) +
    coord_cartesian(
      xlim = c(0, max(df$clade_2, na.rm=TRUE)),
      ylim = c(0, max(df$clade_2, na.rm=TRUE))
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.justification = c(0, 1), 
      legend.background = element_rect(fill = "transparent"),
      legend.title = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    ggpubr::stat_cor(label.x = 0, label.y=12, method="spearman", r.digits = 3, p.digits = 0.00001, cor.coef.name = "rho")
  
  plot_list[[h]] <- p
}

# Combine with patchwork
P_all <- patchwork::wrap_plots(plot_list)
P_all


# Optional: define custom colors for habitats
hab_colors <- c(
  "Wet heath" = "#9932CC",
  "Decalcified Empetrum dunes" = "#E69F00",
  "Humid dune slacks" = "#E69F00",
  "Degraded raised bog" = "red3",
  "Quaking bogs" = "red3",
  "Active raised bogs" = "red3",
  "Bog woodland" = "#228B22",
  "Hydrophilous tall-herb swamp" = "#0072B2",
  "Molinia meadows" = "#0072B2"
)


#mfd_hab3_filter<-c("Decalcified Empetrum dunes", "Humid dune slacks", "Hydrophilous tall-herb swamp", "Molinia meadows", "Wet heath", "Bog woodland", "Degraded raised bog", "Quaking bogs", "Active raised bogs")


# Plot


library(ggtext)

min_double <- .Machine$double.xmin  # R's smallest double-precision value

format_p <- function(p) {
  ifelse(p < min_double,
         paste0("< ", formatC(min_double, format = "e", digits = 2)),
         formatC(p, format = "e", digits = 2))
}

results_filt<-results%>%filter(species1=="clade_2",species2=="RPKM")%>%select(mfd_hab3, fisher_greater_p_BH, similarity)

cor_stats <- df %>%
  left_join(results_filt,by = join_by(mfd_hab3))%>%
  group_by(mfd_hab3) %>%
  summarise(n_samples = n(),
            rho_val = cor(RPKM, clade_2, method = "spearman", use = "complete.obs"),
            rho_p   = cor.test(RPKM, clade_2, method = "spearman")$p.value,
            lm_fit  = list(lm(clade_2 ~ RPKM)),
            coef_val = coef(lm(clade_2 ~ RPKM))[2],
            coef_p   = summary(lm(clade_2 ~ RPKM))$coefficients[2, 4],
            r2       = summary(lm(clade_2 ~ RPKM))$r.squared,
            J        = first(similarity),  
            fisher_greater_p_BH = first(fisher_greater_p_BH)) %>%
  mutate(rho_p_fmt   = format_p(rho_p),
         coef_p_fmt  = format_p(coef_p),
         J_fmt      = round(J, 2),
         fisher_fmt = if_else(fisher_greater_p_BH==1, paste0("1*"), format_p(fisher_greater_p_BH )))%>%
  mutate(label = paste0(mfd_hab3, " (n = ", n_samples, ")\n",
                        "\u03C1: ", round(rho_val, 2), " (p=", rho_p_fmt, ")\n",
                        "R²: ", round(r2, 2), " (p=", coef_p_fmt, ")\n",
                        "J = ", J_fmt, " (p=", fisher_fmt, ")"
                        # "Reg. coeff.: ", round(coef_val, 3), " (p=", coef_p_fmt, ") ",
         )
  )



legend_labels <- setNames(cor_stats$label, cor_stats$mfd_hab3)


line_data <- df %>%
  group_by(mfd_hab3) %>%
  do({
    model <- lm(clade_2 ~ RPKM, data = .)
    data.frame(
      RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200),
      clade_2 = predict(model, newdata = data.frame(RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200))),
      mfd_hab3 = unique(.$mfd_hab3)
    )
  })


df <- df %>%
  mutate(mfd_hab3 = factor(mfd_hab3, 
                           levels = df %>% 
                             arrange(mfd_hab2) %>% 
                             pull(mfd_hab3) %>% 
                             unique()))

p_all <- ggplot(df, aes(x = RPKM, y = clade_2)) +
  # Transparent regression lines
  geom_line(data = line_data, aes(x = RPKM, y = clade_2, color = mfd_hab3, group = mfd_hab3),
            linetype = 1, linewidth = 0.4, show.legend = F)+
  geom_point(aes(fill = mfd_hab3), shape = 21,
             size = 1.1, alpha = 0.66, stroke = 0.3) +
  scale_color_manual(values = hab_colors, labels = legend_labels) +
  scale_fill_manual(values = hab_colors, labels = legend_labels) +
  labs(
    x = "Methylocystis pmoA2 [RPKM]",
    y = "Methylocystis clade 2 [Rel. Tax. Abund]",
    #  title = "Habitat-wise correlation between pmoA2 and clade_2"
  ) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
  theme_classic() +
  theme(
    legend.spacing.y = unit(1, "cm"),
    legend.spacing.x = unit(0.07, "cm"),
    legend.position = c(0.63, 0.7),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5,margin = margin(t = 2, b = 2, l=2)),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.key.size = unit(0.1, "cm"),
    plot.title = element_text(size = 5, face = "bold", vjust = 1)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 1.5, stroke = 0.4, alpha = 0.9))
  )

p_all

ggsave("correlations/Methylocystis_clade_2_pmoA2_wide_25_11_27.png",
       p_all,
       units = c("mm"),
       height = 70,
       width = 100,
       dpi=300)

ggsave("correlations/Methylocystis_clade_2_pmoA2_wide_25_10_16.svg",
       p_all,
       units = c("mm"),
       height = 70,
       width = 100,
       dpi=300)





results%>%filter(species1=="clade_1", species2=="RPKM")%>%arrange(desc(similarity))



############## The "first" clade: ##################




df1 <- clades_of_cystis_sum %>%
  left_join(points_pmoA2) %>%
  filter(!is.na(RPKM))



species_cols <- grep("clade|RPKM", names(df), value = TRUE)
results_list <- list()
habitats<-df1 %>%
  group_by(mfd_hab3) %>%
  filter(n() >= 5)%>%
  select(mfd_hab3)%>%
  distinct()%>%
  filter(!is.na(mfd_hab3))%>%
  pull(unique(mfd_hab3))

for (hab in habitats) {
  subset_df <- df1 %>% filter(mfd_hab3 == hab)
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
  )%>%relocate(mfd_hab3, fisher_greater_p, fisher_greater_p_BH)


results%>%filter(species1=="clade_1", species2=="RPKM")%>%arrange(desc(similarity))

cor_stats <- df1 %>%
  group_by(mfd_hab3) %>%
  filter(n() >= 5) %>% 
  summarise(
    n_samples = n(),
    rho_val = cor(RPKM, clade_1, method = "spearman", use = "complete.obs"),
    rho_p   = cor.test(RPKM, clade_1, method = "spearman")$p.value,
    lm_fit  = list(lm(clade_1 ~ RPKM)),
    coef_val = if (length(unique(RPKM)) > 1) coef(lm(clade_1 ~ RPKM))[2] else NA_real_,
    coef_p   = if (length(unique(RPKM)) > 1) summary(lm(clade_1 ~ RPKM))$coefficients[2, 4] else NA_real_,
    r2       = if (length(unique(RPKM)) > 1) summary(lm(clade_1 ~ RPKM))$r.squared else NA_real_
  ) %>%
  mutate(
    rho_p_fmt  = format_p(rho_p),
    coef_p_fmt = format_p(coef_p),
    label = paste0(
      mfd_hab3, " (n = ", n_samples, ")\n",
      "R²: ", round(r2, 2), " (p=", coef_p_fmt, ")\n",
      "\u03C1: ", round(rho_val, 2), " (p=", rho_p_fmt, ")"
    )
  )

cor_stats%>%arrange(desc(rho_val))

legend_labels <- setNames(cor_stats$label, cor_stats$mfd_hab3)


line_data <- df %>%
  group_by(mfd_hab3) %>%
  do({
    model <- lm(clade_1 ~ RPKM, data = .)
    data.frame(
      RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200),
      clade_1 = predict(model, newdata = data.frame(RPKM = seq(min(.$RPKM), max(.$RPKM), length.out = 200))),
      mfd_hab3 = unique(.$mfd_hab3)
    )
  })



p_all2 <- ggplot(df, aes(x = RPKM, y = clade_1)) +
  # Transparent regression lines
  geom_line(data = line_data, aes(x = RPKM, y = clade_1, color = mfd_hab3, group = mfd_hab3),
            linetype = 1, linewidth = 0.4, show.legend = F)+
  geom_point(aes(fill = mfd_hab3), shape = 21,
             size = 1.1, alpha = 0.66, stroke = 0.3) +
  scale_color_manual(values = hab_colors, labels = legend_labels) +
  scale_fill_manual(values = hab_colors, labels = legend_labels) +
  labs(
    x = "Methylocystis pmoA2 [RPKM]",
    y = "Methylocystis clade 1 [Rel. Tax. Abund]",
    #  title = "Habitat-wise correlation between pmoA2 and clade_2"
  ) +
  coord_cartesian(xlim = c(0, 7), ylim = c(0, 7)) +
  theme_classic() +
  theme(
    legend.spacing.y = unit(1, "cm"),
    legend.spacing.x = unit(0.07, "cm"),
    legend.position = c(0.63, 0.7),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5,margin = margin(t = 2, b = 2, l=2)),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.key.size = unit(0.1, "cm"),
    plot.title = element_text(size = 5, face = "bold", vjust = 1)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 1.5, stroke = 0.4, alpha = 0.9))
  )

p_all2


