# Script for analyzing gene correlations
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



OTU_filtered_long<-OTU_filtered_long%>%
  filter(mfd_hab1 %in% hab_filter)

p_combine<-readRDS("data/palette_mfd_hab2_ISME.rds") 



Methylocella<-OTU_filtered_long%>%
  filter(Tax_short %in% c("Methylocella_USCa", "Methylocella_pxmA"))%>%
  select(SeqId, Tax_short, RPKM, mfd_hab1, mfd_hab2, mfd_hab3)%>%
  pivot_wider(names_from = Tax_short, values_from = RPKM)%>%
  filter()



Tax_filter<-c("TUSC", "Likely_mmoX", "AVCC01_Methylotenera_clade", "Methylomonadaceae; Methyloprofundus_WTBX01_clade",  #"USCg_Methyloglobulus",
              "Methylococcales_unknown", "o_Rhizobiales_pmoA", "Methylococcaceae; Methylogaea", "Methylomagnum", 
              "Rhizobiales; Beijerinckiaceae; Methyloferula_Methylovirgula", "Methylomonadaceae; SXIZ01", "Methylomonadaceae; Methyloprofundus", "pxmA_umbrella", "Methylohalobius", "Methylococcales_mmoX_umbrella", "Beijerinckiaceae_pmoA1", "Methylomonadaceae; Methylomarinum")

comparison_fields_to_rest<-OTU_filtered_long%>%
   filter(mfd_hab2 %in% c("Semi-natural dry grasslands") | 	mfd_hab3 %in% c("Dry heath") | mfd_hab1 %in% c("Forests", "Greenspaces")| mfd_areatype =="Agriculture")%>%
   filter(!mfd_hab3 %in% c("Alluvial woodland", "Bog woodland"))%>%
   filter(!type %in% c("Put. pmoA/mmoX"))%>%
   filter(!Tax_curated%in%Tax_filter)%>%
  group_by(SeqId, mfd_areatype,mfd_hab1, mfd_hab2, mfd_hab3)%>%
  summarise(gene_tot=sum(RPKM))%>%
  group_by(mfd_hab1)%>%
  mutate(genes_sum=sum(gene_tot),
         genes_mean=mean(gene_tot),
         genes_median=median(gene_tot),
         genes_max=max(gene_tot))%>%
  ungroup()

dat<-rstatix::wilcox_test(comparison_fields_to_rest,
                          gene_tot~mfd_areatype,
                          p.adjust.method = "BH", detailed=T)



############ Only Methylocella USCa ##########
comparison_fields_to_rest_USCa<-OTU_filtered_long%>%
  filter(mfd_hab2 %in% c("Semi-natural dry grasslands") | 	mfd_hab3 %in% c("Dry heath") | mfd_hab1 %in% c("Forests", "Greenspaces")| mfd_areatype =="Agriculture")%>%
  filter(!mfd_hab3 %in% c("Alluvial woodland", "Bog woodland"))%>%
  filter(!type %in% c("Put. pmoA/mmoX"))%>%
  filter(Tax_short %in%c("Methylocella_USCa", "Methylocella_pxmA"))%>%
  group_by(SeqId, mfd_areatype,mfd_hab1, mfd_hab2, mfd_hab3)%>%
  summarise(gene_tot=sum(RPKM))%>%
  group_by(mfd_hab1)%>%
  mutate(genes_sum=sum(gene_tot),
         genes_mean=mean(gene_tot),
         genes_median=median(gene_tot),
         genes_max=max(gene_tot))%>%
  ungroup()

dat_USCa<-rstatix::wilcox_test(comparison_fields_to_rest_USCa,
                               gene_tot~mfd_areatype,
  p.adjust.method = "BH", detailed=T)

# greenspaces (Mann–Whitney U-test, U = 1117357, two-sided, P = 8.47 × 10−98) and agriculture (Mann–Whitney U-test, U = 1651720 , two-sided, P = 9.69 × 10−197)


##### 
## Test for normality
Methylocella %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Methylocella_USCa)

Methylocella %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Methylocella_pxmA)


qqnorm(Methylocella$Methylocella_pxmA)
qqline(Methylocella$Methylocella_pxmA)


######### Not at all normally distributed (zero-inflated)

### Jaccard similarity with Fisher's exact test function

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


df<-Methylocella
species_cols <- grep("Methylocella_*", names(df), value = TRUE)
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
  )%>%relocate(mfd_hab1, fisher_greater_p, fisher_greater_p_BH)

results%>%arrange(similarity)


### On hab2 instead

df<-Methylocella
species_cols <- grep("Methylocella_*", names(df), value = TRUE)
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

results%>%arrange(similarity)



############ Tesing on the data #######################

Methylocella2<-Methylocella%>%filter(mfd_hab1 %in% c("Grassland formations", "Temperate heath and scrub", "Forests", "Bogs, mires and fens"))

plot_test <- Methylocella2%>%
  #filter(mfd_hab1 %in% c("Fields"))%>%
  #mutate(label = factor(label, levels = unique(fields$label), ordered = TRUE)) %>%
  ggplot(., aes(x = mfd_hab3, y = Methylocella_pxmA)) + 
  geom_jitter(aes(fill=mfd_hab2), size = 1.3, alpha = 0.8, color = "black", pch = 21,
              position = position_jitter(width = 0.2, height = 0.03) # Adjust jitter spread if needed
  ) +
  geom_violin(aes(fill=mfd_hab2),  # Add median quantile line
              color = "black",     # Median line color
              size = 0.6,
              alpha=0.4
              #fill="transparent"# Median line thickness
  ) +
  scale_fill_manual(values = p_combine)+
  stat_summary(fun = mean, geom = "crossbar", width = 0.8, size = 0.6, color = "red3") +
  # stat_compare_means(method = "wilcox.test", aes(label=..p.adj..),
  #                    ref.group = ".all.", method.args = list(alternative = "less"))+ ###Adjust here for one sided "less" I think
  ggpubr::geom_pwc(method = "wilcox_test",
                   label = 'p.adj.signif',
                   ref.group = "all",
                   hide.ns = T,
                   #angle = 90, vjust = 0.75,
                   p.adjust.method = "fdr",
                   group.by = "x.var",
                   p.adjust.by = "panel",
                   method.args = list(alternative = "greater"),
                   # method.args = list(alternative = "greater"),
                   remove.bracket = T)+
  
  geom_hline(yintercept = mean(Methylocella2%>%pull(Methylocella_pxmA)), linetype = 2)+ # Add horizontal line at base mean
  scale_y_sqrt() +
  facet_grid(.~mfd_hab1, scales="free_x", space = "free")+
  theme_minimal()+
  theme(legend.position="bottom", panel.grid=element_blank(), axis.text.x = element_text(angle=45, hjust = 1))

#label = "p.signif", #aes(label=..p.adj..),
#p.adjust.methods = "fdr",

stats<-thinned_10k%>%
  ungroup()%>%
  filter(mfd_hab1 %in% c("Fields"))%>%
  #group_by(mfd_hab3)%>%
  rstatix::wilcox_test(drep_abundance~label,alternative = "greater", p.adjust.method = "fdr", ref.group = "all", detailed=T)%>%print(n = 24)




##########################################################
######################### USCg ###########################
##########################################################


USCg<-OTU_filtered_long%>%
  filter(grepl("USCg_umbrella|USCg; JACCX01|USCg; USCg_Taylor|USCg|USCa", Tax_curated))%>%
  mutate(Tax_curated=if_else(grepl("USCg_umbrella|USCg; JACCX01|USCg; USCg_Taylor", Tax_curated), "USCg", Tax_curated))%>%
#  filter(grepl("USCg_umbrella|USCg; JACCX01|USCg; USCg_Taylor", Tax_curated))%>%
  group_by(SeqId, Tax_curated) %>%
  summarise(
    abundance_sum = sum(RPKM, na.rm = TRUE),          # or use mean(RPKM) if you prefer
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3),
    SeqId     = first(SeqId)
  ) %>%
  ungroup%>%
  pivot_wider(names_from = Tax_curated, values_from = abundance_sum)




df<-USCg
species_cols <- grep("USC*", names(df), value = TRUE)
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

results%>%filter(!grepl("USCg_Methyloglobulus", species2))%>%arrange(fisher_less_p_BH)

#################################################
######### Going to specific grasslands #########


df<-USCg%>%filter(mfd_hab3 %in% c("Calcareous grassland", "Xeric sand calcareous grasslands"))
species_cols <- grep("USC*", names(df), value = TRUE)
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

results%>%filter(!grepl("USCg_Methyloglobulus", species2))%>%arrange(fisher_less_p_BH)




res<-jaccard_with_fisher(df$`Beijerinckiaceae; Methylocella; Methylocella_USCa`, df$USCg)

#########################
## Across all genes for investigating 
######################

USCg_v2<-OTU_filtered_long%>%
  filter(!grepl("TUSC", Tax_curated))%>%
  filter(!type=="Put. pmoA/mmoX")%>%
  mutate(Tax_curated=if_else(grepl("USCg_umbrella|USCg; JACCX01|USCg; USCg_Taylor", Tax_curated), "USCg", Tax_curated))%>%
  #  filter(grepl("USCg_umbrella|USCg; JACCX01|USCg; USCg_Taylor", Tax_curated))%>%
  group_by(SeqId, Tax_curated) %>%
  summarise(
    abundance_sum = sum(RPKM, na.rm = TRUE),          # or use mean(RPKM) if you prefer
    mfd_hab1  = first(mfd_hab1),
    mfd_hab2  = first(mfd_hab2),
    mfd_hab3  = first(mfd_hab3),
    SeqId     = first(SeqId)
  ) %>%
  ungroup%>%
  pivot_wider(names_from = Tax_curated, values_from = abundance_sum)




df<-USCg_v2
species_cols <- grep("Seq*|mfd*|pxmA_umbrella|TUSC", names(df), value = TRUE, invert=T)
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

exclusive_USCg<-results%>%filter(fisher_less_p_BH<0.05)%>%filter(similarity<0.5)
co_occur_USCg<-results%>%filter(fisher_greater_p<0.05)%>%filter(similarity>0.5)

