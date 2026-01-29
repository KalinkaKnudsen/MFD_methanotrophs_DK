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




#hab_filter<-c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater")

OTU_filtered_long<-OTU_filtered_long%>%
  filter(mfd_hab1 %in% hab_filter)

p_combine<-readRDS("../palette_mfd_hab2_ISME.rds") 


########### Check if data is normally distributed ###########

Methylocella<-OTU_filtered_long%>%
  filter(Tax_short %in% c("Methylocella_USCa", "Methylocella_pxmA"))%>%
  select(SeqId, Tax_short, RPKM, mfd_hab1, mfd_hab2, mfd_hab3)%>%
  pivot_wider(names_from = Tax_short, values_from = RPKM)%>%
  filter()



#Urban greenspaces (Figure 2c) and agricultural soils (Extended Data Figure 4) showed little methanotrophic marker gene diversity,
#and in most samples only Methylocapsa USCα was observed, albeit inconsistently, independent of crop type, 
#and in lower relative abundance compared to natural upland soils (i.e. forests, grasslands, heaths).


### So, I want to compare the relative abundance of upland habitats with parks and fields. This includes the habitats Semi-natural dry grasslands, Dry heath, and 
#mfd_hab1 =="Forests"   filter(!mfd_hab3 %in% c("Alluvial woodland", "Bog woodland"))%>%

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


#lower total methanotrophy marker genes in greenspaces (Mann–Whitney U-test, U = 1,066,787, two-sided, P = 3.21 × 10−70) and agriculture (Mann–Whitney U-test, U = 1,638,616 , two-sided, P = 6.96 × 10−186)

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


#### This means that I should use Spearman instead of Pearson



my_pub_theme <- theme_classic() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    # Axis titles & text
    axis.title.x = element_text(size = 5, margin = margin(t = 0.1, b=1)),
    axis.title.y = element_text(size = 5, margin = margin(r = 0.05)),
    axis.text = element_text(size = 5),
    
    # Force axis lines to be drawn
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    
    # Plot title & subtitle
    plot.subtitle = element_text(
      size = 5,
      #  hjust = 0.05,              # nudge right (0 = far left, 1 = far right)
      vjust = 1,                 # push down (default is ~0.5–1)
    ),
    
    # Legend text & title
    legend.title = element_text(size = 5, margin = margin(b = 1.2), face="bold"),
    legend.text = element_text(size = 5),
    
    # Legend spacing & key size
    legend.key.size = unit(0.1, "cm"),       # overall key box size
    # legend.key.height = unit(0.2, "cm"),     # vertical spacing
    #  legend.key.width = unit(0.2, "cm"),      # horizontal spacing
    legend.spacing.y = unit(0.07, "cm"),      # space between rows
    legend.spacing.x = unit(0.07, "cm"),      # space between columns
    legend.position = c(0, 1),
    legend.justification = c(0, 1),  
    # Legend background
    legend.background = element_rect(fill = "transparent")
  )

# Smallest positive double R can represent
min_double <- .Machine$double.xmin

# List of habitats you want to test
habitats <- list(
  list(
    name = "Temperate heath and scrub",
    filter_col = "mfd_hab1",
    filter_val = "Temperate heath and scrub"
  ),
  list(
    name = "Bogs, mires and fens",
    filter_col = "mfd_hab1",
    filter_val = "Bogs, mires and fens"
  ),
  list(
    name = "Freshwater",
    filter_col = "mfd_hab1",
    filter_val = "Freshwater"
  ),
  list(
    name = "Dunes",
    filter_col = "mfd_hab1",
    filter_val = "Dunes"
  ),
  list(
    name = "Grassland formations",
    filter_col = "mfd_hab1",
    filter_val = "Grassland formations"
  ),

  list(
    name = "Forests",
    filter_col = "mfd_hab1",
    filter_val = "Forests"
  )
)


results <- lapply(habitats, function(h) {
  df <- Methylocella %>%
    filter(.data[[h$filter_col]] %in% h$filter_val)
  
  lm_fit <- lm(`Methylocella_USCa` ~ `Methylocella_pxmA`, data = df)
  coef_p_raw <- summary(lm_fit)$coefficients[2, 4]
  
  spearman_test <- cor.test(df$Methylocella_pxmA, df$Methylocella_USCa, method = "spearman")
  
  tibble(
    habitat = h$name,
    n = n_distinct(df$SeqId),
    coef_val = round(summary(lm_fit)$coefficients[2, 1], 3),
    coef_p_raw = coef_p_raw,
    r2 = round(summary(lm_fit)$r.squared, 3),
    rho_val = round(spearman_test$estimate, 3),
    rho_p_raw = spearman_test$p.value
  )
}) %>% bind_rows()

# Apply Benjamini–Hochberg correction (or change method)
results <- results %>%
  mutate(
    coef_p_adj = p.adjust(coef_p_raw, method = "BH"),
    rho_p_adj  = p.adjust(rho_p_raw,  method = "BH")
  )

format_p <- function(p) {
  ifelse(p < min_double,
         paste0("< ", formatC(min_double, format = "e", digits = 2)),
         formatC(p, format = "e", digits = 2))
}

results <- results %>%
  mutate(
    coef_p_adj = format_p(coef_p_adj),
    rho_p_adj  = format_p(rho_p_adj)
  )






##############################################



# Initialize empty list
plot_list <- list()

# Loop over each habitat in results
for (i in seq_len(nrow(results))) {
  h <- results$habitat[i]
  
  # Filter data
  df_h <- Methylocella %>% 
    filter(mfd_hab1 == h, !is.na(mfd_hab2))
  
  # Extract stats
  coef_val <- results$coef_val[i]
  coef_p <- results$coef_p_adj[i]
  r2 <- results$r2[i]
  rho_val <- results$rho_val[i]
  rho_p <- results$rho_p_adj[i]
  num_samples <- results$n[i]
  
  # Build plot
  p <- df_h %>%
    ggplot(aes(x = Methylocella_USCa, y = Methylocella_pxmA)) +
    geom_point(aes(fill = mfd_hab2), size = 1, alpha = 0.6,
               color = "black", pch = 21, stroke = 0.3) +
    scale_fill_manual(values = p_combine,
                      name = paste0(h, ", n = ", num_samples)) +
    geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2, lwd = 0.5) +
    labs(
      x = "Methylocella_USCa [RPKM]",
      y = "Methylocella_pxmA [RPKM]"
    ) +
    annotate(
      "text",
      x = 0,
      y = 5.3,
      label = paste0(
        "Reg. coeff.: ", coef_val, " (p=", coef_p, ")\n",
        "R²: ", r2, "\n",
        "\u03C1: ", rho_val, " (p=", rho_p, ")"
      ),
      size = 5 / .pt,
      hjust = 0,
      vjust = 1,
      fontface = "plain"
    ) +
    coord_cartesian(
      xlim = c(0, max(Methylocella$Methylocella_USCa, na.rm = TRUE)),
      ylim = c(0, max(Methylocella$Methylocella_pxmA, na.rm = TRUE)+0.6)
    ) +
    my_pub_theme +
    guides(fill = guide_legend(ncol = 2))
  
  plot_list[[h]] <- p
}

# Combine with patchwork
P_all <- wrap_plots(plot_list)


#P_all

#
ggsave("correlations/Correlations_Methylocella_25_10_06.svg",
       P_all,
       units = c("mm"),
       height = 110,
       width = 183,
       dpi=300)







#




## But that is based on the relative abundance, which is somewhat difficult for several reasons. I can, however, try to do the Jaccard instead


### Lets make it clean ##

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



############ making some tests #######################

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
## Final attempt is figuring out if I do this across all genes, then what happens:
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













###########################################################
########## Moving on to Methylomonadaceae #################

Methylomonadaceae<-OTU_filtered_long%>%
  filter(Tax_short %in% c("Methylomonadaceae_pxmA", "Methylomonadaceae"))%>%
  select(SeqId, Tax_short, RPKM, mfd_hab1, mfd_hab2, mfd_hab3)%>%
  pivot_wider(names_from = Tax_short, values_from = RPKM)%>%
  filter()




##### 
## Test for normality
Methylomonadaceae %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Methylomonadaceae_pxmA)

Methylomonadaceae %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Methylomonadaceae)


qqnorm(Methylomonadaceae$Methylomonadaceae)
qqline(Methylomonadaceae$Methylomonadaceae)


#### Grass
Methylomonadaceae_grass<-Methylomonadaceae%>%
  filter(mfd_hab1=="Grassland formations")


cor(Methylomonadaceae_grass$Methylomonadaceae, Methylomonadaceae_grass$Methylomonadaceae_pxmA, method = c("spearman"))
cor.test(Methylomonadaceae_grass$Methylomonadaceae, Methylomonadaceae_grass$Methylomonadaceae_pxmA, method=c("spearman"))


P_grass <-Methylomonadaceae_grass %>%
  filter(!is.na(mfd_hab3))%>%
  ggplot(., aes(
    y = Methylomonadaceae,
    x= Methylomonadaceae_pxmA)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine)+ #name = paste0("Bogs, Mires and Fens,\nn = ",num_samples)
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  #  theme_classic() +
  labs(
    y = "Methylomonadaceae_pmoA [RPKM]",
    x = "Methylomonadaceae_pxmA [RPKM]",
  )+
  coord_cartesian(xlim = c(0, max(Methylomonadaceae_grass$`Methylomonadaceae_pxmA`)), ylim = c(0,max(Methylomonadaceae_grass$Methylomonadaceae)))+
  theme(legend.position = c(0.84,0.15),
        legend.background = element_rect(fill = "transparent"))+
  facet_wrap(.~mfd_hab3)+
  guides(fill = guide_legend(ncol = 2))+
  ggpubr::stat_cor(label.x = 1, method=c("spearman"), r.digits = 3,  cor.coef.name = c("rho"))

P_grass



#### Sediment
Methylomonadaceae_sedi<-Methylomonadaceae%>%
  filter(mfd_hab1=="Freshwater")%>%  filter(Methylomonadaceae<20)



P_sedi <-Methylomonadaceae_sedi %>%
  #filter(!is.na(mfd_hab3))%>%
  ggplot(., aes(
    y = Methylomonadaceae,
    x= Methylomonadaceae_pxmA)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine)+ #name = paste0("Bogs, Mires and Fens,\nn = ",num_samples)
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  #  theme_classic() +
  labs(
    y = "Methylomonadaceae_pmoA [RPKM]",
    x = "Methylomonadaceae_pxmA [RPKM]",
  )+
  coord_cartesian(xlim = c(0, max(Methylomonadaceae_sedi$`Methylomonadaceae_pxmA`)), ylim = c(0,max(Methylomonadaceae_sedi$Methylomonadaceae)))+
  theme(legend.position = c(0.84,0.15),
        legend.background = element_rect(fill = "transparent"))+
  facet_wrap(.~mfd_hab2)+
  guides(fill = guide_legend(ncol = 2))+
  ggpubr::stat_cor(label.x = 2, method=c("spearman"), r.digits = 3,  cor.coef.name = c("rho"))

P_sedi












###########################################################
########## Moving on to Rhodomicrobium #################

#The proposed novel group of methanotrophs within Rhodomicrobium appeared to mimic that of the Methylocystis/Methylosinus clade 
#where pmoA and mmoX co-occur in some bogs and fens, dunes, alluvial- and bog woodlands, across nearly all freshwater sediments, 
#some humid meadows, and a few temperate heath samples, albeit in lower relative abundance compared to canonical methanotrophs (Figure 2). 

Rhodo_cystis<-OTU_filtered_long%>%
  filter(Tax_curated %in% c("Beijerinckiaceae; Methylosinus_Methylocystis_pmoA1", "Rhizobiales; Rhodomicrobium", "Rhodomicrobium_pmoA"))%>%
  select(SeqId, Tax_short, RPKM, mfd_hab1, mfd_hab2, mfd_hab3)%>%
  pivot_wider(names_from = Tax_short, values_from = RPKM)%>%
  rename(Rhodomicrobium_mmoX=Rhodomicrobium)%>%
  filter()




##### 
## Test for normality
Rhodo_cystis %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Rhodomicrobium_mmoX)

Rhodo_cystis %>% ungroup()%>%
  group_by(mfd_hab1) %>%
  rstatix::shapiro_test(Rhodomicrobium_pmoA)
## absolutely not normal, p much lower than 0.05



### 1, show that Rhodo mmoX and pmoA correlates well (in most habitats)


cor.test(Rhodo_cystis$Rhodomicrobium_pmoA, Rhodo_cystis$Rhodomicrobium_mmoX, method=c("spearman"))
#not super great

P_all <-Rhodo_cystis %>%
  filter(!is.na(mfd_hab2))%>%
  ggplot(., aes(
    y = Rhodomicrobium_pmoA,
    x= Rhodomicrobium_mmoX)) +
  geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6, color = "black", pch = 21, stroke=.8)+
  scale_fill_manual(values=p_combine)+ #name = paste0("Bogs, Mires and Fens,\nn = ",num_samples)
  geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
  #  theme_classic() +
  labs(
    y = "Rhodomicrobium_pmoA [RPKM]",
    x = "Rhodomicrobium_mmoX [RPKM]",
  )+
  coord_cartesian(xlim = c(0, max(Rhodo_cystis$Rhodomicrobium_pmoA)), ylim = c(0,max(Rhodo_cystis$Rhodomicrobium_pmoA)))+
  theme(legend.position = c(0.84,0.15),
        legend.background = element_rect(fill = "transparent"))+
  facet_wrap(.~mfd_hab1)+
  guides(fill = guide_legend(ncol = 2))+
  ggpubr::stat_cor(method=c("spearman"), r.digits = 3,  cor.coef.name = c("rho"))

P_all



# Get unique mfd_hab1 groups
habitats <- Rhodo_cystis%>%filter(!mfd_hab1=="Greenspaces")%>%pull(mfd_hab1)%>%unique()

habitats2<-Rhodo_cystis%>%filter(
  mfd_hab1 %in% c("Freshwater") |
  mfd_hab2 %in% c("Calcareous fens") |
  mfd_hab3 %in% c("Humid dune slacks", "Hydrophilous tall-herb swamp", "Molinia meadows", "Alluvial woodland", "Bog woodland", "Wet heath"))%>%pull(mfd_hab1)%>%unique()

# Initialize empty list to store plots
plot_list <- list()

# Loop over habitats
for (h in habitats) {
  p <- Rhodo_cystis %>%
    #filter(!(Rhodomicrobium_pmoA == 0 & Rhodomicrobium_mmoX == 0))  %>% Removing the double-zeroes here influences the correlation like crazy!
    filter(!is.na(mfd_hab2), mfd_hab1 == h) %>%
    ggplot(aes(
      y = Rhodomicrobium_pmoA,
      x = Rhodomicrobium_mmoX)) +
    geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6,
               color = "black", pch = 21, stroke=.8) +
    scale_fill_manual(values = p_combine) +
    geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
    labs(
      y = "Rhodomicrobium_pmoA [RPKM]",
      x = "Rhodomicrobium_mmoX [RPKM]",
      title = h
    ) +
    coord_cartesian(
      xlim = c(0, max(Rhodo_cystis$Rhodomicrobium_mmoX, na.rm=TRUE)),
      ylim = c(0, max(Rhodo_cystis$Rhodomicrobium_mmoX, na.rm=TRUE))
    ) +
    theme_classic() +
    theme(
      legend.position = c(0.6,0.8),
      legend.background = element_rect(fill = "transparent", color="black"),
     # legend.title = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 2)) +
    ggpubr::stat_cor(label.x = 0, label.y=1.8, method="spearman", r.digits = 3, cor.coef.name = "rho")
  
  plot_list[[h]] <- p
}

# Combine with patchwork
P_all <- wrap_plots(plot_list)
P_all


ggsave("correlations/Correlations_Rhodomicrobium.png",
       P_all,
       height = 10,
       width = 18,
       dpi=400)

#


##########################
######### comparing to cystis

# Get unique mfd_hab1 groups
habitats <- Rhodo_cystis%>%filter(!mfd_hab1=="Greenspaces")%>%pull(mfd_hab1)%>%unique()

# Initialize empty list to store plots
plot_list <- list()

# Loop over habitats
for (h in habitats) {
  p <- Rhodo_cystis %>%
    filter(!is.na(mfd_hab2), mfd_hab1 == h) %>%
    ggplot(aes(
      y = Methylosinus_Methylocystis_pmoA1,
      x = Rhodomicrobium_pmoA)) +
    geom_point(aes(fill = mfd_hab2), size = 2, alpha = 0.6,
               color = "black", pch = 21, stroke=.8) +
    scale_fill_manual(values = p_combine) +
    geom_smooth(method = "lm", se = FALSE, col = "darkred", lty = 2) +
    labs(
      y = "Methylosinus_Methylocystis_pmoA1 [RPKM]",
      x = "Rhodomicrobium_pmoA [RPKM]",
      title = h
    ) +
    coord_cartesian(
      xlim = c(0, 10),
      ylim = c(0, 10)
    ) +
    theme_classic() +
    theme(
      legend.position = c(0.6,0.8),
      legend.background = element_rect(fill = "transparent", color="black"),
      # legend.title = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 2)) +
    ggpubr::stat_cor(label.x = 0, label.y=10, method="spearman", r.digits = 3, cor.coef.name = "rho")
  
  plot_list[[h]] <- p
}

# Combine with patchwork
P_all <- wrap_plots(plot_list)
P_all


ggsave("correlations/Correlations_Rhodomicrobium_pmoA_vs_Cystis_pmoA.png",
       P_all,
       height = 10,
       width = 18,
       dpi=400)




## Trying the jaccard 



df<-Rhodo_cystis%>%filter(
  mfd_hab1 %in% c("Freshwater") |
    mfd_hab2 %in% c("Calcareous fens") |
    mfd_hab3 %in% c("Humid dune slacks", "Hydrophilous tall-herb swamp", "Molinia meadows", "Alluvial woodland", "Bog woodland", "Wet heath"))

species_cols <- grep("Seq*|mfd*", names(df), value = TRUE, invert=T) #OBS, inverted!
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

co_occur<-results%>%filter(fisher_greater_p_BH<0.05)%>%arrange(desc(similarity))%>%filter(grepl("Methylocystis*", species2))
