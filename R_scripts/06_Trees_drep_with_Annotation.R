
#!/usr/bin/env Rscript
## ============================================================================
## 06_Trees_drep_with_Annotation.R
## Purpose: Create phylogenetic trees of dereplicated methanotroph MAGs with
##          metabolic annotations and genome quantification (Sylph) overlays.
##          Produces multi-panel phylogenetic visualizations for major lineages
##          showing bootstrap support, metabolic gene presence, and habitat abundance.
## ============================================================================
setwd("path/to/your/repo/MFD_methanotrophs_DK/")

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(stringr)
library(ggh4x)


# ===== PHYLOGENETIC TREE LOADING AND PROCESSING =====
# Load phylogenetic tree in Newick format from dereplicated genome phylogeny
tree_drep <- treeio::read.tree("data/phylogenomic_tree.treefile")

# Remove GTDB reference prefixes (RS_, GB_) from tip labels to leave only genome IDs
tree_drep$tip.label <- gsub("RS_", "", tree_drep$tip.label)
tree_drep$tip.label <- gsub("GB_", "", tree_drep$tip.label)

# ===== TAXONOMY DATA LOADING AND CORRECTIONS =====
# Load taxonomy classification and genome metadata
# Classify genomes by origin (Microflora Danica vs GTDB references)
# Correct misclassified Methylocella species names to proper Methylocapsa genus assignments
tax<-readRDS("data/MFD_renamed_tax_25_03_04.rds")%>%
  mutate(type = if_else(grepl("MFD", user_genome), 'Microflora Danica','GTDB'),
         type=if_else(grepl("LIB", user_genome), "Microflora Danica", type))%>%
  mutate(tip.label=as.character(user_genome))%>%
  relocate(tip.label)%>%
  # Apply species name corrections for three Methylocella genomes reclassified as Methylocapsa
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))%>%
  # Apply same corrections to display labels
  mutate(label_3=gsub("Methylocella sp002890675", "Ca. Methyloaffinis lahnbergensis", label_3))%>%
  mutate(label_3=gsub("Methylocella sp004564215", "Methylocapsa gorgona", label_3))%>%
  mutate(label_3=gsub("Methylocella sp029855125", "Methylocapsa sp. D3K7", label_3))

# ===== DUPLICATE TAXA REMOVAL =====
# Identify duplicate GTDB reference species to avoid redundancy
duplicates <- tax %>%
  group_by(Species) %>%
  filter(n() > 1) %>%
  ungroup()%>%
  filter(type=="GTDB")%>%
  pull(tip.label)

# Remove duplicate GTDB taxa from taxonomy table
tax<-tax%>%filter(!tip.label %in% duplicates)

# ===== TREE ROOTING AND DUPLICATE TIP REMOVAL =====
# Perform midpoint rooting on the phylogenetic tree to place the tree root at the point equidistant from all tips
tr <- phytools::midpoint.root(tree_drep)
# Remove tips corresponding to duplicate GTDB taxa from the tree
pr<-ape::drop.tip(tr, tip = duplicates)
pr_tree <- ggtree::ggtree(pr) 
pr_tree


# ===== SYLPH GENOME QUANTIFICATION LOADING AND HABITAT STANDARDIZATION =====
# Load Sylph quantification data and apply species name corrections matching taxonomy
sylph<-readRDS("data/sylph_gtdb_25_03_06.rds")%>%
  mutate(Species=gsub("s__Methylocella sp002890675", "s__Ca. Methyloaffinis lahnbergensis", Species))%>%
  mutate(Species=gsub("s__Methylocella sp004564215", "s__Methylocapsa gorgona", Species))%>%
  mutate(Species=gsub("s__Methylocella sp029855125", "s__Methylocapsa sp. D3K7", Species))

# Standardize habitat metadata through series of string replacements and conditional mutations:

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
  # Create formatted habitat labels for visualization with line breaks
  mutate(hab1_label=gsub("Bogs, mires and fens", "Bogs, mires\nand fens", mfd_hab1),
         hab1_label=gsub("Grassland formations", "Grassland\nformations", hab1_label),
         hab1_label=gsub("Temperate heath and scrub","Heath and\nscrub", hab1_label),
         hab1_label=gsub("Greenspaces","Green-\nspaces", hab1_label),
         hab1_label=gsub("Coastal","Coast", hab1_label),
         hab1_label=gsub("Freshwater", "Freshwater\nsediment", hab1_label))%>%
  filter(!mfd_hab2 %in% c("Spruce", "Willow"))

# Load pre-computed habitat ordering and color palette for consistent visualization
levels_hab2<-readRDS("output/hab2_sort_order.rds")
palette_mfd_hab2<-readRDS("data/palette_mfd_hab2_ISME.rds")

# Sort samples by habitat for heatmap ordering
hab2_sort<-sylph%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  arrange(SeqId)

# ===== KEGG ANNOTATION DATA LOADING AND PATHWAY LABEL ABBREVIATION =====
# Load KEGG/custom HMM gene annotation table with metabolic classifications
KEGG<-readRDS("output/KEGG_25_12_02.rds")%>%
  # Create tree label combining taxonomic ranks with genome ID for tree tip labeling
  mutate(tree_label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))

# Load color palette for methane oxidation and metabolic enzyme visualization
methane_colors<-readRDS("data/palette_methane_colors.rds")%>%
  unlist()
names(methane_colors)[names(methane_colors) == "-"] <- "Formate oxidation"

# Load KO metadata with gene names and pathway information
# Create abbreviated pathway names and formatted metabolic step labels for compact visualization
KOs<-read_excel("data/KO_KSK_methanotroph_paper.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  # Create gene labels combining gene name with metabolic step
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  # Format metabolic step labels with line breaks for faceted visualization
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

# Normalize Unicode spaces in KO field
KOs$KO <- gsub("\u00A0", " ", KOs$KO)

# Sort KEGG data by genus for consistent ordering
KEGG_sort<-KEGG%>%
  arrange(Genus)

# Load habitat sorting order
hab2_sort_order<-readRDS("output/hab2_sort_order.rds")

library(tidyverse)

# ===== RHODOMICROBIUM LINEAGE SECTION =====
# From here, we create phylogenetic trees and metabolic visualizations for the Rhodomicrobium genus

# Filter taxonomy to Rhodomicrobium genus and extract tip labels for tree subsetting
outgroup_tips <- tax %>%
  filter(Genus=="g__Rhodomicrobium")%>%
  pull(tip.label)

# Extract phylogenetic tree for Rhodomicrobium genomes
Rhodo_tree <- ape::keep.tip(pr, tip = outgroup_tips)

# Create phylogenetic tree plot with tip labels colored by genome origin (GTDB vs Microflora Danica)
Rhodo_plot<-ggtree(Rhodo_tree, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"),guide = "none")+
  theme(legend.position = "none") +
  scale_x_discrete(expand = c(0,0))

# ===== BOOTSTRAP SUPPORT VISUALIZATION =====
# Extract bootstrap confidence values from internal tree nodes (converted from tree labels)
# Assign confidence categories: 95-100%, 85-95%, 0-85%
df1 <- Rhodo_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]

# Add bootstrap confidence node points to phylogenetic tree, color coded by confidence level
Rhodo_plot <- Rhodo_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5)) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),
    keywidth = unit(0.1, "cm"),
    override.aes = list(size = 1.2)
  )) 

# Add outgroup indicator arrow and tree scale
ax<-Rhodo_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Rhodo_plot$data%>%filter(node==min(parent))%>%pull(y)

Rhodo_plot <- Rhodo_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3,
               inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)

Rhodo_plot 

# ===== GENOME ABUNDANCE HEATMAP (SYLPH QUANTIFICATION) =====
# Create heatmap of Sylph quantification (Taxonomic_abundance) for Rhodomicrobium genomes across habitat types
# Filter to specific habitats (Bogs, Grassland, Heath, Dunes, Forests, Freshwater) excluding mire and water samples
heat <- sylph %>%
  filter(user_genome %in% Rhodo_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>%
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # Color gradient scaling: white (absent/low) to black (abundant), with sqrt transformation for better visualization
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),
    trans = "sqrt",
    limits = c(0.00, 2),
    na.value = "black"
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

heat

# ===== METABOLIC GENE HEATMAP =====
# Create heatmap of metabolic genes (RuBisCo and Rhodomicrobium-specific mmoX variants) for Rhodomicrobium genomes
# Filter to genes with presence=1 and relabel Rhodomicrobium mmoX variants for clarity
pl_simple <- KEGG %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(grepl("rbcL - RuBisCo|rbcS - RuBisCo|Rhodomicrobium", Gene))%>%
  mutate(gene_label=if_else(grepl("Rhodomicrobium_umbrella|Rhodomicrobium_like", gene_label), "Rhodomicrobium_mmoX - (1/3)", gene_label))%>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Metabolic_step=gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step))%>%
  # Reorder gene labels according to KO order for consistent visualization
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  ggplot(., aes(x = gene_label, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(methane_colors)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

pl_simple

# ===== PLOT COMPOSITION: TREE + METABOLIC HEATMAP + ABUNDANCE HEATMAP =====
# Compose multi-panel figure combining phylogenetic tree, metabolic gene heatmap, and genome abundance
# Use aplot::insert_left/insert_right to position panels side-by-side with width ratios
tree_met <-pl_simple  %>% aplot::insert_left(Rhodo_plot, width=2)
tree_met_heat<-tree_met %>% aplot::insert_right(heat, width=3.5)

# ===== CONDENSED METABOLIC HEATMAP (SECONDARY VISUALIZATION) =====
# Create condensed version of metabolic gene heatmap by summarizing gene presence across modules
# Filter to major metabolic pathways: CO/H2 oxidation, nitrogen fixation, sulfur metabolism, Complex IV, PHB, denitrification
KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  # Summarize gene presence across each Type (methane oxidation type) for each genome and metabolic step
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),
    genes_present = n_distinct(Gene_collapsed[presence == 1]),
    proportion_present = genes_present / total_genes
  ) %>%
  # Set presence to NA if proportion of genes present is below threshold (60%)
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  # Abbreviate metabolic step labels for compact visualization
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "CO and\nH2", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()

# Create condensed metabolic heatmap for Rhodomicrobium, filtered to major pathways
pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Rhodo_tree$tip.label)%>%
  filter(!Module=="Oxi. Phos. Alternative")%>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "GSH"="#5F9EA0",
    "H4F"="#C7966D",
    "H4MPT"="green",
    "Carbon monoxide oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "Assimilatory sulphate reduction"="#0DB7A8",
    "sat"="#6C961E",
    "aprAB" = "#18B862",
    "DsrABL" = "#77AA79",
    "DsrC"= "#40B645",
    "dsrMKJOP"= "lightblue4",
    "dsrEFH"="#5B815B",
    "SOX"="#92C76C",
    "TCA cycle"="#C7BB6D",
    "PHB"="#BC8615",
    "acetate uptake"="#C7966D",
    "EMC"="#BD5531",
    "PPP"="#AD8ACC",
    "Glycolysis"="#1237B7",
    "cyt. c aa3-type oxidase"="#d97512",
    "High affin. cyt.bd-type"="turquoise4",
    "Cyt. bo3 ubiquinol oxi."="darkred",
    "High affin. cyt c bb3-type"="darkblue"))+
  facet_nested(. ~factor(Metabolic_step, levels=c("From.\noxi","Complex\nIV", "CO and\nH2", "N", "S\noxi./ red.", "S\noxi.","PHB")), scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

options("aplot_guides" = "keep")

tree_met <-pl_simple  %>% aplot::insert_left(Rhodo_plot, width=5.8)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=5.5)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=8.6)


ggsave("./output/Rhodo_met_condense_Sylph_drep_25_11_10.png",tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 186,
       dpi=300)



ggsave("./output/Rhodo_met_condense_Sylph_drep_25_11_10.svg",tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 186,
       dpi=300)



bar_plot <- sylph %>%
  filter(user_genome %in% Rhodo_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_nested(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs, mires\nand fens", "Freshwater\nsediment","Dunes","Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

library(patchwork)
combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))


ggsave("./output/Rhodo_Sylph_bar_25_11_10.svg",combined_plot,
       units = c("mm"),
       height = 110,
       width = 86,
       dpi=300)


#################################################################################################
#################################################################################################
######################      From here, we are working with Methylocella      ####################
#################################################################################################
#################################################################################################



# ===== METHYLOCELLA LINEAGE SECTION =====
# Filter taxonomy to Methylocella genus and extract tip labels

outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocella")%>%
  filter(!label_3=="Methylocella MAG_10")%>%  ### two members of the same species here, quantified as 1 species in sylph
  pull(tip.label)


tax2<-tax%>%
  mutate(clade = if_else(between(genome_number, 1, 23), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 24, 30), "clade_2", clade))%>
  mutate(label_3=gsub("Methylocella", "Methylocapsa", label_3))%>%
  mutate(label_3=gsub("Ca. Methyloaffinis", "Ca. M.", label_3))%>%
  mutate(label_3=gsub("Methylocapsa sp003162995", "M. sp003162995", label_3))

tax2_print_phylo<-tax2%>%
  filter(Genus=="g__Methylocella")%>%
  select(tip.label, user_genome, label_3, clade)


# remove outgroup,I use this instead
Methylocella_tree <- ape::keep.tip(tr, tip = outgroup_tips)
Methylocella_plot<-ggtree(Methylocella_tree, linewidth=0.3) %<+% tax2 + 
  geom_tiplab(aes(label = label_3, fill = clade, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("clade_1" = "#B3A2C7", "clade_2"="#8EB4E3", "other"="grey100"),guide = "none")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0,0))

Methylocella_plot



# put on bootstrap values
df1 <- Methylocella_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


Methylocella_plot <- Methylocella_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5)) +
 guides(color = guide_legend(
       title = "Bootstrap:",
       title.theme = element_text(size=5),
       label.theme = element_text(size=5),
       keyheight = unit(0.1, "cm"),  # Adjust as needed
       keywidth = unit(0.1, "cm"),   # Adjust as needed
       override.aes = list(size = 1.2)
          )) 

ax<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocella_plot$data%>%filter(node==min(parent))%>%pull(y)

Methylocella_plot <- Methylocella_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3,
               inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)


Methylocella_plot 


options("aplot_guides" = "keep")

heat <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

bar_plot <- sylph %>%
  filter(user_genome %in% Methylocella_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_nested(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs, mires\nand fens", "Freshwater\nsediment","Dunes","Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

pl_simple <- KEGG %>%
  filter(genome %in% Methylocella_tree$tip.label) %>%
  #filter(presence==1) %>%
  mutate(presence = if_else(presence == 0, NA, presence)) %>%
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  #mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step)) %>%
  mutate(Gene = if_else(KO %in% c("Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1", "Root; o_Rhizobiales_pmoA; Beijerinckiaceae_pmoA1; Methylocella"), "Beijerinckiaceae_pmoA", Gene)) %>%
  mutate(Gene = if_else(KO %in% c("Root; Likely_mmoX; Rhizobiales_mmoX", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_2",
                                  "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methylocella_1", "Root; Likely_mmoX; Rhizobiales_mmoX; Beijerinckiaceae; Methyloferula_Methylovirgula"), "Beijerinckiaceae mmoX", Gene)) %>%
  filter(Gene %in% c("Beijerinckiaceae_pmoA", "Beijerinckiaceae mmoX", "rbcL - RuBisCo", "rbcS - RuBisCo", "Methylocella_pxmA", "Methylocella_USCa", "Beijerinckiaceae_pxmA", "Homologous_Methylocella", "pmoA_amoA_pxmA",
                     "mxaF", "mxaI", "xoxF")) %>%
  mutate(Gene = factor(Gene, levels = c("Beijerinckiaceae_pmoA", "Beijerinckiaceae_pxmA", "Methylocella_USCa", "Methylocella_pxmA", "pmoA_amoA_pxmA", "Beijerinckiaceae mmoX","Homologous_Methylocella", "mxaF", "mxaI", "xoxF","rbcL - RuBisCo", "rbcS - RuBisCo")), ordered = TRUE) %>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test"))) %>%
  #mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  ggplot(., aes(x = Gene, y = genome)) +  
  #  geom_tile(aes(fill = Type), color="grey90", linewidth=0.09 ) +
  geom_tile(aes(fill = if_else(is.na(presence), "white", Type)), color="grey90", linewidth=0.09 ) +  
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuBis\nCO")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))




KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()


pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Methylocella_tree$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "Denitrification"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative", "Sulphur"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB"))%>%
  mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "Terminal\noxidase")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))


tree_met <-pl_simple  %>% aplot::insert_left(Methylocella_plot, width=2.12)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=0.77)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.6)


ggsave("./output/Methylocapsa_met_condense_Sylph_drep_25_10_30.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 128,
       width = 186,
       dpi=300)




library(patchwork)
combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))





#################################################################################################
###################### From here, we are working with Methylocystis ############################
#################################################################################################


###############

# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Genus=="g__Methylocystis")%>%
  filter(tip.label %in% unique(KEGG$genome))%>%
  pull(tip.label)

tax2<-tax%>%
  mutate(clade = if_else(between(genome_number, 1, 25), "clade_1", "other"))%>%
  mutate(clade = if_else(between(genome_number, 27, 31), "clade_2", clade))%>%
  mutate(clade = if_else(between(genome_number, 37, 44), "clade_3", clade))

# remove outgroup,I use this instead
Methylocystis_tree <- ape::keep.tip(tr, tip = outgroup_tips)


Methylocystis_plot<-ggtree(Methylocystis_tree, linewidth=0.3) %<+% tax2 + 
  geom_tiplab(aes(label = label_3, fill = clade, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(
    values = c("clade_3" = "#C9BFE3", "clade_2" = "#F0D8C0", "clade_1" = "lightblue", "other" = "grey100"),guide = "none")+
  scale_x_discrete(expand = c(0,0))

Methylocystis_plot

# put on bootstrap values
df1 <- Methylocystis_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


Methylocystis_plot <- Methylocystis_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 


ax<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Methylocystis_plot$data%>%filter(node==min(parent))%>%pull(y)


Methylocystis_plot <- Methylocystis_plot +
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.3,
    inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)

Methylocystis_plot 




######## Plotting the heatmap ######

options("aplot_guides" = "keep")


heat <- sylph %>%
  filter(user_genome %in% Methylocystis_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))




bar_plot <- sylph %>%
  filter(user_genome %in% Methylocystis_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
   # legend.position = "bottom",
    legend.position = "none",
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




#

# heat_raster<-ggrastr::rasterise(heat, layers="Tile", dpi=300)
# bar_raster<-ggrastr::rasterise(bar_plot, layers="Tile", dpi=300)
#library(patchwork)


combined_plot <- heat / bar_plot + patchwork::plot_layout(heights = c(20, 0.5))


# ggsave("output/Methylocystis_with_bar_drep_sylph_25_10_25.svg",combined_plot,
#        units = c("mm"),
#        height = 130,
#        width = 100,
#        dpi=300)



########################## importing pmoA2 line ####################################

otu_merged<-readRDS("output/merged_otu_table_24_06_24.rds")
otu_pmoA2<-otu_merged%>%  select(
  matches("pmoA2"),
  matches("Seq|mfd")
)

otu_pmoA2_long<-pivot_longer(otu_pmoA2, starts_with("Root;"), values_to = "RPKM", names_to = "Tax")%>%
  mutate(Tax_short=trimws(str_extract(Tax, "[^;]*$")))%>%
  mutate(type=if_else(grepl("pxmA", Tax), "pxmA", "pmoA"))%>%
  mutate(type=if_else(grepl("mmoX", Tax), "mmoX", type))%>%
  mutate(type=if_else(grepl("pmoA2", Tax), "pmoA2", type))%>%
  mutate(type=if_else(grepl("Nevskiales_Macondimonas|Binatales|Rhodopila", Tax), "Put. pmoA/mmoX", type))

otu_pmoA2_long<-otu_pmoA2_long %>%
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

otu_pmoA2_long<-otu_pmoA2_long%>%
  filter(mfd_hab1 %in% hab_filter)%>%
  filter(!Tax_short %in% c("Methylosinus_pmoA2", "Methylocystis_pmoA2_1"))%>%
  mutate(Tax_short = factor(
    Tax_short,
    levels = c("Methylocystis_pmoA2_3", "Methylocystis_pmoA2_2", "Methylosinus_Methylocystis_pmoA2"),ordered = TRUE))


heat_pmoA2 <- otu_pmoA2_long %>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = Tax_short, x = SeqId, fill = `RPKM`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 5),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Heath and\nscrub", "Bogs\nand fens", "Freshwater\nsediment","Dunes", "Grassland\nformations", "Forests")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
       # axis.text.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))



combined_plot <- heat / heat_pmoA2 / bar_plot + patchwork::plot_layout(heights = c(20, 1.8, 0.5))
#combined_plot


ggsave("output/Methylocystis_with_bar_pmoA2_drep_sylph_25_12_04.svg",combined_plot,
       units = c("mm"),
       height = 130,
       width = 130,
       dpi=300)


################### Making a dataframe with summarised metabolism ###############

KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step))%>%
  ungroup()

############################################################################################
############### Methylocystis condense met similar to the Methylocella #################
############################################################################################


genes_keep<-KEGG%>%filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Metabolic_step %in% c("Custom\nHMM"))%>%
  group_by(Gene)%>%
  filter(presence==1)%>%
  select(Gene)%>%
  distinct()%>%
  filter(!Gene%in%c("Homologous_MO", "Rhodomicrobium_umbrella"))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))


KOs<-readxl::read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)


KOs<-KOs%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA1_", Gene), "Methylocystis_pmoA1", Gene))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))

library(tidyverse)

CM <- KEGG %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Gene %in% c("rbcL - RuBisCo", "rbcS - RuBisCo") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuBis\nCO", Metabolic_step)) %>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA1_", Gene), "Methylocystis_pmoA1", Gene))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pxmA", "Beijerinckiaceae_pxmA", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA1_1", "Methylosinus_Methylocystis_pmoA1", Gene))%>%
  mutate(Gene=gsub("Methylocystis_pmoA2_1", "Methylosinus_Methylocystis_pmoA2", Gene))%>%
  filter(Gene %in% genes_keep$Gene | Gene %in% c("mxaF", "mxaI", "xoxF"))%>% #, "rbcL - RuBisCo", "rbcS - RuBisCo"))%>% (no Rubisco, have checked)
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuBis\nCO")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        strip.text.x = element_text(size = 5, face = "bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

#CM

pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Methylocystis_tree$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "Denitrification"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative", "Sulphur"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB"))%>%
  mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "Terminal\noxidase")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold"),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

#pl_simple_2


setdiff(pl_simple_2$data$genome, Methylocystis_tree$tip.label)

tree_met <-CM  %>% aplot::insert_left(Methylocystis_plot, width=2.1)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=0.85)

#heat_raster<-ggrastr::rasterise(heat, layers="Tile", dpi=400)

tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.8)
#tree_met_heat_2

ggsave("output/Methylocystis_met_condense_Sylph_drep_25_10_24.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 130,
       width = 190,
       dpi=300)

###






#############################################################
###############  Gammaproteobacteria   ######################
#############################################################


# Identify the tip labels corresponding to the outgroup sequences
outgroup_tips <- tax %>%
  filter(Class=="c__Gammaproteobacteria")%>%
  pull(tip.label)

Gamma_tree <- ape::keep.tip(pr, tip = outgroup_tips)
Gamma_plot<-ggtree(Gamma_tree, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend=F),
              geom = "label", align=TRUE, offset = 0.22, size = 5/.pt,label.size = 0.0,hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:',values = c("GTDB" = "#C9BFE3", "Microflora Danica"="#F0D8C0"),guide = "none")+
  scale_x_discrete(expand = c(0,0))

Gamma_plot


# put on bootstrap values
df1 <- Gamma_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



# remove outgroup,I use this instead


# put on bootstrap values
df1 <- Gamma_tree %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]


Gamma_plot <- Gamma_plot %<+% df_boot + geom_nodepoint(aes(color=status), size=1) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.83),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) +
  guides(fill = guide_legend(
    title = "Genome origin:",
    title.theme = element_text(size=7),
    label.theme = element_text(size=6),
    keyheight = unit(0.3, "cm"),  # Adjust as needed
    keywidth = unit(0.3, "cm"),   # Adjust as needed
    override.aes = list(size = 2)
  ))


ax<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(x)
ay<-Gamma_plot$data%>%filter(node==min(parent))%>%pull(y)

Gamma_plot<-Gamma_plot + geom_segment(aes(x = ax, y = ay, xend = -0.05, yend = ay),
                                      arrow = arrow(length = unit(0.2, "cm"), type="closed"), lwd=0.2)+
  annotate("text", x = -0.025, y = ay+1.5, label = "Outgroup", size=2.5)+
  theme(legend.position = c(0.15, 0.83)) +
  geom_treescale(x=0.035, y=9, label="Tree scale", fontsize =2.5)

Gamma_plot 




heat <- sylph %>%
  filter(user_genome %in% Gamma_tree$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Temperate heath and scrub", "Dunes", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
 # mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  geom_tile() +
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  #scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  facet_nested(~ hab1_label, scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "right",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"),
        strip.text.x = element_text(size = 6.5, face="bold"),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))


tree_sylph<-heat%>% aplot::insert_left(Gamma_plot, width=0.7)


############################# Okay,  need to only select some of the interesting genomes and potentially leave the remeaning in supplementary ######

Gamma_plot <- ggtree(Gamma_tree) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align = TRUE, offset = 0.1, size = 2.1, label.size = 0.0,
              hjust = 1, linesize = 0.2, alpha = 1) + 
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0")) + 
  theme(legend.position = c(0.15, 0.83)) + 
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(aes(label = node), hjust = -0.5, size = 2, color = "black")  # Add node IDs

Gamma_plot

tree_sylph<-heat%>% aplot::insert_left(Gamma_plot, width=0.7)



# Remove node 168,. 220
# Single tips: 13, 19
pruned_tips_f <- reduce(offspring(Gamma_tree, c("168", "220", "228", "216", "196", "191", "151", "177")), c)
Single_tips<-c("13", "19", "31", "33", "30", "35", "38", "39", "54", "51", "52", "53", "4")
pruned_tips<-Gamma_plot$data%>%filter(isTip==TRUE)%>%filter(node %in% union(pruned_tips_f, Single_tips))%>%pull(label)

extra <- c("GCA_026399025.1", "GCF_002934365.1",
           "GCF_030680395.1", "GCA_028712225.1","GCF_002005105.1",
           "GCA_024636875.1", "GCA_016715325.1", "GCA_030679035.1", "GCA_028707085.1")

pruned_tips <- unique(c(pruned_tips, extra))


# remove outgroup,I use this instead
Gamma_tree_pruned <- ape::drop.tip(Gamma_tree, tip = pruned_tips)



# Recreate the plot with the pruned tree
Gamma_plot_pruned <- ggtree(Gamma_tree_pruned) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align = TRUE, offset = 0.1, size = 2.1, label.size = 0.0,
              hjust = 1, linesize = 0.2, alpha = 1) + 
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0")) + 
  theme(legend.position = c(0.15, 0.83)) + 
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(aes(label = node), hjust = -0.5, size = 2, color = "black")  # Add node IDs

# Display the plot
Gamma_plot_pruned


setdiff(Gamma_plot$data$label, Gamma_plot_pruned$data$label)


tree_sylph<-heat%>% aplot::insert_left(Gamma_plot_pruned, width=0.7)


Gamma_plot_pruned <- ggtree(Gamma_tree_pruned, linewidth=0.3) %<+% tax + 
  geom_tiplab(aes(label = label_3, fill = type, show.legend = FALSE),
              geom = "label", align=TRUE, offset = 0.56, size = 5/.pt, hjust = 1, linesize = 0.2, alpha=1) +
  scale_fill_manual(name = 'Genome origin:', values = c("GTDB" = "#C9BFE3", "Microflora Danica" = "#F0D8C0"),guide = "none") + 
  scale_x_discrete(expand = c(0, 0)) 



# put on bootstrap values
df1 <- Gamma_tree_pruned %>% as.treedata %>% as_tibble

df_boot <- df1
df_boot$label <- as.numeric(df_boot$label)
df_boot <- df_boot[!is.na(df_boot$label), ]
df_boot$status <- ifelse((df_boot$label >= 95), "95-100 %",
                         ifelse((df_boot$label >= 85),"85-95 %","0-85 %"))
df_boot <- df_boot[, c("node","status")]



Gamma_plot_pruned <- Gamma_plot_pruned %<+% df_boot + geom_nodepoint(aes(color=status), size=0.5) +
  scale_color_manual(values = c("95-100 %" = "grey0", "85-95 %"="gray40", "0-85 %"="gray75")) +
  labs( color="Bootstrap:")+
  theme(legend.position = c(0.15, 0.87),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5) # Adjusts the inner margin of the legend box
  ) +
  guides(color = guide_legend(
    title = "Bootstrap:",
    title.theme = element_text(size=5),
    label.theme = element_text(size=5),
    keyheight = unit(0.1, "cm"),  # Adjust as needed
    keywidth = unit(0.1, "cm"),   # Adjust as needed
    override.aes = list(size = 1.2)
  )) 

ax<-Gamma_plot_pruned$data%>%filter(node==min(parent))%>%pull(x)
ay<-Gamma_plot_pruned$data%>%filter(node==min(parent))%>%pull(y)

Gamma_plot_pruned<-Gamma_plot_pruned+
  geom_segment(data = data.frame(x = ax, y = ay, xend = -0.03, yend = ay),
   aes(x = x, y = y, xend = xend, yend = yend),
   linewidth = 0.3,
   inherit.aes = FALSE)+
  theme(legend.position = c(0.15, 0.9)) +
  geom_treescale(x=0.02, y=30, label="Tree scale", fontsize =5/.pt)

Gamma_plot_pruned 



heat <- sylph %>%
  filter(user_genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(mfd_hab2 = factor(mfd_hab2, levels = hab2_sort_order)) %>%
  ggplot(aes(y = user_genome, x = SeqId, fill = `Taxonomic_abundance`)) +
  ggrastr::geom_tile_rast() +
  # scale_fill_viridis_c(name = "Taxonomic\nAbundance\n(%)", trans = "sqrt", limits = c(0.00, 2), na.value = "#F8E622")+
  scale_fill_gradientn(
    name = "Taxonomic\nAbundance\n(%)",  # The label
    colors = c( "white","#82A3CD", "darkorange", "darkred", "black"),  # Your custom color gradient
    trans = "sqrt",  # Same as in viridis
    limits = c(0.00, 2),  # Set the limits manually
    na.value = "black"  # Handling NA values
  ) +
  facet_grid(. ~ factor(hab1_label, levels=c("Bogs\nand fens", "Freshwater\nsediment", "Grassland\nformations", "Forests", "Green-\nspaces")), scales = "free", space = "free") +
  #theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 8.2),
        legend.position = "none",
        legend.margin=margin(-1,0,0,0),
        legend.box.margin=margin(-1,0,0,0),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.2, "cm"),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        #  strip.text.y = element_text(size = 8.2, angle = 0),  # Size and orientation for y-axis facet labels
        #   strip.background = element_rect(fill = "grey90", color = "grey90"),
        panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
        text = element_text(family = "Arial"),
        plot.background = element_rect(fill = "transparent"),
        panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))





################################################################################
########################### Nice, then it is metabolism time ###################
################################################################################



genes_keep<-KEGG%>%filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(Metabolic_step %in% c("Custom\nHMM"))%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  group_by(Gene)%>%
  filter(presence==1)%>%
  select(Gene)%>%
  distinct()%>%
  filter(!Gene%in%c("Homologous_MO", "Rhodomicrobium_umbrella"))%>%
  # mutate(Gene=if_else(grepl("Methylocystis_pmoA2_", Gene), "Methylocystis_pmoA2", Gene))%>%
  mutate(Gene=gsub("Methylosinus_umbrella", "Methylosinus", Gene))


KOs<-readxl::read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed , ' - ',  Pathway_step))%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)


KOs<-KOs%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))
  



CM <- KEGG %>% 
  filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(Module == 'C1') %>% 
  filter(Type %in% c("RuMP cycle") | Metabolic_step %in% c("Custom\nHMM", "Methanol\noxidation")) %>%
  mutate(Gene=if_else(Type=="RuMP cycle", Gene_collapsed, Gene))%>%
  mutate(Gene=if_else(grepl("Methylomonadaceae_mmoX", KO), "Methylomonadaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Methylococcales_mmoX|Methyloccocaceae_mmoX", Gene), "Methylococcaceae_mmoX", Gene))%>%
  mutate(Gene=if_else(grepl("Root; USCg", KO), "USCg", Gene))%>%
  mutate(Gene=if_else(Gene=="pxmA", "Methylomonadaceae_pxmA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylomonadaceae", KO), "Methylomonadaceae_pmoA", Gene))%>%
  mutate(Gene=if_else(grepl("Root; o_Methylococcales_pmoA; Methylococcaceae", KO), "Methylococcaceae_pmoA", Gene))%>%
  filter(Gene %in% genes_keep$Gene | Gene %in% c("mxaF", "mxaI", "xoxF", "hxlB/hps-phi", "hxlA/fae-hps/hps-phi"))%>% 
  mutate(Gene=if_else(Gene=="Methylocystis_Methylosinus", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylosinus_2", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene=if_else(Gene=="Methylocystis_1", "Methylocystis_Methylosinus_mmoX", Gene))%>%
  mutate(Gene = fct_relevel(Gene, KOs$Gene)) %>%
  mutate(presence=if_else(presence==0, NA, presence))%>%
  mutate(Metabolic_step = gsub("Methanol\noxidation", "MeOH\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Custom\nHMM", "Methane\noxi.", Metabolic_step)) %>%
  mutate(Metabolic_step = gsub("Carbon\nassimilation", "RuMP", Metabolic_step)) %>%
  # mutate(gene_label=if_else(grepl("Methylocystis", gene_label), Gene, gene_label))%>%
  ggplot(aes(x = Gene, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  facet_grid(. ~ factor(Metabolic_step, levels=c("Methane\noxi.", "MeOH\noxi.", "RuMP")), scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(na.value = "transparent", values = c(methane_colors)) +
  #theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),  
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        strip.clip = "off",
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(c(0,0)) +
  scale_y_discrete(c(0,0))

CM



KEGG2<-KEGG%>%
  mutate(Metabolic_step=gsub("Formaldehyde\noxidation/assimilation", "From.\noxi", Metabolic_step))%>%
  filter(Metabolic_step %in% c("CO and H2", "N fixation", "Sulphur oxi / sulphate red.", "Sulphur oxi.", "PHB", "Complex IV", "Denitrification", "DNRA / denitrification", "From.\noxi"))%>%
  group_by(Type, genome, Module, Metabolic_step)%>%
  summarise(
    total_genes = n_distinct(Gene_collapsed),                    # Total number of unique genes within each Type
    genes_present = n_distinct(Gene_collapsed[presence == 1]),              # Count of genes with Presence == 1
    proportion_present = genes_present / total_genes   # Proportion of genes present within each Type
  ) %>%
  mutate(presence=if_else(proportion_present<0.6, NA, 1))%>%
  mutate(Metabolic_step=gsub("Sulphur oxi.", "S\noxi.", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("sulphate red.", "red.", Metabolic_step),
         Metabolic_step=gsub("CO and H2", "Trace\ngas oxi.", Metabolic_step),
         Metabolic_step=gsub("Complex IV", "Complex\nIV", Metabolic_step),
         Metabolic_step=gsub("N fixation", "N\nfix", Metabolic_step),
         Metabolic_step=gsub("Denitrification", "N", Metabolic_step),
         Metabolic_step=gsub("DNRA / denitrification", "N", Metabolic_step))%>%
  ungroup()



pl_simple_2 <- KEGG2 %>% 
  filter(genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(! Type %in% c("NiFe hydrogenase", "membrane-bound hydrogenase", "ferrodoxin hydrogenase", "quinone-reactive Ni/Fe-hydrogenase", "NADP-reducing hydrogenase", "aprAB", "DsrABL", "DsrC", "dsrMKJOP", "dsrEFH"))%>%
  filter(!Module%in% c("Oxi. Phos. Alternative"))%>%
  filter(!Metabolic_step %in% c("From.\noxi", "PHB", "Complex\nIV"))%>%
 # mutate(Metabolic_step=gsub("Complex\nIV", "Terminal\noxidase", Metabolic_step))%>%
  mutate(Metabolic_step=gsub("N\nfix", "N", Metabolic_step))%>%
  mutate(Metabolic_step=if_else(grepl("S\noxi.", Metabolic_step), "S oxi.\n/red.", Metabolic_step))%>%
  # mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  mutate(Type = factor(Type, levels=unique(KOs$Type), ordered = TRUE)) %>%
  mutate(Type = gsub("Carbon monoxide oxidation", "CO oxidation", Type))%>%
  ggplot(., aes(x = Type, y = genome))+ 
  geom_tile(aes(fill = if_else(is.na(presence),  "white", Type)), color="grey90", linewidth=0.09 ) + 
  #facet_grid(. ~Metabolic_step, scales = "free", space = "free") +
  facet_nested(. ~ factor(Metabolic_step, levels=c("Trace\ngas oxi.", "N\nfix", "N", "Terminal\noxidase", "S oxi.\n/red.")), scales = "free", space = "free") +
  scale_fill_manual(na.value="transparent", values=c(
    "CO oxidation"="#52338D",
    "[NiFe] Group 1"="#B58EFC",
    "[NiFe] Group 2"="#85C1E9",
    "[NiFe] Group 3d"="#1F77B4",
    "N fixation"="#B59F62",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "N fixation"="#B59F62",
    "cyt. c aa3-type oxidase"="#8CB88F",
    "High affin. cyt.bd-type"="#407DA8",
    "Cyt. bo3 ubiquinol oxi."="#3FACB8",
    "High affin. cyt c bb3-type"="darkblue",
    "Nitrification"="#EBBC3E",
    "Denitrification"="#EA8D17",
    "DNRA / denitrification"="#CF4F48",
    "DNRA / assim. nitrate reduction"="#EB4D22",
    "Assimilatory nitrate reduction"="darkred",
    "Assimilatory sulphate reduction"="#0DB7A8",
    "sat"="#6C961E",
    "SOX"="#92C76C"
  ))+
  #theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), 
        axis.text.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.ticks.y = element_blank(),
        # axis.text.y=element_text(size=7),
        strip.text.x = element_text(size = 5, face="bold", margin = margin(t = 1, b = 1, unit = "pt")),
        strip.background = element_blank(),
        axis.title = element_blank(),
        strip.clip = "off",
        text = element_text(family = "Arial"),
        legend.position = "none",
        panel.spacing = unit(0.03, "cm", data = NULL),
        legend.title=element_blank(),
        panel.background = element_rect(fill="white")
  )+
  scale_x_discrete(c(0,0))+
  scale_y_discrete(c(0,0))

pl_simple_2


setdiff(Gamma_plot_pruned$data$label, pl_simple_2$data$genome)

tree_met <-CM  %>% aplot::insert_left(Gamma_plot_pruned, width=3.1)
tree_met_2 <-tree_met  %>% aplot::insert_right(pl_simple_2, width=1)
tree_met_heat_2<-tree_met_2 %>% aplot::insert_right(heat, width=3.8)
#tree_met_heat_2


ggsave("./output/Gammas_met_condense_Sylph_drep_25_11_10.png",
       tree_met_heat_2,
       units = c("mm"),
       height = 160,
       width = 190,
       dpi=300)


ggsave("./output/Gammas_met_condense_Sylph_drep_25_11_10.svg",
       tree_met_heat_2,
       units = c("mm"),
       height = 160,
       width = 190,
       dpi=300)




#


bar_plot <- sylph %>%
  filter(user_genome %in% Gamma_tree_pruned$tip.label)%>%
  filter(mfd_hab1 %in% c("Bogs, mires and fens", "Grassland formations", "Forests","Greenspaces", "Freshwater"))%>% ###### Okay, here is where I need to filter the habitats
  filter(!mfd_hab2=="Mire")%>%
  filter(!mfd_sampletype=="Water")%>%
  mutate(hab1_label=gsub("Bogs, mires\nand fens", "Bogs\nand fens", hab1_label))%>%
  mutate(SeqId = factor(SeqId, levels = levels_hab2, ordered = TRUE)) %>%
  ggplot(aes(x = SeqId, fill = mfd_hab2)) +
  ggrastr::geom_tile_rast(aes(y = 1)) +
  facet_grid(. ~ factor(hab1_label, levels=c("Bogs\nand fens", "Freshwater\nsediment", "Grassland\nformations", "Forests", "Green-\nspaces")), scales = "free", space = "free") +
  scale_fill_manual(values = palette_mfd_hab2) + # Specify your colors
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(colour="black", fill=NA, linewidth = 0.2) ,
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "transparent"),
    legend.position = "none",
    panel.spacing = unit(0.03, "cm", data = NULL))+
  scale_y_continuous(expand = c(0,0)) 

combined_plot <- heat / bar_plot + plot_layout(heights = c(20, 0.5))


ggsave("./output/Gamma_Sylph_bar_25_11_10.svg",combined_plot,
       units = c("mm"),
       height = 110,
       width = 86,
       dpi=300)



