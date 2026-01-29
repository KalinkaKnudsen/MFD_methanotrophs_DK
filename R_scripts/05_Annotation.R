#!/usr/bin/env Rscript

.libPaths(c("/home/bio.aau.dk/vj52ou/software/R_packages.v.4.3.2", .libPaths()))
#library(ggplot2)
library(cowplot)
library(vroom)
library(tidyverse)
library(readxl)
#library(patchwork)

setwd("~/scripts/MFD/methanotrophs/R_scripts/")

HMM_ko<-readRDS("HMM_KOs_25_06_06.rds")
# HMM_old<-readRDS("HMM_KOs.rds")
# 
# x<-setdiff(HMM_old$fasta, HMM_ko$fasta)
# write.table(x, "/home/bio.aau.dk/vj52ou/data/MFD/annotation_combined/graftM_MAGs_24_09_12/genomes_missing_24_09_23.txt", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)


HMM_ko$ko_id <- gsub("\u00A0", " ", HMM_ko$ko_id)


KEGG<-vroom("~/data/MFD/annotation_combined/annotations_concat_25_10_10.tsv",delim="\t", col_select=c("fasta", "ko_id")) %>% ##### Change version of KEGG here!
  mutate(fasta=gsub("_genomic","",fasta))%>%
  filter(!is.na(ko_id))%>%
  rbind(HMM_ko)



genomes_KEGG <- unique(KEGG$fasta)


KOs<-read_excel("KO_KSK_25_08_12.xlsx") %>%
  mutate (Pathway = gsub ("Formaldehyde oxidation" , "Form.ald.oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formaldehyde assim. H4F" , "Form.ald. assim. H4F", Pathway))%>%
  mutate (Pathway = gsub ("Methane oxidation" , "Methane oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Methanol oxidation" , "Methanol oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Formate oxidation" , "Formate oxi.", Pathway))%>%
  mutate (Pathway = gsub ("Carbon monoxide oxidation" , "CO oxi.", Pathway))%>%
  mutate (Pathway = gsub ("RuMP cycle" , "RuMP", Pathway))%>%
  mutate (Pathway = gsub ("Serine cycle" , "Serine", Pathway))%>%
  mutate (Pathway = gsub ("CBB cycle" , "CBB", Pathway)) %>%
  mutate (gene_label = paste0(Gene_collapsed, ' - ',  Pathway_step))

KOs$KO <- gsub("\u00A0", " ", KOs$KO)

#Mutating the facets
KOs<-KOs%>%
  mutate(Metabolic_step = gsub("Methane oxidation", "Methane\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Methanol oxidation", "Methanol\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formaldehyde oxidation/assimilation", "Formaldehyde\noxidation/assimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Formate oxidation", "Formate\noxidation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Carbon assimilation", "Carbon\nassimilation", Metabolic_step))%>%
  mutate(Metabolic_step = gsub("Custom HMM", "Custom\nHMM", Metabolic_step))



### Converting into wide format 
### For KEGG
df <- data.frame()
for (i in 1:length(genomes_KEGG)) {
  d <- KOs %>% mutate(genome = genomes_KEGG[i])
  df = rbind(df, d)
}
# make presence column
KEGG <- KEGG %>% mutate(presence = 1)
KEGG <- left_join(df, KEGG, by = c('KO'='ko_id','genome' ='fasta'), relationship = "many-to-many")
# convert NA into 0's
KEGG$presence[is.na(KEGG$presence)] <- 0



########## Sorting the issue of several hits with one being "0" and the other being "1" with the same gene_label.
### An issue for the custom HMM ####



cleaned_KEGG <- KEGG %>%
  filter(presence %in% c(0, 1)) %>%
  group_by(genome, gene_label) %>%
  filter(all(c(0, 1) %in% presence)) %>%   # Keep groups with both 0 and 1
  filter(presence == 1) %>%                # Keep only the rows where presence is 1
  ungroup()


non_mixed <- KEGG %>%
  anti_join(
    KEGG %>%
      filter(presence %in% c(0, 1)) %>%
      group_by(genome, gene_label) %>%
      filter(all(c(0, 1) %in% presence)) %>%
      ungroup(),
    by = c("genome", "gene_label", "presence")
  )

# Combine cleaned + untouched data
final_KEGG <- bind_rows(non_mixed, cleaned_KEGG)







###### Adding taxonomy

# 
# V4<-vroom("/projects/microflora_danica/deep_metagenomes/mags_v4/analysis/gtdb_r220/gtdb_class/classify/gtdbtk.bac120.summary.tsv", delim="\t")%>%
#   select(user_genome, classification)
# SR<-vroom("/projects/microflora_danica/shallow_mags/analysis/gtdb_r220/gtdb_class/classify/gtdbtk.bac120.summary.tsv", delim="\t")%>%
#   select(user_genome, classification)
# gtdb<-vroom("/home/bio.aau.dk/vj52ou/data/gtdb/bac120_taxonomy_r220.tsv", delim="\t", col_names = c("user_genome","classification"))%>%
#   mutate(user_genome=gsub("RS_","", user_genome),
#          user_genome=gsub("GB_","", user_genome))
# DCM<-vroom("/home/bio.aau.dk/vj52ou/data/gtdb/methanotrophs/gtdb_suspected_genomes/gtdb_tk_DMC_bin26/classify/gtdbtk.bac120.summary.tsv", delim="\t")%>%
#   select(user_genome, classification)
# 
# tax<-rbind(V4, SR, gtdb, DCM)%>%
#   separate(classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ';', remove = FALSE)%>%
#   mutate(Species=if_else(user_genome=="DMC_bin26", "s__(Ca. M. kingii)", Species))


saveRDS(tax, "tax_MAGs_gtdb.rds")

tax<-readRDS("tax_MAGs_gtdb.rds")
drep<-read_lines("~/data/MFD/annotation_combined/drep_list.txt")


# tax<-vroom("/projects/microflora_danica/deep_metagenomes/mags_v4/MFD_mags_drep_V4.tsv", delim = "\t")%>%
#   select(bin, gtdb_classification, Completeness_CheckM2)%>%
#   filter(bin %in% genomes_KEGG)%>%
#   separate(gtdb_classification, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ';', remove = FALSE) %>%
#   mutate(Quality = case_when(
#     Completeness_CheckM2 > 90 ~ "91-100",
#     Completeness_CheckM2 > 80 ~ "81-90",
#     Completeness_CheckM2 > 70 ~ "71-80",
#     Completeness_CheckM2 > 60 ~ "61-70",
#     Completeness_CheckM2 > 50 ~ "51-60",
#     TRUE ~ "<=50"
#   ))%>%
#   rename(user_genome=bin)%>%
#   rename(classification=gtdb_classification)%>%
#   select(!Completeness_CheckM2)





KEGG <- left_join(final_KEGG, tax, by =c('genome'='user_genome')) %>%
  mutate(label = paste0(if_else(Genus == "g__", paste0(Family, ";", Genus), if_else(Species == "s__", paste0(Genus,";",Species), Species)),  " | ", genome ))%>%
  mutate(drep=if_else(genome %in% drep, "yes", "no"))%>%
  mutate(drep = if_else(grepl("GCA|GCF", genome), "yes", drep))

#saveRDS(KEGG, "./output/metabolism/KEGG_24_09_13.rds")


################ Importing and sorting sulphur stuff #######################

hydDB<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/annotation_combined/hybDB/hyddb-results_24_09_19_v2.csv", delim=";", col_names = c("sequence", "db_hit"))

hydDB$KO <- str_extract(hydDB$sequence, "ko:K\\d+")
hydDB$genome <- str_extract(hydDB$sequence, "^[^_]+(?:_[^_]+)*(?=_(k127|contig|genomic|NZ|NC))")

hydDB_filt<-hydDB%>%
  mutate(KO=gsub("ko:", "", KO))%>%
 # select(ko_number, db_hit)%>%
  group_by(KO, db_hit) %>%  # Group by ko_number and db_hit
  summarise(count = n(), .groups = 'drop')%>%
  arrange(db_hit)

## From above, we need to keep: 
KO_keep<-c("K00436", "K23549", "K06281")

hydDB<-hydDB%>%
  mutate(KO=gsub("ko:", "", KO))%>%
  mutate(genome=gsub("_genomic", "", genome))%>%
  filter(KO %in% KO_keep)%>%
  filter(!db_hit %in% c("NONHYDROGENASE"))%>%
  select(!sequence)%>%
  distinct()

KEGG2 <- KEGG %>%
  distinct()%>%
  left_join(hydDB, by = c("KO", "genome"))%>%
  mutate(Type=if_else(is.na(db_hit), Type, db_hit))%>%
  select(!db_hit)



KEGG2<-KEGG2%>%
  mutate(Type=if_else(grepl("hox", Gene), "[NiFe] Group 3d", Type))%>%
  mutate(Type=if_else(grepl("hya", Gene), "[NiFe] Group 1", Type))%>%
  mutate(Type=if_else(grepl("hupU|hupV", Gene), "[NiFe] Group 2", Type))



#saveRDS(KEGG, "./output/metabolism/KEGG_24_09_13.rds")




###### DsrMKJOP also #######

DsrMKJOP<-vroom("/home/bio.aau.dk/vj52ou/data/MFD/annotation_combined/dsrMKJOP/e10.tsv", delim="\t", col_names = c("sequence", "db_hit"))

DsrMKJOP$genome <- str_extract(DsrMKJOP$sequence, "^[^_]+(?:_[^_]+)*(?=_(k127|contig|genomic|NZ|NC))")
DsrMKJOP<-DsrMKJOP%>%
  select(!sequence)%>%
  mutate(KO = case_when(
    db_hit == "DsrM" ~ "K27187",
    db_hit == "DsrK" ~ "K27188",
    db_hit == "DsrJ" ~ "K27189",
    db_hit == "DsrO" ~ "K27190",
    db_hit == "DsrP" ~ "K27191",
    TRUE ~ "not"  # default case when none of the conditions match
  ))%>%
  mutate(db_hit=gsub("Dsr", "dsr", db_hit))%>%
  mutate(Gene=db_hit)%>%
  distinct()


KEGG2 <- KEGG2 %>%
  left_join(DsrMKJOP, by = c("genome", "KO", "Gene"))%>%
  mutate(presence=if_else(is.na(db_hit), presence, 1))%>%
  select(!db_hit)



### Change date here below ###
saveRDS(KEGG2, "./output/metabolism/KEGG_25_12_02.rds")
#KEGG<-readRDS("./output/metabolism/KEGG_25_09_02.rds")

count<-KEGG%>%filter(Metabolic_step=="Custom\nHMM")%>%filter(presence>0)%>%
  filter(grepl("MFD|LIB", genome))%>%select(genome, Species, Family, Genus, Gene_collapsed, drep)%>%
  filter(!grepl("Burkholderiaceae|Homologous*|Binatales_JAKAVN01|Nitrosomonadaceae|Nevskiales|Mycobacterium|Nevskiales_Macondimonas_put_pmoA", Gene_collapsed))%>%
  filter(!grepl("f__Acetobacteraceae|f__Acidobacteriaceae|f__Nitrosomonadaceae|f__Gallionellaceae|f__Nitrospiraceae", Family))%>%
  select(genome, Species, Genus, Family, drep)%>%
  distinct()

count%>%count(Family)
count%>%count(drep)
count%>%count(Species, drep)
length(count$genome)



################## RUn above to create new KEGG dataframe #############################





#######################################################################
############# Saving a KEGG df with only methanotrophs ################
#######################################################################


KEGG<-readRDS("./output/metabolism/KEGG_25_12_02.rds")

genomes_keep<-KEGG %>% 
  filter(!is.na(Genus))%>%
  filter(Module == 'C1') %>% 
  filter(Metabolic_step=="Custom\nHMM")%>%
  group_by(Gene_collapsed)%>%
  filter(!Gene_collapsed %in% c("Nitrosomonadaceae", "Betaproteobacteria_amoA", "Nitrosomonas", "Nitrospira_clade_A", "Nitrospira_clade_B", "Actinobacteria", "Nitrosococcus", "Propane_MO_Actino_cluster", "Cycloclasticus"))%>%
  filter(presence==1)%>%
  # filter(!Genus %in% genus_remove_string)%>%
  distinct()



KEGG_filt<-KEGG%>%
  filter(genome %in% genomes_keep$genome)


DRAM_raw<-vroom("~/data/MFD/annotation_combined/annotations_concat_25_10_10.tsv",delim="\t") %>%
  mutate(fasta=gsub("_genomic","",fasta))%>%
  filter(fasta %in% genomes_keep$genome)

DRAM_raw<-DRAM_raw%>%
  rename(seqid='...1')

diff<-setdiff(KEGG_filt$genome, DRAM_raw$fasta)
check<-KEGG_filt%>%filter(genome %in% diff)%>%filter(presence==1)

### The ones missing look ok!

write.table(orgs, file = "C:/Users/orgs_updated.tsv", row.names=FALSE, sep="\t")
#openxlsx::write.xlsx(DRAM_raw, "output/metabolism/DRAM_output.xlsx", sheetName = "Sheet1", rowNames = FALSE)
write.table(DRAM_raw, "output/metabolism/DRAM_output.tsv", row.names = F, sep="\t")
write.table(KEGG_filt, "output/metabolism/DRAM_tidy.tsv", row.names = F, sep="\t")

#######################################################################









################## RUn above to create new KEGG dataframe #############################









#############################
#Getting genomes within Rhodomicrobium to do the dsrAB tree:

dsr<-KEGG%>%
  filter(KO %in% c("K11180", "K11181"))%>%
  filter(presence==1)%>%
  select(genome)%>%
  distinct()
  


readr::write_csv(dsr, "~/data/MFD/annotation_combined/dsrAB/genomes.txt", col_names = FALSE)

######################### Now to the plotting ########################
#custom_colors <-c("#679b60","#5e4fa2" , "#679b60", "turquoise4", "darkred","#d97512","#c46ca1" , "#b3943c", "grey75","darkblue","darkred","darkgreen" )
methane_colors<-readRDS("palette_methane_colors.rds")

KEGG_sort<-KEGG%>%
  arrange(Genus)

pl <- KEGG %>% 
  filter(Module == 'C1') %>% 
  filter(drep=="yes")%>%
  filter(Genus=="g__Rhodomicrobium")%>%
  filter(!Species %in% c("s__Methylocella tundrae",
                         #  "s__Methylocella acidiphila",
                         #   "s__Methylocystis silviterrae",
                         "s__Methylomonas methanica_A",
                         #  "s__Methylomonas methanica_B",
                         #   "s__Methylocella silvestris_A",
                         #  "s__Rhodomicrobium sp003153975",
                         "s__Methylocella aurea",
                         "s__Methylocella sp004564215"))%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
 # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  mutate(label = fct_relevel(label, KEGG_sort$label)) %>%
  ggplot(., aes(x = gene_label, y = label))+ 
  geom_tile(data = . %>% filter(presence > 0),
            aes(fill = Type), color="grey90", linewidth=0.1) +  # Set a fixed width for the tiles (e.g., 0.5)
  facet_grid(~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(values=methane_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_text(size=8.3),
        axis.title = element_blank(),legend.position = "none" # This removes the legend)
  )


pl


pl <- KEGG %>% 
  filter(Module == 'C1') %>% 
  filter(Genus=="g__Methylocella")%>%
  filter(drep=="yes")%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  mutate(label = fct_relevel(label, KEGG_sort$label)) %>%
  ggplot(., aes(x = gene_label, y = label))+ 
  geom_tile(data = . %>% filter(presence > 0),
            aes(fill = Type), color="grey90", linewidth=0.1) +  # Set a fixed width for the tiles (e.g., 0.5)
  facet_grid(~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(values=methane_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_text(size=8.3),
        axis.title = element_blank(),legend.position = "none" # This removes the legend)
  )


pl


# 
# ggsave("output/metabolism/Rhodomicrobium_non_MMO.png",
#        pl,
#        height = 10,
#        width = 15)
# 
# 


##################### exporting the genomes I only have HMM output from ###########
DRAMMED<-read_lines("~/data/MFD/annotation_combined/DRAM_24_07_23/genomes_screened.txt")
tax_filter<-read_lines("/home/bio.aau.dk/vj52ou/data/MFD/deep_metagenomes/methanotrophs/taxonomy_filter_r220.txt")

pl <- KEGG %>% 
  filter(Module == 'C1') %>% 
  #filter(drep=="no")%>%
  filter(Genus=="g__Rhodomicrobium")%>%
  filter(genome %in% DRAMMED)%>%
  mutate(Metabolic_step = fct_relevel(Metabolic_step, KOs$Metabolic_step)) %>%
  # mutate(Quality=fct_relevel(Quality,c("HQ","MQ","Control","Test")))%>%
  mutate(gene_label = fct_relevel(gene_label, KOs$gene_label)) %>%
  mutate(label = fct_relevel(label, KEGG_sort$label)) %>%
  ggplot(., aes(x = gene_label, y = label))+ 
  geom_tile(data = . %>% filter(presence > 0),
            aes(fill = Type), color="grey90", linewidth=0.1) +  # Set a fixed width for the tiles (e.g., 0.5)
 # facet_grid(~Metabolic_step, scales = "free", space = "free") +  # Display Pathway names above the plot
  scale_fill_manual(values=methane_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
        axis.text.y = element_text(size=8.3),
        axis.title = element_blank(),legend.position = "none" # This removes the legend)
  )


pl

pattern_combined <- paste(tax_filter, collapse = "|")

investi <- KEGG %>%
  filter(!genome %in% DRAMMED) %>%
  filter(presence > 0) %>%
  filter(grepl(pattern_combined, classification))%>%
  select(genome)%>%
  distinct()

readr::write_csv(investi, "~/data/MFD/annotation_combined/DRAM_24_07_23/genomes_missing.txt", col_names = FALSE)

############################### Creating a dataframe for the genes I have searched for ###################

tax_remove<-c("Root; pmoA_amoA_pxmA; Homologous_MO", "Root; Homologous_pmoA; Homologous_Binatales", "Root; pmoA_amoA_pxmA; amoA; Nitrospira_clade_B", "Root; pmoA_amoA_pxmA; Nevskiales_Macondimonas", "Root; pmoA_amoA_pxmA; Nevskiales_Macondimonas; Nevskiales",
              "Root; pmoA_amoA_pxmA; amoA; Nitrospira_clade_A", "Root; Homologous_pmoA", "Root; pmoA_amoA_pxmA; amoA; Betaproteobacteria_amoA", "Root; pmoA_amoA_pxmA; Nitrosococcus", "Root; Homologous_pmoA; Actinobacteria", "Root; pmoA_amoA_pxmA","Root")
filter<-read_lines("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/analysis_mags_v4_derep_23_10_17/selected_genomes_e10_pmoA_filter.txt")
V4_pmoA<-vroom("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/analysis_mags_v4_derep_23_10_17/combined_count_table_e10_pmoA.txt", delim="\t")%>%
 # select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% tax_remove)
  
V4_pmoA<-V4_pmoA%>%
pivot_longer(cols=2:length(V4_pmoA), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/mags_v4_all/23_12_05_pmoA/selected_genomes_e10.txt", skip=1)
V4_all<-vroom("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/mags_v4_all/23_12_05_pmoA/combined_count_table_e10.tsv", delim="\t")%>%
  #select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% tax_remove)
V4_all<-V4_all%>%
  pivot_longer(cols=2:length(V4_all), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/MFD/short_read_MAGs/methanotrophs/graftm_output/23_12_04_pmoA/selected_genomes_e10.txt", skip=1)
SR_all<-vroom("~/data/MFD/short_read_MAGs/methanotrophs/graftm_output/23_12_04_pmoA/combined_count_table_e10.txt", delim="\t")%>%
  #select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% tax_remove)

SR_all<-SR_all%>%
  pivot_longer(cols=2:length(SR_all), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/gtdb/methanotrophs/gtdb_r220_24_05_03/pmoA/selected_genomes_e10.txt")
gtdb<-vroom("~/data/gtdb/methanotrophs/gtdb_r220_24_05_03/pmoA/combined_count_table_e10.txt", delim="\t")%>%
  #select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% tax_remove)
gtdb<-gtdb%>%
  pivot_longer(cols=2:length(gtdb), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/gtdb/methanotrophs/gtdb_r214_23_06_22/selected_genomes_pmoA_filter.txt")
gtdb_2<-vroom("~/data/gtdb/methanotrophs/gtdb_r214_23_06_22/combined_count_table_pmoA.txt", delim="\t")%>%
 # select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% tax_remove)
gtdb_2<-gtdb_2%>%
  pivot_longer(cols=2:length(gtdb_2), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

pmoA<-rbind(V4_pmoA, V4_all, SR_all, gtdb, gtdb_2)
rm(V4_pmoA, V4_all, SR_all, gtdb, gtdb_2)

unique(pmoA$ConsensusLineage)

#mmoX#
mmoX_tax_remove<-c("Root; Homologous_mmoX","Root; Homologous_mmoX; Propane_MO", "Root; Homologous_mmoX; Actinobacteria", "Root; Homologous_mmoX; Mycobacterium", "Root")
filter<-read_lines("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/analysis_mags_v4_derep_23_10_17/selected_genomes_e10_pmoA_filter.txt")
V4_mmoX<-vroom("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/analysis_mags_v4_derep_23_10_17/combined_count_table_e10_mmoX.txt", delim="\t")%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)

V4_mmoX<-V4_mmoX%>%
  pivot_longer(cols=2:length(V4_mmoX), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/mags_v4_all/23_12_05_mmoX/selected_genomes_e10.txt", skip=1)
V4_all<-vroom("~/data/MFD/deep_metagenomes/methanotrophs/graftm_output/mags_v4/mags_v4_all/23_12_05_mmoX/combined_count_table_e10.tsv", delim="\t")%>%
#  select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)
V4_all<-V4_all%>%
  pivot_longer(cols=2:length(V4_all), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/MFD/short_read_MAGs/methanotrophs/graftm_output/23_12_04_mmoX/selected_genomes_e10.txt", skip=1)
SR_all<-vroom("~/data/MFD/short_read_MAGs/methanotrophs/graftm_output/23_12_04_mmoX/combined_count_table_e10.txt", delim="\t")%>%
 # select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)

SR_all<-SR_all%>%
  pivot_longer(cols=2:length(SR_all), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()


filter<-read_lines("~/data/gtdb/methanotrophs/gtdb_r220_24_05_03/mmoX/selected_genomes_e10.txt")
gtdb<-vroom("~/data/gtdb/methanotrophs/gtdb_r220_24_05_03/mmoX/combined_count_table_e10.txt", delim="\t")%>%
 # select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)

gtdb<-gtdb%>%
  pivot_longer(cols=2:length(gtdb), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

filter<-read_lines("~/data/gtdb/methanotrophs/gtdb_r214_23_06_22/selected_genomes_mmoX_filter.txt")
gtdb_2<-vroom("~/data/gtdb/methanotrophs/gtdb_r214_23_06_22/combined_count_table_mmoX_filter.txt", delim="\t")%>%
 # select(ConsensusLineage, all_of(filter))%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)
gtdb_2<-gtdb_2%>%
  pivot_longer(cols=2:length(gtdb_2), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()



mmoX<-rbind(V4_mmoX, V4_all, SR_all, gtdb, gtdb_2)
rm(V4_mmoX, V4_all, SR_all, gtdb, gtdb_2)

unique(mmoX$ConsensusLineage)

HMM_ko<-rbind(pmoA, mmoX)%>%
  select(ConsensusLineage, bin)%>%
  rename(fasta=bin)%>%
  rename(ko_id=ConsensusLineage)%>%
  mutate(fasta=gsub("_genomic","", fasta))

saveRDS(HMM_ko, "./HMM_KOs.rds")




######################### Update with new HMMs ##################################



#tax_remove<-c("Root; pmoA_amoA_pxmA; Homologous_MO", "Root; Homologous_pmoA; Homologous_Binatales", "Root; pmoA_amoA_pxmA; amoA; Nitrospira_clade_B", "Root; pmoA_amoA_pxmA; Nevskiales_Macondimonas", "Root; pmoA_amoA_pxmA; Nevskiales_Macondimonas; Nevskiales","Root; pmoA_amoA_pxmA; amoA; Nitrospira_clade_A", "Root; Homologous_pmoA", "Root; pmoA_amoA_pxmA; amoA; Betaproteobacteria_amoA", "Root; pmoA_amoA_pxmA; Nitrosococcus", "Root; Homologous_pmoA; Actinobacteria", "Root; pmoA_amoA_pxmA","Root")
pmoA<-vroom("~/data/MFD/annotation_combined/graftM_MAGs_24_09_12/pmoA/combined_count_table_e10.txt", delim="\t")
pmoA_gtdb<-vroom("~/data/gtdb/methanotrophs/gtdb_r220_25_02_14/pmoA/combined_count_table_e10.txt", delim="\t")
pmoA_gtdb_2<-vroom("/home/bio.aau.dk/vj52ou/data/gtdb/methanotrophs/gtdb_r214_23_06_22/combined_count_table_e10_pmoA.txt", delim="\t")
pmoA_DMC<-vroom("/home/bio.aau.dk/vj52ou/data/gtdb/methanotrophs/gtdb_suspected_genomes/graft_DMC_bin26/combined_count_table.txt", delim="\t")%>%select(!'#ID')

pmoA<-pmoA%>%left_join(pmoA_DMC)%>%mutate(DMC_bin26=if_else(is.na(DMC_bin26), 0, DMC_bin26))

pmoA<-pmoA%>%
  pivot_longer(cols=2:length(pmoA), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

# Find kolonnerne, der skal fjernes
cols_to_remove <- intersect(names(pmoA_gtdb)[-1], names(pmoA_gtdb_2)[-1])

# Fjern kun de kolonner, der findes i pmoA_gtdb_2
pmoA_gtdb_2 <- pmoA_gtdb_2 %>%
  select(-one_of(cols_to_remove))

pmoA_gtdb_2 <- pmoA_gtdb_2 %>%
  pivot_longer(cols=2:length(pmoA_gtdb_2), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

pmoA_gtdb <- pmoA_gtdb %>%
  pivot_longer(cols=2:length(pmoA_gtdb), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

pmoA<-rbind(pmoA, pmoA_gtdb_2, pmoA_gtdb)
unique(pmoA$ConsensusLineage)

#mmoX#
mmoX_tax_remove<-c("Root; Homologous_mmoX","Root; Homologous_mmoX; Propane_MO", "Root; Homologous_mmoX; Actinobacteria", "Root; Homologous_mmoX; Mycobacterium", "Root")
mmoX<-vroom("~/data/MFD/annotation_combined/graftM_MAGs_24_09_12/mmoX/combined_count_table_e10.txt", delim="\t")%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)

mmoX<-mmoX%>%
  pivot_longer(cols=2:length(mmoX), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()

mmoX_gtdb<-vroom("~/data/gtdb/methanotrophs/gtdb_r220_25_02_14/mmoX/combined_count_table_e10.txt", delim="\t")%>%
  filter(!ConsensusLineage %in% mmoX_tax_remove)


mmoX_gtdb<-mmoX_gtdb%>%
  pivot_longer(cols=2:length(mmoX_gtdb), names_to="bin", values_to="count")%>%
  filter(count>0)%>%
  distinct()


mmoX<-mmoX%>%rbind(mmoX_gtdb)


HMM_ko<-rbind(pmoA, mmoX)%>%
  select(ConsensusLineage, bin)%>%
  rename(fasta=bin)%>%
  rename(ko_id=ConsensusLineage)%>%
  mutate(fasta=gsub("_genomic","", fasta))

saveRDS(HMM_ko, "./HMM_KOs_25_06_06.rds")

# 
# 
# 
# 
# methane_colors <- structure(list(`pMMO` = "#b3943c",
#                                  `sMMO` = "darkgreen",
#                                  `calcium (mxa)` = "#5e4fa2",
#                                  `lanthanide (xoxF)` = "#c46ca1",
#                                  putative = "grey75",
#                                  H4MPT = "#d97512",
#                                  GSH = "turquoise4",
#                                  `H4F` = "darkred",
#                                  `Serine cycle` = "darkred",
#                                  `RuMP cycle` = "darkblue" ,
#                                  `CBB cycle` = "#679b60",
#                                  `-` = "#679b60"), class = "character",
#                             row.names = unique_types)
#saveRDS(methane_colors, "palette_methane_colors.rds")


