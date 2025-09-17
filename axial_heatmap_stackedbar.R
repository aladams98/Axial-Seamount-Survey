metadata_input <- read_csv("~/R/Axial2023/datainput/axial_compiledmetadata_revised_samplestokeep.csv") %>%
  rename(sample = NEW_FASTQ_NAME) #%>% column_to_rownames("sample")

library(readxl)
axial_decontam_DF <- read_excel("R/Axial2023/datainput/axial_decontam_DF.xlsx") 

#Filter for samples in axial_decontam_DF that are only in metadata_input
axial_samplestokeep_DF <- axial_decontam_DF %>% filter(sample %in% metadata_input$sample) %>% 
  select(-Taxon) %>% filter(Domain == "Eukaryota") %>% mutate(phylum_class = paste(Phylum, Class, sep = "_")) %>% 
  mutate(assoc_vent = case_when(
    (Location_Name == "ASHES") ~ "ASHES Plume",
    (Location_Name == "International District") ~ "International District Plume",
    (Location_Name == "Diva") ~ "International District Vent",
    (Location_Name == "Gollum") ~ "ASHES Vent",
    (Location_Name == "Marker 33") ~ "International District Vent",
    (Location_Name == "Marker 113") ~ "International District Vent",
    (Location_Name == "Anemone") ~ "ASHES Vent",
    (Location_Name == "Background") ~ "Background",
    (Location_Name == "Bag City") ~ "International District Vent",
    (Location_Name == "Vixen") ~ "International District Vent")) 

#Bring over depth category information to new DF from old metadata_input
axial_samplestokeep_DF <- axial_samplestokeep_DF %>% 
  left_join(metadata_input %>% select(sample, DEPTH_CATEGORY),
            by = "sample")

axial_RelAbun_byVent <- axial_samplestokeep_DF %>%
  # Organize taxonomy
  mutate(assoc_vent_wDepth = case_when(
    grepl("Vent", assoc_vent) ~ assoc_vent,
    grepl("Plume", assoc_vent) ~ paste(assoc_vent, DEPTH_CATEGORY, sep = "_"),
    grepl("Background", assoc_vent, ignore.case = TRUE) ~ assoc_vent),
    Class = case_when(
      grepl("_X", Class) ~ "unknown",
      TRUE ~ Class),
    phylum_class = paste(Phylum, Class, sep = "_")) %>%
  filter(!is.na(Phylum), !is.na(Class)) %>%
  # Step 1: collapse raw counts per sample × phylum_class
  group_by(assoc_vent_wDepth, phylum_class) %>%
  summarise(seq_sum = sum(SEQUENCE_COUNT, na.rm = TRUE), .groups = "drop") %>%
  # Step 2: get total per sample
  group_by(assoc_vent_wDepth) %>%
  mutate(
    sample_seq_total = sum(seq_sum),
    taxa_relAbun = seq_sum / sample_seq_total) %>%
  ungroup()

VENT_SITE_ORDER <- c("Background", "Vent", "Under Plume", "Peak Plume","Plume Top","30m Above", "Oxycline", "Euphotic")

axial_RelAbun <- axial_RelAbun_byVent %>%
  select(assoc_vent_wDepth, phylum_class, taxa_relAbun) %>% mutate(vent = case_when(
    grepl("ASHES", assoc_vent_wDepth) ~ "ASHES",
    grepl("International District", assoc_vent_wDepth) ~ "International District",
    grepl("Background", assoc_vent_wDepth) ~ "Background")) %>% 
  mutate(Depth_Category = case_when(
      grepl("30m", assoc_vent_wDepth) ~ "30m Above",
      grepl("oxycline", assoc_vent_wDepth) ~ "Oxycline",
      grepl("top", assoc_vent_wDepth) ~ "Plume Top",
        grepl("peak", assoc_vent_wDepth) ~ "Peak Plume",
        grepl("under", assoc_vent_wDepth) ~ "Under Plume",
        grepl("euphotic", assoc_vent_wDepth) ~ "Euphotic",
      grepl("Vent", assoc_vent_wDepth) ~ "Vent",
      grepl("Background", assoc_vent_wDepth) ~ "Background")) %>% 
 mutate(VENT_SITE_ORDER = factor(Depth_Category, levels = VENT_SITE_ORDER))

unique(axial_RelAbun_byVent$assoc_vent_wDepth)
  
  
#Find the 5th percentile cutoff
#cutoff <- quantile(taxa_abundance$mean_rel, probs = 0.05)
  
# Tag rare taxa as "Other"
#axial_RelAbun <- axial_RelAbun %>%
  mutate(phylum_class = if_else(phylum_class %in% 
                                  taxa_abundance$phylum_class[taxa_abundance$mean_rel <= cutoff],
                                "Other", phylum_class))

#Heatmap
ggplot(axial_RelAbun, aes(x = assoc_vent_wDepth, y = phylum_class, fill = taxa_relAbun)) +
  geom_tile(color = "black") +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10)) +
  labs(x = "Location", y = "Phylum_Class", fill = "Rel. Abundance")


##Stacked bar plot
ggplot(axial_RelAbun, aes(x = VENT_SITE_ORDER, y = taxa_relAbun, fill = phylum_class)) +
  geom_bar(stat = "identity", color = "black", position = "fill") +
  coord_flip() +
  facet_grid(cols = vars(vent), space = "free", scales = "free") +
  theme_minimal() +
  labs(
    title = "Relative Sequence Abundance",
    x = "", y = "")
