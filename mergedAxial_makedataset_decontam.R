#Setup Axial Dataset and merge with metadata

library(tidyverse)

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)
library(phyloseq)

# Change working directory to locate file
setwd("/scratch/group/hu-lab/deepsea-18S-merged-outputs")

list.files("/scratch/group/hu-lab/deepsea-18S-merged-outputs/merged-taxonomy-0.9_0.8_PR2-5.0.0/taxonomy-merged-05022025")


#Import file
deepsea_18S_merged_asv <- read.delim("/scratch/group/hu-lab/deepsea-18S-merged-outputs/deepsea-18s-merged-asv-table-05022025.tsv", header = TRUE, sep = "\t", skip=1)
deepsea_18S_PR2_tax <- read.delim("/scratch/group/hu-lab/deepsea-18S-merged-outputs/merged-taxonomy-0.9_0.8_PR2-5.0.0/taxonomy-merged-05022025/taxonomy-05022025.tsv", sep = "\t")


deepsea_18S_asv_wtax_NO_FILTER <- deepsea_18S_merged_asv %>% 
  select(Feature.ID = 'X.OTU.ID', everything()) %>% 
  pivot_longer(cols = !Feature.ID, 
               names_to = "SAMPLE", values_to = "value")%>% 
  left_join(deepsea_18S_PR2_tax, by = c("Feature.ID" = "Feature.ID")) %>% 
  separate(Taxon, c("Domain", "Supergroup",
                    "Phylum", "Class", "Order",
                    "Family", "Genus", "Species"), sep = ";", remove = FALSE) 

axial_18S_asv_wtax_NOFILTER <- subset(deepsea_18S_asv_wtax_NO_FILTER, grepl("AXIAL", SAMPLE, ignore.case = TRUE)) %>% 
  select(-Taxon) 
#71341 obs of 12 variables



#Save files to cluster
write.table(axial_18S_asv_wtax_NOFILTER, 
            "/scratch/user/alexis98adams/deepsea-18S-datasets/axial_18S_asv_wtax_NOFILTER.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

##Decontam
axial_18S_asv_wtax_NOFILTER <- read.delim("/scratch/user/alexis98adams/deepsea-18S-datasets/axial_18S_asv_wtax_NOFILTER.tsv", header = TRUE, sep = "\t")
axial_compiledmetadata <- read.csv("/scratch/user/alexis98adams/deepsea-18S-datasets/axial_compiledmetadata.csv") 

axial_18S_compiled <- axial_18S_asv_wtax_NOFILTER %>% 
  mutate(sample_or_control = case_when(
    grepl("Control", SAMPLE, ignore.case = TRUE) ~ "control",
    TRUE ~ "Sample")) %>% 
  left_join(axial_compiledmetadata, by = c("SAMPLE" = "NEW_FASTQ_NAME"))

#Import datasets as phyloseq objects
tax_matrix <- axial_18S_compiled %>% select(Feature.ID, Domain, Supergroup, 
                                            Phylum, Class, Order, Family, Genus, Species) %>% 
  distinct() %>% column_to_rownames(var = "Feature.ID") %>% as.matrix

asv_matrix <- axial_18S_compiled %>% select(Feature.ID, SAMPLE, value) %>% 
  pivot_wider(names_from = "SAMPLE", values_fill = 0, values_from = value) %>% 
  column_to_rownames(var = "Feature.ID") %>% as.matrix

rownames(tax_matrix) <- row.names(asv_matrix)

axial_metadata_cones <- axial_18S_compiled %>% 
  select(SAMPLE, SAMPLE_TYPE, Location_Name, sample_or_control) %>% 
  distinct() %>%  column_to_rownames(var = "SAMPLE")

#Import asv and tax matrices
ASV = otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
phylo_obj <- phyloseq(ASV, TAX)

#Import metadata as sample data in phyloseq
samplenames <- sample_data(axial_metadata_cones)
ls 
#Join as phyloseq object
physeq_wnames = merge_phyloseq(phylo_obj, samplenames)

#Identify contaminant ASVs
# When "Control" appears in "Sample_or_Control column, this is a negative control"
sample_data(physeq_wnames)$is.neg <- sample_data(physeq_wnames)$sample_or_control == "control"

#ID contaminants using Prevalence information
contam_prev <- isContaminant(physeq_wnames, 
                             method="prevalence", 
                             neg="is.neg", 
                             threshold = 0.5, normalize = TRUE) 

# Report number of ASVs IDed as contamintants
table(contam_prev$contaminant)


##Remove problematic ASVs
#Subset contaminant ASVs
contams <- filter(contam_prev, contaminant == "TRUE")
list_of_contam_asvs <- as.character(row.names(contams))

taxa_contam <- as.data.frame(tax_matrix) %>% 
  rownames_to_column(var = "Feature.ID") %>% 
  filter(Feature.ID %in% list_of_contam_asvs)

#View
axial_asv_wtax_decon <- axial_18S_compiled %>% 
  filter(!(Feature.ID %in% list_of_contam_asvs)) %>% 
  filter(!(sample_or_control == "control"))

tmp_orig <- (axial_18S_compiled %>% filter(!(sample_or_control == "control")))

#Stats on lost ASVs
x <- length(unique(tmp_orig$Feature.ID)); x
#27559
y <- length(unique(axial_asv_wtax_decon$Feature.ID)); y #26911
100*((y-x)/x) #2.35% of ASVs lost
a <- sum(tmp_orig$value);a #8060496 --> 8.1 million
b <- sum(axial_asv_wtax_decon$value);b #6944159
100*((b-a)/a) #Lost 13.85% of sequences from whole dataset


## Subsample to clean ASVs
axial_asv_wtax_wstats_NOFILTER <- axial_18S_compiled %>% 
  mutate(DECONTAM = case_when(
    Feature.ID %in% list_of_contam_asvs ~ "FAIL",
    TRUE ~ "PASS"
  ))

axial_18S_asv_FinalProduct <- axial_asv_wtax_wstats_NOFILTER %>% 
  filter(sample_or_control == "Sample") %>% 
  filter(DECONTAM == "PASS") %>% 
  select(SAMPLE, Feature.ID, value, Domain, Supergroup, Phylum, 
         Class, Order, Family, Genus, Species, Location_Name, SAMPLE_TYPE, DEPTH_CATEGORY,
         DEPTH, DOC..um., POC..um., TDN.um., Prok_avg)

write.csv(axial_18S_asv_FinalProduct, file = "/scratch/user/alexis98adams/deepsea-18S-datasets/axial_18S_asv_FinalProduct.csv", row.names = FALSE)

