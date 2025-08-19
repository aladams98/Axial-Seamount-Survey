library("phyloseq")
library("dplyr")
library("tidyverse")
library("ggplot2")

# 1. Three tables are needed: OTU (or ASV), taxonomy and samples
ASV <- read.delim("C:/Users/alexis98adams/Documents/R/Axial2023/datainput/deepsea-18s-merged-asv-table-05022025.tsv", header = TRUE, sep = "\t", skip = 1)
tax <- read.delim("C:/Users/alexis98adams/Documents/R/Axial2023/datainput/taxonomy.tsv", header = TRUE, sep = "\t")
metadata_input <- read.csv("C:/Users/alexis98adams/Documents/R/Axial2023/datainput/axial_compiledmetadata.csv", header = TRUE)

#Phyloseq objects need to have row.names. Define row names from otu column.
#Then trransform into matrixes. Metadata can be left as DF. 
asv_mat <- ASV %>%  tibble::column_to_rownames("X.OTU.ID") %>%  as.matrix()
tax_mat <- tax %>% tibble::column_to_rownames("Feature.ID") %>% as.matrix()


metadata_input <- metadata_input[1:90, ]

samples_df <- metadata_input %>% tibble::column_to_rownames("NEW_FASTQ_NAME")

#Transform into phyloseq objects
ASV_phy = otu_table(asv_mat, taxa_are_rows = TRUE)
tax_phy = tax_table(tax_mat)
samples = sample_data(samples_df)

Axial_phy <- phyloseq(ASV_phy, tax_phy, samples)

#Heatmap
plot_heatmap(Axial_phy, method = "NMDS", distance = "bray")
