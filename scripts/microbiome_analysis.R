##Load and install the required packages
BiocManager::install(c("microbiome", "picante", "dendextend", 
                       "selbal", "rms", "breakaway"))
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(vegan)
library(ALDEx2)
library(metagenomeSeq)
library(microbiome)
library(picante)
library(dendextend)
library(rms)
library(breakaway)
library(selbal)

##Load the phyloseq object
ps <- readRDS("data/ps_giloteaux_2016.rds")
ps

##Access the OTU table, sample metadata and taxonomy table files
otu_table(ps)
dim(otu_table(ps))
sample_data(ps)
dim(sample_data(ps))
tax_table(ps)
dim(tax_table(ps))
colnames(sample_data(ps))

##Sort the samples by total read count
sort(phyloseq::sample_sums(ps))

##Perform sample filtering: remove samples with less than 5k reads
ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 5000)
ps

##Remove OTUs seen in those samples
ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
ps

##Assign a new sample metadata field
phyloseq::sample_data(ps)$Status <- ifelse(phyloseq::sample_data(ps)$Subject == "Patient",
                                           "Chronic Fatigue", "Control")
sample_data(ps)$Status
class(sample_data(ps)$Status)

##Change the data type of Status to factor
phyloseq::sample_data(ps)$Status <- factor(phyloseq::sample_data(ps)$Status, 
                                           levels = c("Control", "Chronic Fatigue"))
class(sample_data(ps)$Status)

ps %>% 
  sample_data %>% 
  dplyr::count(Status)
table(sample_data(ps)$Status)

##Visualise the relative abundance
#Get the count of phyla
table(phyloseq::tax_table(ps)[, "Phylum"])

#Get the count of orders
table(phyloseq::tax_table(ps)[, "Order"])

#Get the count of genera
table(phyloseq::tax_table(ps)[, "Genus"])

#Convert the counts into relative abundance
ps_rel_abund = phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})
ps_rel_abund

#Inspect the count/abundance data before transformation
phyloseq::otu_table(ps)[1:5, 1:5]

#Inspect the abundance data after transformation
phyloseq::otu_table(ps_rel_abund)[1:5, 1:5]

#Visualise the relative abundance
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative abundance") +
  facet_wrap(~Status, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Agglomerate to the phylum level and rename
ps_phylum <- phyloseq::tax_glom(ps, "Phylum")
ps_phylum
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
phyloseq::psmelt(ps_phylum) %>% 
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(aes(color = OTU, height = 0, width = .2)) +
  labs(x = "", y = "Abundance") +
  facet_wrap(~OTU, scales = "free")



















