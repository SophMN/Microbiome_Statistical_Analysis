##Load and install the required packages
BiocManager::install(c("microbiome", "picante", "dendextend", 
                       "selbal", "rms", "breakaway"))
BiocManager::install("HMP")
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
library(HMP)

##Load the data
ps <- readRDS("output/ps.rds")
ps
ps_rel_abund <- readRDS("output/ps_rel_abund_giloteaux_2016.rds")
ps_rel_abund
ps_rare <- readRDS("output/ps_rare.rds")
ps_rare
alpha_diversity <- readRDS("output/alpha_diversity.rds")
alpha_diversity
ps_clr <- readRDS("output/ps_clr.rds")

##Access the OTU table, sample metadata and taxonomy table files
otu_table(ps)
head(otu_table(ps))
dim(otu_table(ps))
colnames(otu_table(ps))
sample_data(ps)
colnames(sample_data(ps))
dim(sample_data(ps))
tax_table(ps)
dim(tax_table(ps))
colnames(tax_table(ps))

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
saveRDS(ps, file = "output/ps.rds")

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

#Save the output after transformation
saveRDS(ps_rel_abund, file = "output/ps_rel_abund_giloteaux_2016.rds")

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
saveRDS(ps_phylum, "output/ps_phylum_giloteaux_2016.rds")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]

#Melt and plot
phyloseq::psmelt(ps_phylum) %>% 
  ggplot(data = ., aes(x = Status, y = Abundance)) +
  geom_boxplot() +
  geom_jitter(aes(color = OTU, height = 0, width = .2)) +
  labs(x = "", y = "Abundance") +
  facet_wrap(~OTU, scales = "free")

##Evaluate whether taxa frequencies observed in both groups are equal
#Subset the groups
controls <- phyloseq::subset_samples(ps_phylum, Status == "Control")
cf <- phyloseq::subset_samples(ps_phylum, Status == "Chronic Fatigue")
cf

#Convert the OTU tables of the groups into data frames for statistical testing
controls_otu <- data.frame(phyloseq::otu_table(controls))
cf_otu <- data.frame(phyloseq::otu_table(cf))
class(controls_otu)
View(controls_otu)

#Group the rare phyla together to improve testing
controls_otu <- controls_otu %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes
         + Verrucomicrobia + Fusobacteria) %>%
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes
                , -Verrucomicrobia, -Fusobacteria)

cf_otu <- cf_otu %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes 
         + Verrucomicrobia + Fusobacteria) %>% 
  dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes 
                , -Verrucomicrobia, -Fusobacteria)

#Perform the Dirichlet-Multinomial test comparison across several samples
group_data <- list(controls_otu, cf_otu)
group_data
class(group_data)
xdc <- HMP::Xdc.sevsample(group_data)
xdc

##Hierarchical clustering
#Extract the OTU table from the transformed phyloseq object with relative abundances
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))

#Transform the OTU table
ps_rel_otu <- t(ps_rel_otu)
ps_rel_otu

#Compute the Bray-Curtis dissimilarity
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
bc_dist
as.matrix(bc_dist)[1:5, 1:5]
class(bc_dist)

#Save as a dendogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
ward

#Convert the sample metadata into a data frame
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
glimpse(meta)

#Provide the colour codes
unique(meta$Status)
colorCode <- c(Control = "red", `Chronic Fatigue` = "blue")
labels_colors(ward) <- colorCode[meta$Status][order.dendrogram(ward)]
plot(ward)

##Alpha diversity
#Examine the correlation between the observed OTUs and total read count
ggplot(data = data.frame("total_reads" = phyloseq::sample_sums(ps),
                         "observed" = phyloseq::estimate_richness
                         (ps, measures = "Observed")[, 1]), aes(x = total_reads, y = observed)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Total reads", y = "Observed richness")

#Subsample reads from each sample and rarefy to a consistent sequencing depth
ps_rare <- phyloseq::rarefy_even_depth(ps, rngseed = 123, replace = FALSE)
ps_rare
saveRDS(ps_rare, file = "output/ps_rare.rds")
head(phyloseq::sample_sums(ps_rare))

#Examine the correlation between observed OTUs and total reads after rarefaction
ggplot(data = data.frame("total_reads" = phyloseq::sample_sums(ps_rare),
                         "observed" = phyloseq::estimate_richness
                         (ps_rare, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Total reads", y = "Observed richness")

#Generate a data.frame with alpha diversity metrics
alpha_diversity <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  "Status" = phyloseq::sample_data(ps_rare)$Status
)
head(alpha_diversity)
View(alpha_diversity)
saveRDS(alpha_diversity, file = "output/alpha_diversity.rds")

#Plot alpha diversity metrics
alpha_diversity %>% 
  ggplot(aes(x = Status, y = Observed, color = Status, fill = Status)) +
  geom_boxplot(color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("Control" = "red", "Chronic Fatigue" = "blue")) +
  labs(x = "Status", y = "Observed richness")

alpha_diversity %>% 
  ggplot(aes(x = Status, y = Shannon, color = Status, fill = Status)) +
  geom_boxplot(color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("Control" = "red", "Chronic Fatigue" = "blue")) +
  labs(x = "Status", y = "Shannon's diversity index")

#Calculate the median of the alpha diversity metrics
alpha_diversity %>% 
  group_by(Status) %>% 
  summarise(median_observed = median(Observed),
            median_shannon = median(Shannon))

#Statistical testing
wilcox.test(Observed ~ Status, data = alpha_diversity, exact = FALSE, conf.int = TRUE)
wilcox.test(Shannon ~ Status, data = alpha_diversity, conf.int = TRUE)

##Beta diversity
#Centred log ratio (CLR) transformation of abundance data
ps_clr <- microbiome::transform(ps, "clr")
ps_clr
saveRDS(ps_clr, file = "output/ps_clr.rds")
nrow(otu_table(ps_clr))

#Examine the OTU table after CLR transformation
phyloseq::otu_table(ps_clr)[1:5, 1:5]

#Principal component analysis
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

#Examine the eigen values and %proportion of variance explained
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig, function(x) x / sum(ord_clr$CA$eig))

#Scale the axes and plot the ordination
clr_1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr_2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ps, ord_clr, type = "samples", color = "Status") +
  geom_point(size = 2) +
  coord_fixed(clr_2 / clr_1) +
  stat_ellipse(aes(group = Status), linetype = 2)

#PERMANOVA
#Generate a distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean")

#ADONIS test
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Status)

#Perform a dispersion test since the ADONIS test can be confounded by differences in dispersion
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Status)
dispr
permutest(dispr)
plot(dispr, main = "Ordination and Dispersion Labels: Aitchison Distance", sub = "")

#Beta diversity analysis using Bray-Curtis dissimilarity
#
bray_dist <- phyloseq::distance(ps_rel_abund, method = "bray")
bray_dist
as.matrix(bray_dist)[1:5, 1:5]

#Ordination
bray_ord <- phyloseq::ordinate(ps_rel_abund, method = "PCoA", distance = bray_dist)
bray_ord

#Generate PCoA plot
plot_ordination(ps_rel_abund, bray_ord, color = "Status") +
  stat_ellipse(aes(group = Status), linetype = 2)

#ADONIS test
vegan::adonis2(bray_dist ~ phyloseq::sample_data(ps_rel_abund)$Status)

#Differential abundance testing
##Wilcoxon rank sum test
#Generate a data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox
colnames(ps_wilcox)
ps_wilcox$Status <- phyloseq::sample_data(ps_clr)$Status

#Define functions to pass to map
wilcox_model <- function(df){wilcox.test(abund ~ Status, data = df)}
wilcox_pval <- function(df){wilcox.test(abund ~ Status, data = df)$p.value}

#Create nested data frames by OTU and loop over each using map
wilcox_results <- ps_wilcox %>% 
  gather(key = OTU, value = abund, -Status) %>% 
  group_by(OTU) %>% 
  nest() %>% 
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))

#Show results
head(wilcox_results)
head(wilcox_results$data[[1]])
wilcox_results$wilcox_test[[1]]

#Unnest
wilcox_results <- wilcox_results %>% 
  dplyr::select(OTU, p_value) %>% 
  unnest()
head(wilcox_results)

#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% 
  rownames_to_column(var = "OTU")

#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>% 
  full_join(taxa_info) %>% 
  arrange(p_value) %>% 
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>% 
  filter(BH_FDR < 0.05) %>% 
  dplyr::select(OTU, p_value, BH_FDR, everything())

#Print results
print.data.frame(wilcox_results)

#Save results
write.csv(wilcox_results, file = "output/wilcox_results.csv")





