## This scrip exemplified the steps followed on DNA-BSAM
setwd("U:/Analysis/R/Soil") # Change with porject directory
# Loading libraries 
set.seed(100)
packages <- c('tidyverse', 'ggpubr', 'rstatix', 'readr', 'qiime2R', 'SRS', 'vegan', 'microbiome',
              'phyloseq', 'plyr', 'dplyr', 'pairwiseAdonis', 'cowplot', 'ggvenn', 'microbial', 'UpSetR',
              'mltools', 'data.table')
sapply(packages, require, character.only = TRUE)

#Importing metadata files
metadata_16S <- read_q2metadata("metadata-16S.txt") # Change with the metadata file name
metadata_16S$Treatment <- as.character(metadata_16S$Treatment)
metadata_ITS <- read_q2metadata("metadata-ITS.txt")
metadata_ITS$Treatment <- as.character(metadata_ITS$Treatment)

#Importing numbers-feautres and shannon qza files
numbers_16S <- read_qza("observed_features-16S.qza")
numbers_16S <- numbers_16S$data
numbers_16S <- cbind(numbers_16S, Treatment = metadata_16S$Treatment)

numbers_ITS <- read_qza("observed_features-ITS.qza")
numbers_ITS <- numbers_ITS$data
numbers_ITS <- cbind(numbers_ITS, Treatment = metadata_ITS$Treatment)

#Adding a barcode variable before merging diversity data
numbers_16S$Barcode <- "16S"
numbers_ITS$Barcode <- "ITS"
numbers_ASVs <- rbind(numbers_16S, numbers_ITS)

#Stats
#Testing the effect of Treatment using Kruskal-Wallis
KW.numbers.Treatment <- numbers_ASVs %>%
  group_by(Barcode) %>%
  kruskal_test(observed_features ~ Treatment) %>%
  add_significance()

#Calculating the effect size
KW.numbers.Treatment.E2 <- numbers_ASVs %>%
  group_by(Barcode) %>%
  kruskal_effsize(observed_features ~ Treatment)

#Merging both diversity results into one object
KW.numbers.Treatment <- cbind(KW.numbers.Treatment,
                              KW.numbers.Treatment.E2[(colnames(KW.numbers.Treatment.E2) %in% c("effsize", "magnitude"))])

#Pairwise tests using Wilcoxon with Benjamini Hochberg adjusted p value 
pairwise.numbers.Treatment <- numbers_ASVs %>%
  group_by(Barcode) %>%
  wilcox_test(observed_features ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() 

#Exort statistical result 
dir.create("Stats")
write.csv(KW.numbers.Treatment, file = "Stats/KW.numbers.Treatment.csv", row.names = F)
write.csv(pairwise.numbers.Treatment, file = "Stats/pairwise.numbers.Treatment.csv", row.names = F)

#Creating box-plot of diversity results
#Adding xy positions to represent significant pairwise result with brackets
pairwise.numbers.Treatment <- pairwise.numbers.Treatment %>% add_xy_position(x = "Treatment") 
pairwise.numbers.Treatment$p.adj <- round(pairwise.numbers.Treatment$p.adj, 3)

#Calculating media value for geom_hline 
numbers_mean_o <- numbers_ASVs %>%
  group_by(Barcode) %>%
  dplyr::summarize(mean_val = mean(observed_features))

ggboxplot(numbers_ASVs, x = "Treatment", y = "observed_features", ylab = "Number of Taxa",
          scales = "free", color = "Treatment", palette = "jco", legend = "none",
          add = "jitter", short.panel.labs = TRUE, ggtheme = theme_light()) + 
  rotate_x_text(angle = 45) + ylim(0, 2250) +
  geom_hline(data= numbers_mean_o, aes(yintercept=mean_val), linetype = 2) + # Add horizontal line at base mean
  stat_compare_means(label.x = 1.2) +
  stat_pvalue_manual(pairwise.numbers.Treatment, hide.ns = TRUE, step.group.by = "Barcode",
                     bracket.nudge.y = -50, label = "p.adj", step.increase = 0.01, size = 2.2, colour = "group1") +
  facet_wrap(~ Barcode, scales = "free")

#BETA-DIVERSITY
#Import taxonomy and feature tables for creating a phyloseq object
taxonomy_16S <- read_qza("taxonomy-16S.qza")
taxonomy_ITS <- read_qza("taxonomy-ITS.qza")

taxtable_16S <- taxonomy_16S$data %>% as.tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version
taxtable_ITS <- taxonomy_ITS$data %>% as.tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) #convert the table into a tabular split version

table_16S <- read_qza("table-16S.qza")
table_16S <- as.data.frame(table_16S$data)
ps_16S <- phyloseq(otu_table(table_16S, taxa_are_rows = T),
                       tax_table(as.data.frame(taxtable_16S) %>% dplyr::select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
                       sample_data(metadata_16S %>% as.data.frame() %>% column_to_rownames("SampleID")))

table_ITS <- read_qza("table-ITS.qza")
table_ITS <- as.data.frame(table_ITS$data)
ps_ITS <- phyloseq(otu_table(table_ITS, taxa_are_rows = T),
                       tax_table(as.data.frame(taxtable_ITS) %>% dplyr::select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
                       sample_data(metadata_ITS %>% as.data.frame() %>% column_to_rownames("SampleID")))

#Calculating Jacard distances
ujaccard_16S <- phyloseq::distance(ps_16S, "jaccard", binary = TRUE)
ujaccard_ITS <- phyloseq::distance(ps_ITS, "jaccard", binary = TRUE)
wjaccard_16S <- phyloseq::distance(ps_16S, "jaccard", binary = FALSE)
wjaccard_ITS <- phyloseq::distance(ps_ITS, "jaccard", binary = FALSE)


u.permanova_16S_Treatment <- adonis(ujaccard_16S ~ Treatment, data = metadata_16S)
w.permanova_16S_Treatment <- adonis(wjaccard_16S ~ Treatment, data = metadata_16S)
u.permanova_ITS_Treatment <- adonis(ujaccard_ITS ~ Treatment, data = metadata_ITS)
w.permanova_ITS_Treatment <- adonis(wjaccard_ITS ~ Treatment, data = metadata_ITS)
permanova <- rbind(u.permanova_16S_Treatment$aov.tab[c("R2","Pr(>F)")],
                   w.permanova_16S_Treatment$aov.tab[c("R2","Pr(>F)")],
                   u.permanova_ITS_Treatment$aov.tab[c("R2","Pr(>F)")],
                   w.permanova_ITS_Treatment$aov.tab[c("R2","Pr(>F)")])

pairwiseAdonis_ujaccard_16S <- pairwise.adonis(ujaccard_16S, metadata_16S$Treatment, p.adjust.m = "BH")
pairwiseAdonis_ujaccard_16S$Distance <- "Types - 16S"
pairwiseAdonis_wjaccard_16S <- pairwise.adonis(wjaccard_16S, metadata_16S$Treatment, p.adjust.m = "BH")
pairwiseAdonis_wjaccard_16S$Distance <- "Abundances - 16S"
pairwiseAdonis_ujaccard_ITS <- pairwise.adonis(ujaccard_ITS, metadata_ITS$Treatment, p.adjust.m = "BH")
pairwiseAdonis_ujaccard_ITS$Distance <- "Types - ITS"
pairwiseAdonis_wjaccard_ITS <- pairwise.adonis(wjaccard_ITS, metadata_ITS$Treatment, p.adjust.m = "BH")
pairwiseAdonis_wjaccard_ITS$Distance <- "Abundances - ITS"
pairwiseAdonis <- rbind(pairwiseAdonis_ujaccard_16S,
                        pairwiseAdonis_wjaccard_16S,
                        pairwiseAdonis_ujaccard_ITS,
                        pairwiseAdonis_wjaccard_ITS)

write.csv(permanova, "Stats/permanova.csv") 
write.csv(pairwiseAdonis, "Stats/pairwiseAdonis.csv", row.names = F)

#Ordinating Jaccard distance matrices
PCoA_ujaccard_16S <- ape::pcoa(ujaccard_16S)
PCoA_ujaccard_ITS <- ape::pcoa(ujaccard_ITS)
PCoA_wjaccard_16S <- ape::pcoa(wjaccard_16S)
PCoA_wjaccard_ITS <- ape::pcoa(wjaccard_ITS)

#Creating PCoA
plot_ujaccard_16S <- plot_ordination(ps_16S, PCoA_ujaccard_16S, color = "Treatment",
                                     title = "Types - 16S") +
  stat_ellipse(type = "norm", linetype = 2, aes(group = Treatment)) +
  stat_ellipse(type = "t", aes(group = Treatment)) +
  scale_colour_brewer(type="qual", palette="Set1")

plot_ujaccard_ITS <- plot_ordination(ps_ITS, PCoA_ujaccard_ITS, color = "Treatment",
                                     title = "Types - ITS") +
  stat_ellipse(type = "norm", linetype = 2, aes(group = Treatment)) +
  stat_ellipse(type = "t", aes(group = Treatment)) +
  scale_colour_brewer(type="qual", palette="Set1")

plot_wjaccard_16S <- plot_ordination(ps_16S, PCoA_wjaccard_16S, color = "Treatment",
                                     title = "Abundance - 16S") +
  stat_ellipse(type = "norm", linetype = 2, aes(group = Treatment)) +
  stat_ellipse(type = "t", aes(group = Treatment)) +
  scale_colour_brewer(type="qual", palette="Set1")

plot_wjaccard_ITS <- plot_ordination(ps_ITS, PCoA_wjaccard_ITS, color = "Treatment",
                                     title = "Abundance - ITS") +
  stat_ellipse(type = "norm", linetype = 2, aes(group = Treatment)) +
  stat_ellipse(type = "t", aes(group = Treatment)) +
  scale_colour_brewer(type="qual", palette="Set1")

#Creating an object conatining the four plots without legends
grouped_plots <- plot_grid(plot_ujaccard_16S + theme(legend.position="none"),
                           plot_ujaccard_ITS + theme(legend.position="none"),
                           plot_wjaccard_16S + theme(legend.position="none"),
                           plot_wjaccard_ITS + theme(legend.position="none"),
                           ncol = 2, align = "hv")

#Substracting the lengend from the first plot
legend <- cowplot::get_legend(plot_ujaccard_16S)

#Combining plots object and legend
plot_grid(grouped_plots, legend, rel_widths = c(3, .4))
