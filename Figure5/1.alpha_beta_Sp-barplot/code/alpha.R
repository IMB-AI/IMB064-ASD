# Call phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(DESeq2)
library(ape)


ROOT_DIR <- "set up root directory"
setwd(ROOT_DIR)

# import abundance table
ab_table <- read.csv("./Target_ASV_result.tsv", sep = "\t")

# sample list
g1_sample <- c('MG1.1','MG1.2','MG1.3','MG1.4','MG1.5','MG1.6')
g2_sample <- c('MG2.1','MG2.2','MG2.3','MG2.4','MG2.5','MG2.6')
g_sample <- c(g1_sample, g2_sample)
group_info <- c(rep("PBS", length(g1_sample)), rep("IMB064", length(g2_sample)))

# metadata
metadata <- data.frame(
  Sample = g_sample,
  Group = group_info
)
rownames(metadata) <- metadata$Sample
metadata$Group <- factor(metadata$Group, levels = c("PBS", "IMB064"))

# abundance only
ab_matrix <- ab_table[g_sample]
rownames(ab_matrix) <- ab_table$X
ab_matrix <- ab_matrix[rowSums(ab_matrix) != 0, ]

# Figure style setting
theme_set(theme_bw())

# Figure save function
save_pdf <- function(plot_expr, filename, width, height) {
  pdf(filename, width=width, height=height)
  print(plot_expr)
  dev.off()
}



# ====================== #
# beta Diversity #
# ====================== #

# abundance 넣기
abun <- otu_table(ab_matrix, taxa_are_rows = TRUE)
tr <- read.tree("tree.nwk")

# metadata 넣기
sd <- sample_data(metadata)
ps <- phyloseq(abun, sd, phy_tree(tr))
ps <- prune_samples(g_sample, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)          # 전부 0인 taxa 제거

# rarefaction
rarecurve(t(ab_matrix), step=50, cex=0.5) # rarefaction curves
ps.rarefied <- rarefy_even_depth(ps, rngseed=1, sample.size=min(sample_sums(ps)), replace=F)
# Rarefaction processed: sample size = 34327 (minimum sample depth)
# `set.seed(1)` was used to initialize repeatable random subsampling.
# 10 ASVs were removed because they are no longer present in any sample after random subsampling
sample_sums(ps)
sample_sums(ps.rarefied)


## alpha ##
alpha_df <- estimate_richness(ps.rarefied, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
alpha_df <- cbind(metadata, alpha_df)
write.table(alpha_df, "Alpha_diversity.tsv", sep="\t", row.names = FALSE)
p <- plot_richness(ps.rarefied, x="Sample", measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
save_pdf(p, "Alpha_diversity.pdf", 10, 4) # save plot 10*4
p <- plot_richness(ps.rarefied, x="Group", measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
save_pdf(p, "Alpha_diversity_group.pdf", 8, 4) # save plot 8*4
p <- p + geom_boxplot()
save_pdf(p, "Alpha_diversity_group_box.pdf", 8, 4) # save plot 8*4

# =============================================================================