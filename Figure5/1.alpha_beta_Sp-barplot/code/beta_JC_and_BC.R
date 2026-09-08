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


## beta - bray ##
bray_dist <- phyloseq::distance(ps.rarefied, method="bray")

set.seed(42)
a <- adonis2(bray_dist ~ sample_data(ps.rarefied)$Group, by = "terms")
a
a.res <- c(a$`R2`[1], a$`F`[1], a$`Pr(>F)`[1])
lab <- sprintf("**PERMANOVA**:   R²=%.3f,  F=%.3f,  *p*=%s",
               round(a.res[1], 3), round(a.res[2], 3),
               formatC(a.res[3], format="f", digits=3))

ord_pcoa <- ordinate(ps.rarefied, method="PCoA", distance=bray_dist)
p <- plot_ordination(ps.rarefied, ord_pcoa, color="Group") + geom_point(size=3.5)
p <- p + coord_cartesian(xlim = c(-0.4, 0.3), ylim = c(-0.5, 0.2))
p <- p + scale_color_manual(values = c("IMB064" = "#1f78b4", "PBS" = "#e31a1c"))
p <- p + labs(color="Group", title="PCoA (Bray-Curtis)")
p <- p + ggtext::geom_richtext(label.color = NA, size = 3.5, fill=NA,
                               color = "grey40", inherit.aes = FALSE,
                               hjust = 0, vjust =0, x = -Inf, y = -Inf,
                               label= lab)
print(p)
save_pdf(p, "Beta_diversity_PCoA_bray.pdf", 6, 5)



## beta - Jaccard ##
bray_dist <- phyloseq::distance(ps.rarefied, method="jaccard", binary = TRUE)

set.seed(42)
a <- adonis2(bray_dist ~ sample_data(ps.rarefied)$Group, by = "terms")
a
a.res <- c(a$`R2`[1], a$`F`[1], a$`Pr(>F)`[1])
lab <- sprintf("**PERMANOVA**:   R²=%.3f,  F=%.3f,  *p*=%s",
               round(a.res[1], 3), round(a.res[2], 3),
               formatC(a.res[3], format="f", digits=3))

ord_pcoa <- ordinate(ps.rarefied, method="PCoA", distance=bray_dist)
p <- plot_ordination(ps.rarefied, ord_pcoa, color="Group") + geom_point(size=3.5)
p <- p + coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.25, 0.4)) #jaccard
p <- p + scale_color_manual(values = c("IMB064" = "#1f78b4", "PBS" = "#e31a1c"))
p <- p + labs(color="Group", title="PCoA (Jaccard)")
p <- p + ggtext::geom_richtext(label.color = NA, size = 3.5, fill=NA,
                               color = "grey40", inherit.aes = FALSE,
                               hjust = 0, vjust =0, x = -Inf, y = -Inf,
                               label= lab)
print(p)
save_pdf(p, "Beta_diversity_PCoA_jaccard.pdf", 6, 5)
# =============================================================================