# Load necessary libraries
library(data.table) #fread
library("SNPRelate")
library("ggplot2")
library(car)        # ANOVA and mixed models for factorial experiments
library(nlme)       # gls
library(emmeans)     # Estimated marginal means (least-squares means)
library(effects)     # Effect plots for regression models

dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP_25_PCA_karyotyping"
setwd(dir)

# Find region associating with WL
gemma_results_NP25 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_gwas/NP25_gwas_WL_mm80/gemma_lmm_results.assoc.txt")
gwas_wl = subset(gemma_results_NP25, p_wald < bonferroni_threshold)
pos_min = min(gwas_wl$ps)
pos_max = max(gwas_wl$ps)
pos_min
pos_max

# Define input and output files
vcfFile <- "P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr2_phased_regionWL.vcf.gz"
gdsFile <- "P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr2_phased_regionWL.gds"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcfFile, gdsFile, method = "biallelic.only")

# Open GDS file for analysis
genofile <- snpgdsOpen(gdsFile)

# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile, autosome.only = TRUE)

# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".BarSW.filtered.sorted.bam", "", sNames)

popInfo <- data.frame(
  V1 = c("NP25")
)

# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

# divide in three clusters based on PC1
plot_data$karyotype <- cut(plot_data$PC1, breaks = 3, labels = c("type1", "type2", "type3"))
plot_data$IID <- ccm_pca$sample.id

# PCA Plot Colored by Population
ggplot(plot_data, aes(x = PC1, y = PC2, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Population", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# PC2 and PC3
ggplot(plot_data, aes(x = PC2, y = PC3, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Population", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# IID lists of different karyotypes
iid_lists <- split(plot_data$IID, plot_data$karyotype)
iid_lists$type1
iid_lists$type2
iid_lists$type3

write.table(iid_lists$type3, "homozygotesSW.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# test hardy weinberg equilibrium
aa = length(iid_lists$type1)
aA = length(iid_lists$type2)
AA = length(iid_lists$type3)
N=192
freq_aa = (aa*2+aA)/N*2 
freq_AA = (AA*2+aA)/N*2 

exp_aa = freq_a^2*N # 9.0 expected, 11 found
exp_aA = 2*freq_a*freq_A*N # 65 expected, 61 found
exp_AA = freq_A^2*N # 118.0 expected, 120 found
(83/384)^2+2*(83/384)*(301/384) #0.3855726
(301/384)^2 #0.6144274

# Plot relMRWS of different clusters
phenotype <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_gwas/phenotype.txt")

plot_data$relMRWS <- phenotype$relMRWS[match(plot_data$IID, phenotype$IID)]

ggplot(plot_data, aes(x = karyotype, y = relMRWS, fill = karyotype)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Cluster", y = "relMRWS", title = "relMRWS per PCA Cluster")

fit <- lm(relMRWS ~ karyotype, data=plot_data)
summary(fit)
Anova(fit)
shapiro.test(residuals(fit))  # W > 0.9, so no deviation from normality
ncvTest(fit) #the variances do significantly diverge from homogeneous
spreadLevelPlot(fit)

# so use gls to allow different variances
fit2 <- gls(relMRWS ~ karyotype,weights = varIdent(form=~1|karyotype), data=plot_data)
summary(fit2)
Anova(fit2)
plot(allEffects(fit2), lty=0)

pairwise_results <- pairs(emmeans(fit2, ~ karyotype), adjust = "tukey")
summary(pairwise_results)

# Plot EL per type
plot_data$EL <- phenotype$EL[match(plot_data$IID, phenotype$IID)]

ggplot(plot_data, aes(x = karyotype, y = EL, fill = karyotype)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Cluster", y = "EL", title = "EL per PCA Cluster")

fit2 <- gls(EL ~ karyotype,weights = varIdent(form=~1|karyotype), data=plot_data)
summary(fit2)
Anova(fit2)
plot(allEffects(fit2), lty=0)

pairwise_results <- pairs(emmeans(fit2, ~ karyotype), adjust = "tukey")
summary(pairwise_results)


# Plot time5trials per type
plot_data$time5trials <- phenotype$time5trials[match(plot_data$IID, phenotype$IID)]

ggplot(plot_data, aes(x = karyotype, y = time5trials, fill = karyotype)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(x = "Cluster", y = "time5trials", title = "time5trials per PCA Cluster")

fit2 <- gls(time5trials ~ karyotype,weights = varIdent(form=~1|karyotype), data=plot_data)
summary(fit2)
Anova(fit2)
plot(allEffects(fit2), lty=0)

pairwise_results <- pairs(emmeans(fit2, ~ karyotype), adjust = "tukey")
summary(pairwise_results)
