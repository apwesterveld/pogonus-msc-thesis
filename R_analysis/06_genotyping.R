# Load necessary libraries
library(data.table)  #fread
library(SNPRelate)
library(ggplot2)
library(car)         # ANOVA and mixed models for factorial experiments
library(nlme)        # gls
library(emmeans)     # Estimated marginal means (least-squares means)
library(effects)     # Effect plots for regression models
library(scales)      # for label_percent()
library(dplyr)
library(plotly)      # 3D PCA plot
library(htmlwidgets) # to save 3D plot



dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef"
setwd(dir)

gemma_results_NP25 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80_homoSW/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80_AA_bb/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25 <- fread("NP25_gwas/NP25_gwas_EL_cov/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25 <- fread("NP25_gwas/NP25_gwas_cov_time3trials/gemma_results/gemma_lmm_results.assoc.txt")

## Region for PCA ####
num_tests <- nrow(gemma_results_NP25)  # Total number of SNPs tested
bonferroni_threshold <- 0.05 / num_tests  # Genome-wide significance threshold
gwas_wl = subset(gemma_results_NP25, p_wald < bonferroni_threshold*100)
pos_min = min(gwas_wl$ps)
pos_max = max(gwas_wl$ps)

# vcf/gds files ####
# Define input and output files
# chr2
vcfFile <- "NP25_PCA_karyotyping/NP25_PCA_chr2/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr2_phased_regionWL.vcf.gz"
gdsFile <- "NP25_PCA_karyotyping/NP25_PCA_chr2/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr2_phased_regionWL.gds"
# chr3
vcfFile <- "NP25_PCA_karyotyping/NP25_PCA_chr3/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3_regionWL_phased.vcf.gz"
gdsFile <- "NP25_PCA_karyotyping/NP25_PCA_chr3/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3_regionWL_phased.gds"
# chr3_2
vcfFile <- "NP25_PCA_karyotyping/NP25_PCA_chr3_2/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3inv.vcf.gz"
gdsFile <- "NP25_PCA_karyotyping/NP25_PCA_chr3_2/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3inv.gds"
# chr2 AA chr3 bb
vcfFile <- "NP25_PCA_karyotyping/NP25_PCA_AA_bb/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3inv_AA_bb.vcf.gz"
gdsFile <- "NP25_PCA_karyotyping/NP25_PCA_AA_bb/P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit_chr3inv_AA_bb.gds"

## Convert VCF to GDS format
#snpgdsVCF2GDS(vcfFile, gdsFile, method = "biallelic.only")

# Open GDS file for analysis
genofile <- snpgdsOpen(gdsFile)
# check chr region in gds
chroms <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
positions <- read.gdsn(index.gdsn(genofile, "snp.position"))
ranges <- tapply(positions, chroms, range)
# number of samples
length(read.gdsn(index.gdsn(genofile, "sample.id")))

## PCA on highest SNP ####
# choose snp with highsest significance
min_p_idx <- which.min(gemma_results_NP25$p_wald)
target_pos <- gemma_results_NP25$ps[min_p_idx]
target_chr <- gemma_results_NP25$chr[min_p_idx]  

positions <- read.gdsn(index.gdsn(genofile, "snp.position"))
chroms <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp_idx <- which(chroms == target_chr & positions == target_pos)

# Find matching SNP in GDS
snp_idx <- which(chroms == target_chr & positions == target_pos)
# do PCA on just one SNP
ccm_pca <- snpgdsPCA(genofile, autosome.only = TRUE, snp.id = snp_idx)
snpgdsClose(genofile)

## PCA on region around highest SNP ####
# Step 1: Find the most significant SNP
min_p_idx <- which.min(gemma_results_NP25$p_wald)
target_pos <- gemma_results_NP25$ps[min_p_idx]
target_chr <- gemma_results_NP25$chr[min_p_idx]

# Step 2: Define a window around the SNP (e.g., ±100 kb)
window_size <- 500  # adjust as needed
start_pos <- target_pos - window_size
end_pos <- target_pos + window_size

# Step 3: Read positions and chromosomes from the GDS
positions <- read.gdsn(index.gdsn(genofile, "snp.position"))
chroms <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

# Step 4: Filter SNPs in that region
region_idx <- which(
  chroms == target_chr &
    positions >= start_pos &
    positions <= end_pos
)

# Optional: Check how many SNPs are selected
length(region_idx)

# Step 5: Run PCA on SNPs in this region
ccm_pca <- snpgdsPCA(genofile, snp.id = region_idx, autosome.only = FALSE)
snpgdsClose(genofile)


## PCA on all SNPs ####
# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile, autosome.only = TRUE)
snpgdsClose(genofile)

# compiling
var_pct <- data.frame(
  PC = paste0("PC", 1:5),
  Variance = (ccm_pca$varprop * 100)[1:5]
)


# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".BarSW.filtered.sorted.bam", "", sNames)

# Screeplot ####
screeplot <- ggplot(var_pct, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 16) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Screeplot",
       y = "Variance Explained",
       x = NULL)
screeplot
#ggsave(screeplot, file = "NP25_PCA_chr3/screeplot.svg", width = 3, height = 2.5)

# Plotting ####
# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

pca_1_2 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3,alpha=0.3) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal(base_size = 16) + 
  theme(legend.position = c(0.8,0.8),legend.title = element_blank())
pca_1_2

# Clustering ####
### cut() ####
# chr2
plot_data$cluster <- cut(plot_data$PC1, breaks = 3, labels = c('chr2: aa','chr2: aA','chr2: AA'))
# chr3
plot_data$cluster <- cut(plot_data$PC1, breaks = 3, labels = c('chr3:bb','chr3:bB','chr3:BB'))
plot_data$IID <- ccm_pca$sample.id

### kmeans() ####
set.seed(123)
kmeans_result <- kmeans(plot_data[, c("PC1", "PC2", "PC3")], centers = 3, nstart = 25)
plot_data$cluster <- as.factor(kmeans_result$cluster)
levels(plot_data$cluster) <- c("BB", "bB", "bb")


### plot ####
pca_1_2 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3,alpha=0.3) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal(base_size = 18) + 
  stat_ellipse(aes(color = cluster), type = "norm", level = 0.95, linetype = "solid", size = 1.5) +
  theme(legend.position = c(0.8,0.2),legend.title = element_blank())
pca_1_2

ggsave(pca_1_2, file = "NP25_PCA_karyotyping/NP25_PCA_chr3_2/PCA_PC1andPC2.svg", width = 4, height = 3.5)

# PC2 and PC3
pca_2_3 <- ggplot(plot_data, aes(x = PC2, y = PC3)) +
  geom_point(size = 2,alpha=0.3) +
  labs(x = "Principal Component 2", y = "Principal Component 3") +
  theme_minimal() +
  stat_ellipse(aes(color = cluster), type = "norm", level = 0.95, linetype = "solid", size = 1.5) +
  theme(legend.position = "none")
pca_2_3
ggsave(pca_2_3, file = "NP25_PCA_karyotyping/NP25_PCA_chr3/PCA_PC2andPC3.svg", width = 4, height = 3.5)


### 3D plot ####
pca_3d <- plot_ly(
  data = plot_data,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  type = "scatter3d",
  mode = "markers",
  color = ~cluster,  # Automatically assigns one color per cluster
  marker = list(size = 4, opacity = 0.8)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    legend = list(title = list(text = "Cluster"))
  )
pca_3d

saveWidget(pca_3d, "NP25_PCA_karyotyping/NP25_PCA_chr3_2/3D_PCA_chr3.html", selfcontained = TRUE)

# Export genotypes ####
iid_lists <- split(plot_data$IID, plot_data$cluster)

write.table(iid_lists$bb, "NP25_PCA_karyotyping/NP25_PCA_chr3_2/karyotype_1SNPchr3_bb1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(iid_lists$bB, "NP25_PCA_karyotyping/NP25_PCA_chr3_2/karyotype_1SNPchr3_bB2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(iid_lists$BB, "NP25_PCA_karyotyping/NP25_PCA_chr3_2/karyotype_1SNPchr3_BB3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(
  plot_data[, c("IID", "cluster")],
  file = "PCA_chr3_WLpeak_AAbb_7clusters.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

# HW test ####
aa = length(iid_lists$type1)
aa = 11
aA = length(iid_lists$type2)
aA = 61
AA = length(iid_lists$type3)
AA = 120
N=192
freq_a = (aa*2+aA)/(N*2) 
freq_A = (AA*2+aA)/(N*2) 

exp_aa = freq_a^2*N # 9.0 expected, 11 found
exp_aA = 2*freq_a*freq_A*N # 65 expected, 61 found
exp_AA = freq_A^2*N # 118.0 expected, 120 found

observed <- c(aa, aA, AA)
observed
expected <- c(exp_aa, exp_aA, exp_AA)
expected
chisq.test(observed, p = expected / sum(expected))

HWtable <- data.frame(
  karyotype = c("Homozygote LW", "Heterozygote", "Homozygote SW"),
  sample = c(length(iid_lists$type1), length(iid_lists$type2), length(iid_lists$type3)),
  HW_equilibrium = c(exp_aa, exp_aA, exp_AA),
  stringsAsFactors = FALSE
)
print(HWtable)
write.table(
  HWtable, 
  file = "HWtable.txt", 
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)
bb = sum(phenotypes$chr3 == "chr3:bb")
bB = sum(phenotypes$chr3 == "chr3:bB")
BB = sum(phenotypes$chr3 == "chr3:BB")
N=192
freq_b = (bb*2+bB)/(N*2) 
freq_B = (BB*2+bB)/(N*2) 

exp_bb = freq_b^2*N # 9.0 expected, 11 found
exp_bB = 2*freq_b*freq_B*N # 65 expected, 61 found
exp_BB = freq_B^2*N # 118.0 expected, 120 found

observed <- c(bb, bB, BB)
observed
expected <- c(exp_bb, exp_bB, exp_BB)
expected
chisq.test(observed, p = expected / sum(expected))

# Phenotype of clusters ####
phenotype <- fread("NP25_phenotypes/phenotype.txt")

plot_data$relMRWS <- phenotype$relMRWS[match(plot_data$IID, phenotype$IID)]

ggplot(plot_data, aes(x = cluster, y = relMRWS, fill = cluster)) +
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


# statistical testing

pairwise_results <- pairs(emmeans(fit2, ~ karyotype), adjust = "tukey")
summary(pairwise_results)

phenotypes_clean <- fread("NP_25_PCA_karyotyping/phenotypes_clean.txt")
str(phenotypes_clean)
data_long_time <- melt(data, id.vars = c("IID",'EL','Esize','WL','WEratio','WEsurfaceratio','relMRWS','sex','emergeCount','emergeCount3trials'), measure.vars = c('time1','time2','time3','time4','time5'),
                       variable.name = "trial", value.name = "time")
phenotype_long <- melt(phenotypes_clean,id.vars = c("relMRWS"), measure.vars = c('chr2_type','chr3_type'),variable.name = 'trial',value.name = 'chr_combination')


ggplot(phenotypes_clean, aes(x = chr2_type_chr3_type, y = relMRWS)) +
  geom_violin(fill = "lightgray", color = "black") +
  geom_boxplot(width = 0.1, outlier.size = 0.8) +
  labs(x = "chr2 / chr3 combination", y = "relMRWS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Genotype 1 SNP ####
# 1. Identify most significant SNP
min_idx <- which.min(gemma_results_NP25$p_wald)
chr <- gemma_results_NP25$chr[min_idx]
pos <- gemma_results_NP25$ps[min_idx]

# 2. Load SNP and sample data
alleles <- read.gdsn(index.gdsn(genofile, "snp.allele"))
chroms <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
positions <- read.gdsn(index.gdsn(genofile, "snp.position"))
samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

# 3. Match SNP index in GDS
idx <- which(chroms == chr & positions == pos)
allele_pair <- strsplit(alleles[idx], "/")[[1]]  # e.g. A/G

# 4. Get genotypes (0 = ref/ref, 1 = ref/alt, 2 = alt/alt)
geno <- snpgdsGetGeno(genofile, snp.id = idx)

# 5. Map numeric genotypes to allele strings
genostr <- c(
  paste(allele_pair[1], allele_pair[1], sep = "/"),
  paste(allele_pair[1], allele_pair[2], sep = "/"),
  paste(allele_pair[2], allele_pair[2], sep = "/")
)
geno_alleles <- ifelse(is.na(geno), NA, genostr[geno + 1])

# 6. Combine with sample names
df <- data.frame(sample = samples, genotype = geno_alleles, stringsAsFactors = FALSE)
head(df)
