####
# put 1 or 2 manhatten plots together with Fst plot
# manhatten plots from lmm (using p_wald) or bslmm (bayesian, using PIP)

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)   # Arranging multiple grid-based plots
library(cowplot) # can align multi plots with plot_grid()

dir <- "..."
setwd(dir)

# Chromosome mapping ####
chromosome_mapping <- data.frame(
  chromosome = 1:10,  # Assign chromosome numbers 1 to 10
  scaffold = c('CHR1', 'CHR2', 'CHR3', 'CHR4', 'CHR5', 'CHR6', 'CHR7', 'CHR8', 'CHR9', 'CHR10'),
  chromStarts = c(
    1,
    103436458,
    213016840,
    271369969,
    320927038,
    386693873,
    460225492,
    516670785,
    576006880,
    609931117
  ),
  chromMid = c(
    47218228.5,
    153726648,
    237693403.5,
    291648502.5,
    349310454.5,
    418959681.5,
    483948137.5,
    541838831.5,
    588468997.5,
    623145269
  )) 

chromosome_mapping <- chromosome_mapping %>%
  mutate(scaffold = as.character(scaffold))
chromosome_mapping <- chromosome_mapping %>%
  mutate(chromosome = as.factor(chromosome))

# Fst plot####
FstSpain <- fread("Fst scans/new_genome_stats_Spain.tsv")
FstSpain$FstWC_S_HUE_T_S_COT_S_egg[FstSpain$FstWC_S_HUE_T_S_COT_S_egg == "None"] <- NA
FstSpain$FstWC_S_HUE_T_S_COT_S_egg <- as.numeric(FstSpain$FstWC_S_HUE_T_S_COT_S_egg)
FstSpain$FstWC_S_HUE_T_S_COT_S_egg[FstSpain$FstWC_S_HUE_T_S_COT_S_egg<0] <- 0

Fst <- ggplot(FstSpain,aes(x = chromPos,y = FstWC_S_HUE_T_S_COT_S_egg)) +
  geom_point(aes(color = factor(chromosome))) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]  
    }
  )+
  theme_minimal() +
  labs(x = "Chromosome",
       y = expression(F[st])) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))


# Manhatten p_wald 1 ####
gemma_results_NP25_1 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25_1 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80_homoSW/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25_1 <- fread("NP25_gwas/NP25_gwas_cov_time3trials/gemma_results/gemma_lmm_results.assoc.txt")
gemma_results_NP25_1 <- fread("NP25_gwas/NP25_gwas_EL_cov/gemma_results/gemma_lmm_results.assoc.txt")

gemma_results_NP25_1$chr <- as.factor(gemma_results_NP25_1$chr)  # Make sure chromosome is a factor
gemma_results_NP25_1$p_wald <- as.numeric(gemma_results_NP25_1$p_wald)
gemma_results_NP25_1$ps <- as.numeric(gemma_results_NP25_1$ps)

gemma_results_NP25_1 <- gemma_results_NP25_1 %>%
  left_join(chromosome_mapping, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# bonferroni treshold corrects for number of tests (conservative)
num_tests_1 <- nrow(gemma_results_NP25_1)  # Total number of SNPs tested
bonferroni_threshold_1 <- 0.05 / num_tests_1  # Genome-wide significance threshold
bonferroni_threshold_1_log <- -log10(bonferroni_threshold_1)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_1 <- length(unique(gemma_results_NP25_1$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

# Plot1
manhatten_pwald_1 <- ggplot(gemma_results_NP25_1, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald > bonferroni_threshold_1),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald < bonferroni_threshold_1),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_1_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_1_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  theme_minimal() +
  labs(title = "%MRWS (n=192, sex as covariate)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
gwas1

# Combined plot ####

manhattenANDfst <- plot_grid(
  manhatten_pwald_1,
  Fst,
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(1, 0.4)
)
manhattenANDfst
output_file <- "NP25_gwas/NP25_gwas_relMRWS_cov/gwas_relMRWS_cov.svg"

if (!file.exists(output_file) || readline("File exists. Overwrite? [y/n]: ") == "y") {
  ggsave(manhattenANDfst, file = output_file, width = 9, height = 4)
}

# Manhatten p_wald 2 ####
gemma_results_NP25_2 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80_Chr2InvAA_chr3SnpBB/gemma_results/gemma_lmm_results.assoc.txt")

gemma_results_NP25_2$chr <- as.factor(gemma_results_NP25_2$chr)  # Make sure chromosome is a factor
gemma_results_NP25_2$p_wald <- as.numeric(gemma_results_NP25_2$p_wald)
gemma_results_NP25_2$ps <- as.numeric(gemma_results_NP25_2$ps)

gemma_results_NP25_2 <- gemma_results_NP25_2 %>%
  left_join(chromosome_mapping, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# bonferroni treshold corrects for number of tests (conservative)
num_tests_2 <- nrow(gemma_results_NP25_2)  # Total number of SNPs tested
bonferroni_threshold_2 <- 0.05 / num_tests_2  # Genome-wide significance threshold
bonferroni_threshold_2_log <- -log10(bonferroni_threshold_2)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_2 <- length(unique(gemma_results_NP25_2$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_2)

# Plot2
manhatten_pwald_2 <- ggplot(gemma_results_NP25_2, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald > bonferroni_threshold_2),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald < bonferroni_threshold_2),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_2_log, color = "red", linetype = "dashed", size = 0.5) +  # Bonferroni threshold
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  # suggestive treshold
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_2_log + 1, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.873, 
           y = -log10(1e-5) -1.5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  theme_minimal() +
  labs(title = "%MRWS - chr2:InvAA_chr3:1SnpBB (n=64)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))



# Combined plot ####
manhattenANDfst <- plot_grid(
  manhatten_pwald_1,
  manhatten_pwald_2,
  Fst,
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(1, 1, 0.4)
)
manhattenANDfst

output_file <- "NP25_gwas/NP25_gwas_relMRWS_mm80_homoSW/gwas_relMRWS_homoSW.svg"

if (!file.exists(output_file) || readline("File exists. Overwrite? [y/n]: ") == "y") {
  ggsave(manhattenANDfst, file = output_file, width = 9, height = 6)
}

# Manhatten PIP ####
# BSLMM - Bayesian Sparse Linear Mixed Model
# based on PIP - posterior inclusion probability

gemma_results_NP25_gamma <- fread("NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gemma_results/gemma_lmm_results.gamma.txt")
gemma_results_NP25_bv <- fread("NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gemma_results/gemma_lmm_results.bv.txt")
gemma_results_NP25_hyp <- fread("NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gemma_results/gemma_lmm_results.hyp.txt")
gemma_results_NP25_param <- fread("NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gemma_results/gemma_lmm_results.param.txt")


# Load param.txt
param <- fread("NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gemma_results/gemma_lmm_results.param.txt")

colnames(param) <- c("chr", "SNP", "ps", "n_miss", "alpha", "beta", "PIP")
param$chr <- as.factor(param$chr)

# Optional: add cumulative position for Manhattan plot
param <- param %>%
  group_by(CHR) %>%
  mutate(chr_pos = BP + min(row_number() == 1) * 1e6) %>%
  ungroup()

param <- param %>%
  left_join(chromosome_mapping, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# Plot
ggplot(param, aes(x = pos_cumulative, y = PIP, color = as.factor(chr))) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = rep(c("grey30", "grey60"), length.out = 10)) +
  labs(x = "Genomic Position", y = "Posterior Inclusion Probability (PIP)", title = "BSLMM Manhattan Plot") +
  theme_minimal() +
  theme(legend.position = "none")

manhattan_PIP <- ggplot(param, aes(x = pos_cumulative, y = PIP)) +
  geom_point(data = subset(param, PIP < 0.1),
             aes(x = pos_cumulative, y = PIP, color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(param, PIP > 0.1),
             aes(x = pos_cumulative, y = PIP),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
#  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = 0.1, color = "blue", linetype = "dashed", size = 0.5) +  
#  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.88, 
#           y = bonferroni_threshold_1_log + 1.5, label = "Bonferroni threshold", color = "red", size = 4) +
#  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.873, 
#           y = -log10(1e-5) - 1.5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  theme_minimal() +
  labs(title = "BSLMM: Emerge before 100s or after 200s (n=192)",
       x = NULL,
       y = 'PIP') +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan

# Combined plot ####
manhattenANDfst <- grid.arrange(
  manhattan_PIP, Fst,
  nrow = 2,
  heights = c(1, 0.4))
manhattenANDfst
output_file <- "NP25_gwas/NP25_gwas_time_binary_bef100-aft200/gwas_bslmm_time_binary.svg"

if (!file.exists(output_file) || readline("File exists. Overwrite? [y/n]: ") == "y") {
  ggsave(manhattenANDfst, file = output_file, width = 9, height = 4)
}

