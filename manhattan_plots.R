library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)   # Arranging multiple grid-based plots


dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef"
setwd(dir)

# Create chromosome mapping
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

# GWAS 1 ####
gemma_results_NP25_1 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80/gemma_results/gemma_lmm_results.assoc.txt")

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
gwas1 <- ggplot(gemma_results_NP25_1, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald > bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald < bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_log, color = "blue", linetype = "dashed", size = 0.5) +  # Bonferroni threshold
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 0.5) +  # p = 0.05 treshold
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.89, 
           y = bonferroni_threshold_log + 1.5, label = "Bonferroni threshold", color = "blue", size = 4) +
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.90, 
           y = -log10(0.05) + 1.5, label = "p < 0.05 threshold", color = "red", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  theme_minimal() +
  labs(title = "Nieuwpoort 2025 (n=192) - %MRWS",
       x = NULL,
       y = expression(-log[10](p-value))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 16))
gwas1

# GWAS 2 ####
gemma_results_NP25_2 <- fread("NP25_gwas/NP25_gwas_relMRWS_mm80_homoSW/gemma_results/gemma_lmm_results.assoc.txt")

gemma_results_NP25_2$chr <- as.factor(gemma_results_NP25_2$chr)  # Make sure chromosome is a factor
gemma_results_NP25_2$p_wald <- as.numeric(gemma_results_NP25_2$p_wald)
gemma_results_NP25_2$ps <- as.numeric(gemma_results_NP25_2$ps)

gemma_results_NP25_2 <- gemma_results_NP25_2 %>%
  left_join(chromosome_mapping, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# bonferroni treshold corrects for number of tests (conservative)
num_tests_1 <- nrow(gemma_results_NP25_2)  # Total number of SNPs tested
bonferroni_threshold_1 <- 0.05 / num_tests_1  # Genome-wide significance threshold
bonferroni_threshold_1_log <- -log10(bonferroni_threshold_1)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_2 <- length(unique(gemma_results_NP25_2$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_2)

# Plot2
gwas2 <- ggplot(gemma_results_NP25_2, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald > bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald < bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_log, color = "blue", linetype = "dashed", size = 0.5) +  # Bonferroni threshold
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 0.5) +  # p = 0.05 treshold
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.89, 
           y = bonferroni_threshold_log + 1, label = "Bonferroni threshold", color = "blue", size = 4) +
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.90, 
           y = -log10(0.05) + 1, label = "p < 0.05 threshold", color = "red", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome) +
  theme_minimal() +
  labs(title = "Nieuwpoort 2025 homozygotes SW (n=120) - %MRWS",
       x = NULL,
       y = expression(-log[10](p-value))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 16))
gwas2

# Fst ####
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
  theme_minimal() +
  labs(x = "Chromosome",
       y = 'Fst') +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))

gwasANDfst <- grid.arrange(
  gwas1, gwas2, Fst,
  nrow = 3,
  heights = c(1, 1, 0.8))

gwasANDfst

output_file <- "NP25_gwas/NP25_gwas_relMRWS_mm80_homoSW/gwas_relMRWS_homozygotes.svg"

if (!file.exists(output_file) || readline("File exists. Overwrite? [y/n]: ") == "y") {
  ggsave(gwasANDfst, file = output_file, width = 9, height = 8)
}
