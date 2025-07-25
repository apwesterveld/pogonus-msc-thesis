library(data.table) #fread
library(ggplot2)
library(car)        # ANOVA and mixed models for factorial experiments
library(emmeans)
library(MuMIn)       # Model selection and model averaging
library(nlme)
library(effects)

dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef"
setwd(dir)

## data ####
phenotypes <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_phenotypes/phenotype.txt")
phenotypes <- na.omit(phenotypes)
# chr2
genotype_chr2_aa1 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr2/karyotype_chr2_aa1.txt", header = FALSE)
genotype_chr2_aA2 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr2/karyotype_chr2_aA2.txt", header = FALSE)
genotype_chr2_AA3 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr2/karyotype_chr2_AA3.txt", header = FALSE)

#choose region of chr3
# CHR3:38310157-38330157
genotype_chr3_bb1 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3/karyotype_chr3_bb1.txt", header = FALSE)
genotype_chr3_bB2 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3/karyotype_chr3_bB2.txt", header = FALSE)
genotype_chr3_BB3 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3/karyotype_chr3_BB3.txt", header = FALSE)

# chr3_2:34874684 to 34909607
genotype_chr3_bb1 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3_2/karyotype_chr3_bb1.txt", header = FALSE)
genotype_chr3_bB2 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3_2/karyotype_chr3_bB2.txt", header = FALSE)
genotype_chr3_BB3 <- fread("/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/NP25_PCA/NP25_PCA_chr3_2/karyotype_chr3_BB3.txt", header = FALSE)

# chr3_2:1SNP
genotype_chr3_bb1 <- fread("NP25_PCA/NP25_PCA_chr3_2/karyotype_1SNPchr3_bb1.txt", header = FALSE)
genotype_chr3_bB2 <- fread("NP25_PCA/NP25_PCA_chr3_2/karyotype_1SNPchr3_bB2.txt", header = FALSE)
genotype_chr3_BB3 <- fread("NP25_PCA/NP25_PCA_chr3_2/karyotype_1SNPchr3_BB3.txt", header = FALSE)

## name genotypes ####
# genotype chr2 gwas peak
phenotypes$chr2 <- NA
phenotypes$chr2[phenotypes$IID %in% genotype_chr2_aa1$V1] <- "chr2:aa"
phenotypes$chr2[phenotypes$IID %in% genotype_chr2_aA2$V1] <- "chr2:aA"
phenotypes$chr2[phenotypes$IID %in% genotype_chr2_AA3$V1] <- "chr2:AA"

# genotype chr3 gwas peak
phenotypes$chr3 <- NA
phenotypes$chr3[phenotypes$IID %in% genotype_chr3_bb1$V1] <- "chr3:bb"
phenotypes$chr3[phenotypes$IID %in% genotype_chr3_bB2$V1] <- "chr3:bB"
phenotypes$chr3[phenotypes$IID %in% genotype_chr3_BB3$V1] <- "chr3:BB"

# Combine chr2 and chr3 into a single factor
phenotypes$genotype <- interaction(phenotypes$chr2, phenotypes$chr3, drop = TRUE)

# Clean up labels (e.g., "aA.AA" → "aA x AA")
levels(phenotypes$genotype) <- gsub("\\.", "_x_", levels(phenotypes$genotype))

## define order
# phenotypes$karyo_combo <- factor(
#   phenotypes$karyo_combo,
#   levels = c("aa x bb", "aa x bB", "aa x BB",
#              "aA x bb", "aA x bB", "aA x BB",
#              "AA x bb", "AA x bB", "AA x BB")
# )
phenotypes$genotype <- factor(
  phenotypes$genotype,
  levels = c("chr2:aa_x_chr3:bb", "chr2:aa_x_chr3:bB", "chr2:aa_x_chr3:BB",
             "chr2:aA_x_chr3:bb", "chr2:aA_x_chr3:bB", "chr2:aA_x_chr3:BB",
             "chr2:AA_x_chr3:bb", "chr2:AA_x_chr3:bB", "chr2:AA_x_chr3:BB"))

## epistatis other way around
# phenotypes$genotype <- factor(
#   phenotypes$genotype,
#   levels = c("chr2:aa_x_chr3:bb", "chr2:aA_x_chr3:bb", "chr2:AA_x_chr3:bb",
#              "chr2:aa_x_chr3:bB", "chr2:aA_x_chr3:bB", "chr2:AA_x_chr3:bB",
#              "chr2:aa_x_chr3:BB", "chr2:aA_x_chr3:BB", "chr2:AA_x_chr3:BB"))

# export
#write.table(phenotypes, "phenotypes_karyotypes.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)


## plotting ####
plot_epistasis <- ggplot(phenotypes, aes(x = genotype, y = relMRWS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    x = "genotype (chr2 x chr3)",
    y = "Relative MRWS",
    title = "relMRWS across genotype Combinations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
plot_epistasis
ggsave(plot_epistasis, file = "NP25_PCA_chr3_2/NP25_epistasis_chr2_1SNPchr3.svg", width = 8, height = 5)

iid_lists <- split(phenotypes$IID, phenotypes$genotype)
write.table(iid_lists[["chr2:AA_x_chr3:BB"]]
, "NP25_PCA_karyotyping/NP25_PCA_chr3_2/genotype_InvChr2AA_1SnpChr3BB.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


fit <- lm(relMRWS ~ genotype, data = phenotypes)
Anova(fit, type='III')
ncvTest(fit) # non equal variance!

fit <- gls(relMRWS ~ genotype, data = phenotypes)
Anova(fit, type='III')

# Get estimated means and 95% CIs
emm <- emmeans(fit, ~ genotype)
emm_df <- as.data.frame(emm)

pairs(emm, adjust = "tukey")

# Effect plot
ggplot(emm_df, aes(x = genotype, y = emmean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  labs(
    x = "genotype (chr2 x chr3)",
    y = "Estimated relMRWS",
    title = "Effect Plot of genotype Combinations"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# formal test of epistasis
fit <- lm(relMRWS ~ chr2 * chr3, data=phenotypes)
Anova(fit, type='III')
ncvTest(fit) # non equal variance!

fit1 <- gls(relMRWS ~ chr3 * chr2, data=phenotypes)
fit2 <- gls(relMRWS ~ chr3 + chr2, data=phenotypes)
fit3 <- gls(relMRWS ~ chr2, data=phenotypes)
fit4 <- gls(relMRWS ~ chr3, data=phenotypes)
AICc(fit1,fit2,fit3,fit4)
Anova(fit1, type='III')
plot(allEffects(fit1))

Anova(fit2, type='III')
Anova(fit3, type='III')
summary(fit3)
plot(allEffects(fit3))
Anova(fit4, type='III')
plot(allEffects(fit4))

summary(fit1)
table(phenotypes$chr2, phenotypes$chr3)
chisq.test(phenotypes$chr2, phenotypes$chr3)
cor.test(phenotypes$chr2, phenotypes$chr3, method="pearson")
fisher.test(table(phenotypes$chr2, phenotypes$chr3))

phenotypes$fitted <- fitted(fit1)

# Plot interaction
ggplot(phenotypes, aes(x = chr3, y = fitted, color = chr2, group = chr2)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(title = "Fitted Interaction: chr2 × chr3", y = "Fitted relMRWS") +
  theme_minimal()

emm <- emmeans(fit1, ~ chr2 * chr3)
emm
counts <- phenotypes %>%
  group_by(chr2, chr3) %>%
  summarise(n = n(), .groups = "drop")

# 4. Merge counts with EMMs
emm_df <- left_join(emm_df, counts, by = c("chr2", "chr3"))

# 5. Plot with sample sizes as labels
relMRWSchr3chr2 <- ggplot(emm_df, aes(x = chr3, y = emmean, color = chr2, group = chr2)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_line(position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.5),
            vjust = -1.2, size = 3.5, show.legend = FALSE) +
  labs(
    y = "%MRWS",
    x = NULL,
    title = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8,0.8),legend.title = element_blank())

ggsave(relMRWSchr3chr2, file = "epistatis_mrws_chr3chr2.svg", width = 4, height = 3.5)

relMRWSchr2chr3 <- ggplot(emm_df, aes(x = chr2, y = emmean, color = chr3, group = chr3)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_line(position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.5),
            vjust = -1.2, size = 3.5, show.legend = FALSE) +
  labs(
    y = "%MRWS",
    x = NULL,
    title = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8,0.8),legend.title = element_blank())
relMRWSchr2chr3
ggsave(relMRWSchr2chr3, file = "epistatis_mrws_chr2chr3.svg", width = 4, height = 3.5)
