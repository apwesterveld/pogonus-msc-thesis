# collapse all: alt+cmd+O (mac)
# Initialise r ####
# loading packages, working directory and set sum contrast 

library(openxlsx)    # Read/write Excel files
library(dplyr)
library(afex)        # ANOVA and mixed models for factorial experiments
library(lme4)        # Linear and generalized linear mixed models
library(MuMIn)       # Model selection and model averaging
library(emmeans)     # Estimated marginal means (least-squares means)
library(effects)     # Effect plots for regression models
library(ggplot2)     # Data visualization with grammar of graphics (in tidyverse)
library(ggthemes)    # Additional themes and scales for ggplot2
library(grid)        # Low-level graphics system
library(gridExtra)   # Arranging multiple grid-based plots
library(scales)      # for percent_format()
library(ggpubr)      # in combined plots, ex rremove("ylab")
library(readxl)      # for old .XLS format


# Set working directory
dir <- "..."
setwd(dir)

# to compare levels (when using factors) to the grand mean and coefficients sum to zero
set_sum_contrasts()

# Load data ####
data <- read.xlsx("NP25_phenotypes/Pogonus_Nieuwpoort_2025_local.xlsx")

# transformations
data$IID <- as.factor(data$IID)
data$sex <- as.factor(data$sex)
data$action1 <- factor(data$action1, levels = c("submerge", "emerge"))
data$action2 <- factor(data$action2, levels = c("submerge", "emerge"))
data$action3 <- factor(data$action3, levels = c("submerge", "emerge"))
data$action4 <- factor(data$action4, levels = c("submerge", "emerge"))
data$action5 <- factor(data$action5, levels = c("submerge", "emerge"))
data$position300s1 <- as.factor(data$position300s1)
data$position300s2 <- as.factor(data$position300s2)
data$position300s3 <- as.factor(data$position300s3)
data$position300s4 <- as.factor(data$position300s4)
data$position300s5 <- as.factor(data$position300s5)

# Calculate %MRWS ####
# calculate %MRWS (Maximum Realizable Wing Size) 
# group 2: + and - functional flight muscles
# A_females = 1.2146 
# A_males = 1.1924 
# B_females = 0.8497 
# B_males = 0.8516 
# group 5: standard group (+ functional flight muscles and regularly flying)
A_females = 1.4654 
A_males = 1.4879 
B_females = 0.8454 
B_males = 0.8385 
data$Esize <- data$EL * data$EW
data$Wsize <- data$WL * data$WW

data$maxWsize <- NA

for (i in 1:nrow(data)) {
  ifelse(data$sex[i]=='F',data$maxWsize[i]<-exp(A_females) * data$Esize[i]^B_females,
         ifelse(data$sex[i]=='M',data$maxWsize[i]<-exp(A_males) * data$Esize[i]^B_males,NA))
}
data$relMRWS <- data$Wsize/data$maxWsize

# sex ratio
#m <- sum(data$sex == "M" , na.rm = TRUE) #all males
m <- count(subset(data, data$sex == "M" & relMRWS>0)) # males with %MRWS measurement
#f <- sum(data$sex == "F", na.rm = TRUE) #all females
f <- count(subset(data, data$sex == "F" & relMRWS>0)) # females with %MRWS measurement
m
f
m/f

# Summary ####

summary(data[c('EL','EW','WL','WW','relMRWS')])

# SD

sd(na.omit(data$relMRWS))

# standard error of the mean
sd(na.omit(data$WL)) / sqrt(length(na.omit(data$WL)))
sd(na.omit(data$EL)) / sqrt(length(na.omit(data$EL)))
sd(na.omit(data$relMRWS)) / sqrt(length(na.omit(data$relMRWS)))

mean(na.omit(data$WEratio))
sd(na.omit(data$WEratio)) / sqrt(length(na.omit(data$WEratio)))
mean(na.omit(data$WEsurfaceratio))
sd(na.omit(data$WEsurfaceratio)) / sqrt(length(na.omit(data$WEsurfaceratio)))

# Correlations ####

# check pearson correlation between explanatory variables
cor(na.omit(data[c('EL','WL','relMRWS')]), method = "pearson")
#cor(na.omit(data[c('EL','Esize','WL','Wsize','relMRWS')]), method = "pearson")

# EL ~ WL
# there is a correlation between EL and WL 
plot(data$EL, data$WL)
cor.test(data$EL, data$WL, method="pearson")

# EL ~ %MRWS
# using %MRWS as wing size removes the correlation with EL (and therefor probably also with sex)
plot(data$EL, data$relMRWS)
cor.test(data$EL, data$relMRWS, method="pearson")

plot(data_long_trial12345$relMRWS, data_long_trial12345$emergeCount)

plot(data$relMRWS, data$emergeCount)
cor.test(data$relMRWS, data$emergeCount, method="pearson")



cor.test(data$sex, data$emergeCount, method="pearson")


# Violin plots ####
# make plots
y_min <- 2.5 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 5 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_EL <- ggplot(data, aes(x = factor(1), y = EL)) +
  geom_violin(fill = "brown", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black",alpha = 0.7) +
  labs(x = "Elytra length", y = "mm") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank()) 

plot_WL <- ggplot(data, aes(x = factor(1), y = WL)) +
  geom_violin(fill = "lightyellow", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "Wing length", y = "mm") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

y_min <- 1 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 2.5 #max(c(data$EL, data$EW, data$WL, data$WW))

plot_EW <- ggplot(data, aes(x = factor(1), y = EW)) +
  geom_violin(fill = "brown", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "Elytra width", y = "mm") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

plot_WW <- ggplot(data, aes(x = factor(1), y = WW)) +
  geom_violin(fill = "lightyellow", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "Wing width", y = "mm") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

y_min <- 3 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 10 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_Esize <- ggplot(data, aes(x = factor(1), y = Esize)) +
  geom_violin(fill = "brown", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "Elytra size", y = "mm^2") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

y_min <- 3 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 10 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_Wsize <- ggplot(data, aes(x = factor(1), y = Wsize)) +
  geom_violin(fill = "lightyellow", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "Wing size", y = "mm^2") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

y_min <- 0.5
y_max <- 1.6

plot_WEratio <- ggplot(data, aes(x = factor(1), y = WEratio)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "WEratio", y = "mm") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

plot_WEsurfaceratio <- ggplot(data, aes(x = factor(1), y = WEsurfaceratio)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "WEsurfaceratio", y = "mm") +
  coord_cartesian(ylim = c(y_min,y_max)) +  # Use coord_cartesian to fix y-axis range
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

y_min <- 0.20
y_max <- 0.80

plot_relMWRS <- ggplot(data, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "", y = "%MRWS") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())

# show plots

grid.arrange(plot_EL, plot_WL, plot_EW, plot_WW, nrow = 2)
ggsave(
  grid.arrange(plot_EL, plot_WL, plot_EW, plot_WW, nrow = 2), 
  file="plot_morphology.svg", width=7, height=5)

#grid.arrange(plot_EL, plot_WL, plot_EW, plot_WW, plot_WEratio, plot_WEsurfaceratio, nrow = 3)

plot_relMWRS
#ggsave(plot_relMWRS, file="plot_relMWRS.svg", width=3.5, height=4.5)
ggsave(plot_relMWRS, file="plot_relMWRS.svg", width=3.5, height=3.5)

# Compare sexes ####
fit <- lm(EL ~ sex, data=data)
shapiro.test(residuals(fit))  # W > 0.9, so no deviation from normality
ncvTest(fit) #the variances do not significantly diverge from homogeneous
summary(fit) #p<0.001
emmeans(fit, ~ sex)

fit <- lm(WL ~ sex, data=data)
shapiro.test(residuals(fit))  # W > 0.9, so no deviation from normality
ncvTest(fit) #the variances do not significantly diverge from homogeneous
summary(fit) #p<0.05
emmeans(fit, ~ sex)

fit <- lm(relMRWS ~ sex, data=data)
shapiro.test(residuals(fit))  # W > 0.9, so no deviation from normality
ncvTest(fit) #the variances do not significantly diverge from homogeneous
summary(fit) # NS
emmeans(fit, ~ sex)


y_min <- 2.5 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 5 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_EL_sex <- ggplot(data, aes(x = sex, y = EL)) +
  geom_violin(fill = "brown", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
  labs(title = "Elytra length (mm)",x = NULL, y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5,size=18))
plot_EL_sex

y_min <- 2.5 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 5 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_WL_sex <- ggplot(data, aes(x = sex, y = WL)) +
  geom_violin(fill = "lightyellow", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
  labs(x = NULL, title = "Wing length (mm)",y=NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  theme_minimal(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5,size=18))

y_min <- 0 #min(c(data$EL, data$EW, data$WL, data$WW))
y_max <- 1 #max(c(data$EL, data$EW, data$WL, data$WW))
plot_relMRWS_sex <- ggplot(data, aes(x = sex, y = relMRWS)) +
  geom_violin(fill = "#4B8CEE", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.7) +
  scale_x_discrete(labels = c("F" = "Female", "M" = "Male")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = NULL,title = "%MRWS",y=NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  theme_minimal(base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5,size=18))

plot_morph_sex <- grid.arrange(plot_EL_sex, plot_WL_sex, plot_relMRWS_sex,
             nrow=1)
ggsave(
  plot_morph_sex, 
  file="results/plot_morph_sex.svg", width=10.5, height=3.5)

# Compare temporal ####
# compare %MRWS with other localities 

# import ±2011 data
olddata <- read.xlsx("resources/POGONUS_Genome_sequencing_MASTER_FILE 1.xlsx")
olddata$Locality <- as.factor(olddata$Locality)
olddata$Sex <- as.factor(olddata$Sex)
olddata$Habitat <- as.factor(olddata$Habitat)

olddata$Esize <- olddata$EL * olddata$EW
olddata$Wsize <- olddata$WL * olddata$WW
olddata$maxWsize <- NA

for (i in 1:nrow(olddata)) {
  ifelse(olddata$Sex[i]=='F',olddata$maxWsize[i]<-exp(A_females) * olddata$Esize[i]^B_females,
         ifelse(olddata$Sex[i]=='M',olddata$maxWsize[i]<-exp(A_males) * olddata$Esize[i]^B_males,NA))
}
olddata$relMRWS <- olddata$Wsize/olddata$maxWsize

ggplot(olddata, aes(x = Locality, y = relMRWS)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  labs(title = "Violin Plot of relMRWS by Locality",
       x = "Locality",
       y = "relMRWS") +
  theme_minimal()

# manually add Nieuwpoort 1995 data
NP_95_relMRWS <- c(
  40.04, 51.48, 38.98, 50.03, 57.34, 32.35, 43.16, 29.40, 29.39, 38.30,
  36.67, 58.97, 49.20, 33.92, 40.72, 54.43, 41.04, 37.90, 46.75, 66.09,
  49.64, 36.03, 48.86, 34.84, 44.90, 39.89, 52.78, 35.61, 47.67, 47.89,
  43.35, 41.67, 52.76, 41.88, 38.55, 69.11, 31.34, 49.38, 54.50, 52.90,
  58.30, 25.75, 46.36, 59.14, 56.30, 49.98, 54.62, 63.69, 41.58, 42.79,
  42.32, 37.02, 49.88, 43.30, 36.84, 51.52, 38.54, 54.83, 60.26, 37.14,
  42.82, 69.13, 52.65, 52.33
)
NP_95_relMRWS<-NP_95_relMRWS/100
data_NP_95 <- data.frame(relMRWS = NP_95_relMRWS)

# manually add Dudzele/Lissewege 2002 data
LIS_02_relMRWS <- c(
  80.12, 77.66, 83.23, 74.28, 80.27, 82.93, 81.49, 83.52, 88.44, 79.71,
  71.09, 78.43, 70.92, 73.52, 74.28, 72.56, 77.49, 76.46, 69.58, 79.80,
  65.16, 67.77, 81.79, 82.94, 78.20, 81.55, 69.80, 67.50, 78.48
)
LIS_02_relMRWS<-LIS_02_relMRWS/100
data_LIS_02 <- data.frame(relMRWS = LIS_02_relMRWS)

# test differences between NP '95, '11 and '25
sd(NP_95_relMRWS) # 0.09827353
sd(data_nieuwpoort_tidal$relMRWS) # 0.06324164
sd(na.omit(data$relMRWS)) # 0.06666338

NP_relMRWS_long <- data.frame(
  relMRWS = c(NP_95_relMRWS, data_nieuwpoort_tidal$relMRWS, na.omit(data$relMRWS)),
  group = factor(c(
    rep("NP_95", length(NP_95_relMRWS)),
    rep("NP_11", length(data_nieuwpoort_tidal$relMRWS)),
    rep("NP_25", length(na.omit(data$relMRWS)))
  ))
)
# NP_relMRWS_long <- data.frame(
#   relMRWS = c(NP_95_relMRWS, data$relMRWS),
#   group = factor(c(
#     rep("NP_95", length(NP_95_relMRWS)),
#     rep("NP_25", length(data$relMRWS))
#   ))
# )

fit <- lm(relMRWS ~ group, data=NP_relMRWS_long)
shapiro.test(residuals(fit))  # W > 0.9, so no deviation from normality
ncvTest(fit) #the variances do significantly diverge from homogeneous
spreadLevelPlot(fit)
# so use gls to allow different variances
library(nlme)
fit2 <- gls(relMRWS ~ group,weights = varIdent(form=~1|group), data=NP_relMRWS_long)
summary(fit2)
emmeans(fit, ~ group)
plot(allEffects(fit2), lty=0)

# type III testing not possible with gls, but doesn't matter with only one categorical predictor variable
# Pairwise comparisons (Tukey-adjusted)
pairwise_results <- pairs(emmeans(fit2, ~ group), adjust = "tukey")
pairwise_results
summary(pairwise_results)
# contrast      estimate     SE   df t.ratio p.value
# NP_11 - NP_25   0.0361 0.0138 29.8   2.619  0.0356
# NP_11 - NP_95  -0.0959 0.0178 64.3  -5.384  <.0001
# NP_25 - NP_95  -0.1320 0.0132 83.3 -10.007  <.0001

y_min <- 0
y_max <- 1

plot_relMWRS_NP_tidal_25 <- ggplot(data, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "#4B8CEE", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: tidal\nYser est. '25\n n=192", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMWRS_NP_tidal_25

nrow(data_NP_95)
plot_relMWRS_NP_tidal_95 <- ggplot(data_NP_95, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "#ADD8E6", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: tidal\nYser est. '95\n n=64", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMWRS_NP_tidal_95

data_nieuwpoort_tidal <- subset(olddata,Locality=="Belgium: Nieuwpoort")
nrow(data_nieuwpoort_tidal)
plot_relMRWS_nieuwpoort_tidal_11 <- ggplot(data_nieuwpoort_tidal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "#7CB5E9", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: tidal\nYser est. '11\n n=24", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMRWS_nieuwpoort_tidal_11

nrow(data_LIS_02)
plot_relMWRS_dudzele_seasonal_02 <- ggplot(data_LIS_02, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "red", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: seasonal\nDudzele '02\n n=29", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMWRS_dudzele_seasonal_02

data_dudzele_seasonal <- subset(olddata,Locality=="Belgium: Dudzele")
nrow(data_dudzele_seasonal)
plot_relMRWS_dudzele_seasonal_11 <- ggplot(data_dudzele_seasonal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "red", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: seasonal\nDudzele '11\n n=24", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMRWS_dudzele_seasonal_11

data_guérande_tidal <- subset(olddata,Locality=="France: Guérande"&Habitat=='tidal')
nrow(data_guérande_tidal)
plot_relMRWS_guérande_tidal <- ggplot(data_guérande_tidal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "blue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "FR: tidal\nGuérande '11\n n=24", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_relMRWS_guérande_tidal

data_guérande_seasonal <- subset(olddata,Locality=="France: Guérande"&Habitat=='seasonal')
plot_relMRWS_guérande_seasonal <- ggplot(data_guérande_seasonal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "red", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "FR: seasonal\nGuérande '11\n n=24", y = '%MRWS') +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
plot_relMRWS_guérande_seasonal

relMRWS_comp <- grid.arrange(
  plot_relMRWS_guérande_seasonal,
  plot_relMRWS_guérande_tidal, 
  plot_relMWRS_dudzele_seasonal_02,
  plot_relMRWS_dudzele_seasonal_11, 
  plot_relMWRS_NP_tidal_95,
  plot_relMRWS_nieuwpoort_tidal_11,
  plot_relMWRS_NP_tidal_25,
  nrow = 1,widths = c(2.8, 2, 2, 2, 2, 2, 2))

relMRWS_comp
ggsave(
  relMRWS_comp, 
  file="results/plot_relMRWS_comp.svg", width=12, height=4)

