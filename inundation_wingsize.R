# collapse all: alt+cmd+O (mac)
# Initialise r ####
# loading packages, working directory and set sum contrast 

library(openxlsx)    # Read/write Excel files
library(tidyverse)   # loads dplyr, ggplot2, purrr, tidyr
#library(dplyr)       # Data manipulation and transformation (in tidyverse)
library(afex)        # ANOVA and mixed models for factorial experiments
library(reshape2)    # Data reshaping (wide <-> long format)
library(irr)         # Inter-rater reliability metrics for repeatability checking
library(lattice)     # for xyplot
library(nnet)        # For glmer. Multinomial logistic regression and neural networks
library(lme4)        # Linear and generalized linear mixed models
library(MuMIn)       # Model selection and model averaging
library(car)         # Companion to applied regression; includes Anova for glmer
library(emmeans)     # Estimated marginal means (least-squares means)
library(effects)     # Effect plots for regression models
#library(ggplot2)     # Data visualization with grammar of graphics (in tidyverse)
library(ggthemes)    # Additional themes and scales for ggplot2
library(grid)        # Low-level graphics system
library(gridExtra)   # Arranging multiple grid-based plots
library(survival)    # Survival analysis
library(flexsurv)    # Flexible parametric survival models
library(scales)      # for percent_format()
library(ggfortify)   # Plot survival analysis
library(ggpubr)      # in combined plots, ex rremove("ylab")



# Set working directory
dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/results/"
setwd(dir)

# to compare levels (when using factors) to the grand mean and coefficients sum to zero
set_sum_contrasts()

# Load data ####
#data <- read.xlsx("CLEANED.xlsx")
data <- read.xlsx("Pogonus_Nieuwpoort_2025_local.xlsx")

for (i in 1:5) {
  action_col <- paste0("action", i)
  time_col <- paste0("time", i)
  new_col <- paste0("actionInduced", i)
  
  data[[new_col]] <- ifelse(
    data[[time_col]] == 0, "emerge",
    ifelse(data[[time_col]] == 300, "submerge", data[[action_col]])
  )
}

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

# calculate %MRWS (Maximum Realizable Wing Size)
A_females = 1.2146 
A_males = 1.1924 
B_females = 0.8497 
B_males = 0.8516 
data$Esize <- data$EL * data$EW
data$Wsize <- data$WL * data$WW
data$maxWsize[data$Sex == "F"] <- exp(A_females) * data$Esize^B_females
data$maxWsize[data$Sex == "M"] <- exp(A_males) * data$Esize^B_males
data$relMRWS <- data$Wsize/data$maxWsize

# calculate some other possible variables
data$WEratio <- data$WL/data$EL
data$Wsurface <- 3.14 * data$WL/2 * data$WW/2 # assuming oval wing
data$Esurface <- sqrt(data$EL^2 - data$EW^2) * data$EW # assuming rectangular elytra
data$WEsurfaceratio <- data$Wsurface/data$Esurface

data$emergeCount <- rowSums(data[
  c('actionInduced1', 'actionInduced2', 'actionInduced3', 'actionInduced4', 'actionInduced5')] == "emerge",
  na.rm = TRUE)
data$emergeCount3trials <- rowSums(data[
  c('actionInduced3', 'actionInduced4', 'actionInduced5')] == "emerge",
  na.rm = TRUE)


# put data in long format and transform some variables
data_long_time <- melt(data, id.vars = c("IID",'EL','Esize','WL','WEratio','WEsurfaceratio','relMRWS','sex','emergeCount','emergeCount3trials'), measure.vars = c('time1','time2','time3','time4','time5'),
                  variable.name = "trial", value.name = "time")
data_long_action <- melt(data, id.vars = c("IID",'EL','Esize','WL','WEratio','WEsurfaceratio','relMRWS','sex','emergeCount','emergeCount3trials'), measure.vars = c('action1','action2','action3','action4','action5'),
                           variable.name = "trial", value.name = "action")
data_long_position <- melt(data, id.vars = c("IID",'EL','Esize','WL','WEratio','WEsurfaceratio','relMRWS','sex','emergeCount','emergeCount3trials'), measure.vars = c('position300s1','position300s2','position300s3','position300s4','position300s5'),
                  variable.name = "trial", value.name = "position")

data_long_time$trial <- gsub("time", "", data_long_time$trial)  # Convert 'time3' to '3'
data_long_action$trial <- gsub("action", "", data_long_time$trial)  # Convert 'action3' to '3'
data_long_position$trial <- gsub("position300s", "", data_long_position$trial)  # Convert 'position300s3' to '3'

data_long <- merge(data_long_time, data_long_position, by = c("IID",'EL','Esize','WL', 'WEratio','WEsurfaceratio','relMRWS','sex','trial','emergeCount','emergeCount3trials'))
data_long <- merge(data_long, data_long_action, by = c("IID",'EL','Esize','WL', 'WEratio','WEsurfaceratio','relMRWS','sex','trial','emergeCount','emergeCount3trials'))

data_long$action[data_long$time == 0] <- "emerge"
data_long$action[data_long$time == 300] <- "submerge"

data_long$trial <- as.factor(data_long$trial)
data_long$position <- factor(data_long$position,levels = c("submerge", "emerge"))
data_long$action <- factor(data_long$action,levels = c("submerge", "emerge"))

data_long$position_numeric <- ifelse(data_long$position == 'emerge', 1, 0)
data_long$action_numeric <- ifelse(data_long$action == 'emerge', 1, 0)
data_long$positionStart <- ifelse(data_long$time==0,1,0)
data_long$censored <- ifelse(data_long$time == 300, 0, 1)

# make two subsets: (1) use all inundation trials or (2) only trial 3, 4 and 5 (because of better repeatability)
data_long_trial12345 <- na.omit(subset(data_long, select = c(-position, -position_numeric)))
data_long_trial345 <- na.omit(data_long %>% filter(trial %in% c("3", "4", "5")))

#compare with old data
olddata <- read.xlsx("POGONUS_Genome_sequencing_MASTER_FILE 1.xlsx")
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

# MORPHOLOGY ####

## Summary ####

summary(data[c('EL','EW','WL','WW','relMRWS')])

# standard error of the mean
sd(na.omit(data$WL)) / sqrt(length(na.omit(data$WL)))
sd(na.omit(data$EL)) / sqrt(length(na.omit(data$EL)))
sd(na.omit(data$relMRWS)) / sqrt(length(na.omit(data$relMRWS)))

mean(na.omit(data$WEratio))
sd(na.omit(data$WEratio)) / sqrt(length(na.omit(data$WEratio)))
mean(na.omit(data$WEsurfaceratio))
sd(na.omit(data$WEsurfaceratio)) / sqrt(length(na.omit(data$WEsurfaceratio)))

# sex ratio
#m <- sum(data$sex == "M" , na.rm = TRUE) #all males
m <- count(subset(data, data$sex == "M" & relMRWS>0)) # males with %MRWS measurement
#f <- sum(data$sex == "F", na.rm = TRUE) #all females
f <- count(subset(data, data$sex == "F" & relMRWS>0)) # females with %MRWS measurement
m
f
m/f

## Correlations ####

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

## Violin plots ####
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

# compare %MRWS with other localities
ggplot(olddata, aes(x = Locality, y = relMRWS)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  labs(title = "Violin Plot of relMRWS by Locality",
       x = "Locality",
       y = "relMRWS") +
  theme_minimal()

y_min <- 0
y_max <- 1.2

plot_relMWRS_NP_tidal <- ggplot(data, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: tidal (n=192)", y = "%MRWS") +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
plot_relMWRS_NP_tidal

data_dudzele_seasonal <- subset(olddata,Locality=="Belgium: Dudzele")
nrow(data_dudzele_seasonal)
plot_relMRWS_dudzele_seasonal <- ggplot(data_dudzele_seasonal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "BE: seasonal (n=24)", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
plot_relMRWS_dudzele_seasonal

data_guérande_tidal <- subset(olddata,Locality=="France: Guérande"&Habitat=='tidal')
nrow(data_guérande_tidal)
plot_relMRWS_guérande_tidal <- ggplot(data_guérande_tidal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "FR: tidal (n=24)", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
plot_relMRWS_guérande_tidal

data_guérande_seasonal <- subset(olddata,Locality=="France: Guérande"&Habitat=='seasonal')
plot_relMRWS_guérande_seasonal <- ggplot(data_guérande_seasonal, aes(x = factor(1), y = relMRWS)) +
  geom_violin(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(width = 0.1, color = "black") +
  labs(x = "FR: seasonal (n=24)", y = NULL) +
  coord_cartesian(ylim = c(y_min, y_max)) +  # Use coord_cartesian to fix y-axis range
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank())
plot_relMRWS_guérande_seasonal

relMRWS_comp <- grid.arrange(plot_relMWRS_NP_tidal, plot_relMRWS_dudzele_seasonal, plot_relMRWS_guérande_tidal, plot_relMRWS_guérande_seasonal, nrow = 1)
ggsave(
  relMRWS_comp, 
  file="plot_relMRWS_comp.svg", width=12, height=4)


# INUNDATION ####

# many more females have 300s!
data_input <- subset(data_long_trial12345, time ==300)
sum(data_input$sex == "M", na.rm = TRUE)
sum(data_input$sex == "F", na.rm = TRUE)
data_input <- subset(data_long_trial12345, time ==0)
sum(data_input$sex == "M", na.rm = TRUE)
sum(data_input$sex == "F", na.rm = TRUE)


## Repeatability ####

# Use "twoway" if all raters rate all subjects.
# agreement: beetles have to perform identical in every experiment
# consistency: beetles should keep the same rank, but means can vary (because eq. weather?)
# icc doesn't work on binary data
icc(na.omit(data[c('time1','time2',"time3",'time4','time5')]),model="twoway",type='agreement')
icc(na.omit(data[c('time1','time2',"time3",'time4','time5')]),model="twoway",type='consistency')

icc(na.omit(data[c("time3",'time4','time5')]),model="twoway",type='agreement')
icc(na.omit(data[c("time3",'time4','time5')]),model="twoway",type='consistency')

# use Fleiss' kappa for binary data with more than 2 raters/trials
kappam.fleiss(na.omit(data[c('action1','action2','action3','action4','action5')]))
kappam.fleiss(na.omit(data[c('action3','action4','action5')]))

kappam.fleiss(na.omit(data[c('actionInduced1','actionInduced2','actionInduced3','actionInduced4','actionInduced5')]))
kappam.fleiss(na.omit(data[c('actionInduced3','actionInduced4','actionInduced5')]))

kappam.fleiss(na.omit(data[c('position300s3','position300s4','position300s5')]))

# plot raw data
data_input <- data_long_trial12345
plot_time <- ggplot(data_input, aes(x = trial, y = time, fill = trial)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_boxplot(aes(x = trial, y = time), 
               width = 0.15, 
               color = "black", 
               alpha = 0.7) +
  scale_x_discrete(labels = c("Trial 1", "Trial 2", "Trial 3", "Trial 4", "Trial 5")) +
  labs(x = NULL, y = "Time until emerging (s)") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")
plot_time
ggsave(plot_time, file="plot_time.svg", width=6, height=4)

data_input <- data_long_trial12345
plot_action <- ggplot(data_input, aes(x = trial, fill = action)) +
  geom_bar() +
  scale_x_discrete(labels = c("Trial 1", "Trial 2", "Trial 3", "Trial 4", "Trial 5")) +
  labs(x = NULL, y = "Beetles") +
  theme_minimal(base_size = 18) 
plot_action
ggsave(plot_action, file="plot_action.svg", width=6, height=3)


data_input <- data_long_trial345
plot_position <- ggplot(data_input, aes(x = trial, fill = position)) +
  geom_bar() +
  scale_x_discrete(labels = c("Trial 3", "Trial 4", "Trial 5")) +
  labs(x = NULL, y = "Beetles") +
  theme_minimal(base_size = 18) 
plot_position
ggsave(plot_position, file="plot_position.svg", width=5, height=3)

## Ratio 5 trials ####
# barplot: how often beetles emerged
# count for each beetle how often they emerged
data$emergeCount <- rowSums(data[
  c('actionInduced1', 'actionInduced2', 'actionInduced3', 'actionInduced4', 'actionInduced5')] == "emerge",
  na.rm = TRUE)

# make frequency distribution
emergeDistribution <- data %>%
  count(emergeCount = factor(emergeCount, levels = 0:5)) %>%
  mutate(proportion = n / sum(n)) 

# calculate SE for binomial distribution, while estimating proportions
emergeDistribution <- emergeDistribution %>%
  mutate(
    se = sqrt((proportion * (1 - proportion)) / sum(n)),
    ymin = proportion - se,
    ymax = proportion + se
  )

emergeDistribution$emergeLabel <- paste0(emergeDistribution$emergeCount, "/5")

ggplot(emergeDistribution, aes(x = emergeLabel, y = proportion)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Proportion escape",
    y = "Beetles (%)"
  ) +
  theme_minimal(base_size = 18)

# Two barplots that compare actual data with randomized data

data_shuffled <- data.frame(matrix(nrow = nrow(data), ncol = 0))

# make 100 randomized datasets, each with 5 repeats stored temporarely in data_shuffled
# randomize per trial. So for randomized s1, pull only values from trial1, etc.
for (seed in 0:99) {
  set.seed(seed)
  
  s1 <- data$actionInduced1
  s1[!is.na(s1)] <- sample(na.omit(s1), replace = FALSE)
  
  s2 <- data$actionInduced2
  s2[!is.na(s2)] <- sample(na.omit(s2), replace = FALSE)
  
  s3 <- data$actionInduced3
  s3[!is.na(s3)] <- sample(na.omit(s3), replace = FALSE)
  
  s4 <- data$actionInduced4
  s4[!is.na(s4)] <- sample(na.omit(s4), replace = FALSE)
  
  s5 <- data$actionInduced5
  s5[!is.na(s5)] <- sample(na.omit(s5), replace = FALSE)
  
  shuffled_df <- data.frame(s1, s2, s3, s4, s5)
  
  emerge_count <- rowSums(shuffled_df == "emerge", na.rm = TRUE)
  
  data_shuffled[[paste0("randomized", seed)]] <- emerge_count
}

# make frequency distribution
# Initialize a matrix to hold the counts for each level (0 to 5) across 100 randomizations
emergeCountRandom <- matrix(0, nrow = 6, ncol = ncol(data_shuffled))
rownames(emergeCountRandom) <- 0:5  # levels of emerge count
colnames(emergeCountRandom) <- colnames(data_shuffled)

# Loop over each randomized column
for (i in seq_along(data_shuffled)) {
  counts <- table(factor(data_shuffled[[i]], levels = 0:5))
  emergeCountRandom[, i] <- counts
}

# Convert to data frame and summarize
emergeDistributionRandom <- data.frame(
  emergeCount = 0:5,
  mean = rowMeans(emergeCountRandom),
  se = apply(emergeCountRandom, 1, function(x) sd(x) / sqrt(length(x))),
  ymin = NA,
  ymax = NA
)

# Calculate confidence bands (mean ± SE)
emergeDistributionRandom$ymin <- emergeDistributionRandom$mean - emergeDistributionRandom$se
emergeDistributionRandom$ymax <- emergeDistributionRandom$mean + emergeDistributionRandom$se

# Calculate 95% confidence intervals (mean ± z*SE)
z=1.96
emergeDistributionRandom$yminCI <- emergeDistributionRandom$mean - z * emergeDistributionRandom$se
emergeDistributionRandom$ymaxCI <- emergeDistributionRandom$mean + z * emergeDistributionRandom$se

# Convert mean and bounds to proportions
emergeDistributionRandom$mean_prop <- emergeDistributionRandom$mean / nrow(data_shuffled)
emergeDistributionRandom$ymin_prop <- emergeDistributionRandom$ymin / nrow(data_shuffled)
emergeDistributionRandom$ymax_prop <- emergeDistributionRandom$ymax / nrow(data_shuffled)
emergeDistributionRandom$yminCI_prop <- emergeDistributionRandom$yminCI / nrow(data_shuffled)
emergeDistributionRandom$ymaxCI_prop <- emergeDistributionRandom$ymaxCI / nrow(data_shuffled)

# Add labels for x-axis
emergeDistributionRandom$label <- paste0(emergeDistributionRandom$emergeCount, "/5")

plot_data <- bind_rows(
  emergeDistribution %>%
    transmute(
      label = paste0(emergeCount, "/5"),
      percent = proportion * 100,
      ymin = ymin * 100,
      ymax = ymax * 100,
      type = "Actual"
    ),
  emergeDistributionRandom %>%
    transmute(
      label = paste0(emergeCount, "/5"),
      percent = mean_prop * 100,
      ymin = ymin_prop * 100,
      ymax = ymax_prop * 100,
#      ymin = yminCI_prop * 100,
#      ymax = ymaxCI_prop * 100,
      type = "Randomized"
    )
)

plot_action_randomized <- 
  ggplot(plot_data, aes(x = label, y = percent, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Proportion emerge",
    y = "Percentage of Beetles",
    fill = ""
  ) +
  scale_fill_manual(values = c("Actual" = "Dark Red", "Randomized" = "Light Coral")) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.2, 0.95))
plot_action_randomized
ggsave(plot_action_randomized, file="plot_action_randomized.svg", width=6, height=4)


# test for significant difference between distributions
# not sure if this is correct
expected_counts <- rowMeans(emergeCountRandom)  
observed_counts <- emergeDistribution$n
chisq.test(x = observed_counts, p = expected_counts / sum(expected_counts))

## Ratio 3 trials ####
# barplot: how often beetles emerged
# count for each beetle how often they emerged
data$emergeCount3trials <- rowSums(data[
  c('actionInduced3', 'actionInduced4', 'actionInduced5')] == "emerge",
  na.rm = TRUE)

# make frequency distribution
emergeDistribution <- data %>%
  count(emergeCount3trials = factor(emergeCount3trials, levels = 0:3)) %>%
  mutate(proportion = n / sum(n)) 

# calculate SE for binomial distribution, while estimating proportions
emergeDistribution <- emergeDistribution %>%
  mutate(
    se = sqrt((proportion * (1 - proportion)) / sum(n)),
    ymin = proportion - se,
    ymax = proportion + se
  )

emergeDistribution$emergeLabel <- paste0(emergeDistribution$emergeCount, "/5")

ggplot(emergeDistribution, aes(x = emergeLabel, y = proportion)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Proportion escape",
    y = "Beetles (%)"
  ) +
  theme_minimal(base_size = 18)

# Two barplots that compare actual data with randomized data

data_shuffled <- data.frame(matrix(nrow = nrow(data), ncol = 0))

# make 100 randomized datasets, each with 5 repeats stored temporarely in data_shuffled
# randomize per trial. So for randomized s1, pull only values from trial1, etc.
for (seed in 0:99) {
  set.seed(seed)
  
  #s1 <- data$actionInduced1
  #s1[!is.na(s1)] <- sample(na.omit(s1), replace = FALSE)
  
  #s2 <- data$actionInduced2
  #s2[!is.na(s2)] <- sample(na.omit(s2), replace = FALSE)
  
  s3 <- data$actionInduced3
  s3[!is.na(s3)] <- sample(na.omit(s3), replace = FALSE)
  
  s4 <- data$actionInduced4
  s4[!is.na(s4)] <- sample(na.omit(s4), replace = FALSE)
  
  s5 <- data$actionInduced5
  s5[!is.na(s5)] <- sample(na.omit(s5), replace = FALSE)
  
  shuffled_df <- data.frame(s3, s4, s5)
  
  emerge_count <- rowSums(shuffled_df == "emerge", na.rm = TRUE)
  
  data_shuffled[[paste0("randomized", seed)]] <- emerge_count
}

# make frequency distribution
# Initialize a matrix to hold the counts for each level (0 to 5) across 100 randomizations
emergeCountRandom <- matrix(0, nrow = 4, ncol = ncol(data_shuffled))
rownames(emergeCountRandom) <- 0:3  # levels of emerge count
colnames(emergeCountRandom) <- colnames(data_shuffled)

# Loop over each randomized column
for (i in seq_along(data_shuffled)) {
  counts <- table(factor(data_shuffled[[i]], levels = 0:3))
  emergeCountRandom[, i] <- counts
}

# Convert to data frame and summarize
emergeDistributionRandom <- data.frame(
  emergeCount3trials = 0:3,
  mean = rowMeans(emergeCountRandom),
  se = apply(emergeCountRandom, 1, function(x) sd(x) / sqrt(length(x))),
  ymin = NA,
  ymax = NA
)

# Calculate confidence bands (mean ± SE)
emergeDistributionRandom$ymin <- emergeDistributionRandom$mean - emergeDistributionRandom$se
emergeDistributionRandom$ymax <- emergeDistributionRandom$mean + emergeDistributionRandom$se

# Calculate 95% confidence intervals (mean ± z*SE)
z=1.96
emergeDistributionRandom$yminCI <- emergeDistributionRandom$mean - z * emergeDistributionRandom$se
emergeDistributionRandom$ymaxCI <- emergeDistributionRandom$mean + z * emergeDistributionRandom$se

# Convert mean and bounds to proportions
emergeDistributionRandom$mean_prop <- emergeDistributionRandom$mean / nrow(data_shuffled)
emergeDistributionRandom$ymin_prop <- emergeDistributionRandom$ymin / nrow(data_shuffled)
emergeDistributionRandom$ymax_prop <- emergeDistributionRandom$ymax / nrow(data_shuffled)
emergeDistributionRandom$yminCI_prop <- emergeDistributionRandom$yminCI / nrow(data_shuffled)
emergeDistributionRandom$ymaxCI_prop <- emergeDistributionRandom$ymaxCI / nrow(data_shuffled)

# Add labels for x-axis
emergeDistributionRandom$label <- paste0(emergeDistributionRandom$emergeCount, "/5")

plot_data <- bind_rows(
  emergeDistribution %>%
    transmute(
      label = paste0(emergeCount3trials, "/5"),
      percent = proportion * 100,
      ymin = ymin * 100,
      ymax = ymax * 100,
      type = "Actual"
    ),
  emergeDistributionRandom %>%
    transmute(
      label = paste0(emergeCount3trials, "/5"),
      percent = mean_prop * 100,
      ymin = ymin_prop * 100,
      ymax = ymax_prop * 100,
      #      ymin = yminCI_prop * 100,
      #      ymax = ymaxCI_prop * 100,
      type = "Randomized"
    )
)

plot_action_randomized <- 
  ggplot(plot_data, aes(x = label, y = percent, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Proportion emerge",
    y = "Percentage of Beetles",
    fill = ""
  ) +
  scale_fill_manual(values = c("Actual" = "Dark Red", "Randomized" = "Light Coral")) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.2, 0.95))
plot_action_randomized
ggsave(plot_action_randomized, file="plot_action_randomized_3trials.svg", width=6, height=4)


# test for significant difference between distributions
# not sure if this is correct
expected_counts <- rowMeans(emergeCountRandom)  
observed_counts <- emergeDistribution$n
chisq.test(x = observed_counts, p = expected_counts / sum(expected_counts))

# MODELLING ####

# Can morphology predict behaviour?

# Remarks:
# Use EL as proxy for body size
# Use %MRWS as 'wing size', proxy for adaptation to tidal habitat
# take effect of sex into account
# correct for random effect IID
# don't think it makes sense to add nonlinear models for WL or EL?

# to do: try random effect of date and time of day

# ACTION ####
# can action (submerging/emerging after first breath) be predicted by morphology?

## Data filtering ####
# use all trials
data_input <- data_long_trial12345
# OPTIONAL: use all trials, but exclude beetles that emerge ~50% of the time
data_input <- subset(data_long_trial12345, emergeCount != 2 & emergeCount != 3)
# OPTIONAL: use all trials, but only include beetles with 'extreme' times 
data_input <- subset(data_long_trial12345, time < 50 | time > 250)
# OPTIONAL: use only trial 3, 4 and 5
data_input <- data_long_trial345
# OPTIONAL: use only trial 3, 4 and 5, and only include beetles with 'extreme' times 
data_input <- subset(data_long_trial345, time < 50 | time > 250)
# OPTIONAL: only males
data_input <- subset(data_long_trial12345, sex == 'M')

# check trends by eye
xyplot(action ~ WL, data=data_input ,type=c("p","r"))
xyplot(action ~ relMRWS, data=data_input ,type=c("p","r"))
xyplot(action ~ EL, data=data_input ,type=c("p","r"))
xyplot(action ~ sex, data=data_input ,type=c("p","r"))

## WL  ####
fit1a=glmer(action ~ EL*WL*sex+(1|IID),family=binomial,data=data_input)
fit1b=glmer(action ~ EL*WL+sex+(1|IID),family=binomial,data=data_input)
fit1c=glmer(action ~ EL+WL+sex+(1|IID),family=binomial,data=data_input)
fit1d=glmer(action ~ EL*WL+(1|IID),family=binomial,data=data_input)
fit1e=glmer(action ~ EL+WL+(1|IID),family=binomial,data=data_input)
fit1f=glmer(action ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit1g=glmer(action ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit1h=glmer(action ~ WL*sex+(1|IID),family=binomial,data=data_input)
fit1i=glmer(action ~ WL+sex+(1|IID),family=binomial,data=data_input)
fit1j=glmer(action ~ EL+(1|IID),family=binomial,data=data_input)
fit1k=glmer(action ~ WL+(1|IID),family=binomial,data=data_input)
fit1l=glmer(action ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g,fit1h,fit1i,fit1j,fit1k,fit1l)
Anova(fit1g,type="III") # significant effect of sex
Anova(fit1i,type="III") # significant effect of sex
Anova(fit1l,type="III") # significant effect of sex
plot(allEffects(fit1g), type='response') 
plot(allEffects(fit1i), type='response') 
plot(allEffects(fit1l), type='response') 

## WEratio (skip) ####
fit2a=glmer(action ~ EL*WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2b=glmer(action ~ EL*WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2c=glmer(action ~ EL+WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2d=glmer(action ~ EL*WEratio+(1|IID),family=binomial,data=data_input)
fit2e=glmer(action ~ EL+WEratio+(1|IID),family=binomial,data=data_input)
fit2f=glmer(action ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit2g=glmer(action ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit2h=glmer(action ~ WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2i=glmer(action ~ WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2j=glmer(action ~ EL+(1|IID),family=binomial,data=data_input)
fit2k=glmer(action ~ WEratio+(1|IID),family=binomial,data=data_input)
fit2l=glmer(action ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit2a,fit2b,fit2c,fit2d,fit2e,fit2f,fit2g,fit2h,fit2i,fit2j,fit2k,fit2l)
Anova(fit2l,type="III") # significant effect of sex

## WEsurfaceratio (skip) ####
fit3a=glmer(action ~ EL*WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3b=glmer(action ~ EL*WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3c=glmer(action ~ EL+WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3d=glmer(action ~ EL*WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3e=glmer(action ~ EL+WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3f=glmer(action ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit3g=glmer(action ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit3h=glmer(action ~ WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3i=glmer(action ~ WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3j=glmer(action ~ EL+(1|IID),family=binomial,data=data_input)
fit3k=glmer(action ~ WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3l=glmer(action ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit3a,fit3b,fit3c,fit3d,fit3e,fit3f,fit3g,fit3h,fit3i,fit3j,fit3k,fit3l)
Anova(fit3g,type="III") # significant effect of sex
Anova(fit3i,type="III") # significant effect of sex
Anova(fit3l,type="III") # significant effect of sex

## %MRWS ####
fit4a=glmer(action ~ EL*relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4b=glmer(action ~ EL*relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4c=glmer(action ~ EL+relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4d=glmer(action ~ EL*relMRWS+(1|IID),family=binomial,data=data_input)
fit4e=glmer(action ~ EL+relMRWS+(1|IID),family=binomial,data=data_input)
fit4f=glmer(action ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit4g=glmer(action ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit4h=glmer(action ~ relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4i=glmer(action ~ relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4j=glmer(action ~ EL+(1|IID),family=binomial,data=data_input)
fit4k=glmer(action ~ relMRWS+(1|IID),family=binomial,data=data_input)
fit4l=glmer(action ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit4a,fit4b,fit4c,fit4d,fit4e,fit4f,fit4g,fit4h,fit4i,fit4j,fit4k,fit4l)
Anova(fit4i,type="III") # significant effect of sex
Anova(fit4j,type="III")
Anova(fit4k,type="III") 
Anova(fit4l,type="III") # significant effect of sex 
plot(allEffects(fit4i), type='response') 
plot(allEffects(fit4l), type='response') 


# POSITION ####
# Position: were beetles emerged/submerged after 300s
# data only available for trial 3, 4 and 5
data_input <- data_long_trial345
# OPTIONAL: only include beetles with 'extreme' times 
data_input <- subset(data_long_trial345, time < 50 | time > 250)



# check trends by eye
xyplot(position ~ WL, data=data_input ,type=c("p","r"))
xyplot(position ~ relMRWS, data=data_input ,type=c("p","r"))
xyplot(position ~ EL, data=data_input ,type=c("p","r"))
xyplot(position ~ sex, data=data_input ,type=c("p","r"))

## WL ####
fit1a=glmer(position ~ EL*WL*sex+(1|IID),family=binomial,data=data_input)
fit1b=glmer(position ~ EL*WL+sex+(1|IID),family=binomial,data=data_input)
fit1c=glmer(position ~ EL+WL+sex+(1|IID),family=binomial,data=data_input)
fit1d=glmer(position ~ EL*WL+(1|IID),family=binomial,data=data_input)
fit1e=glmer(position ~ EL+WL+(1|IID),family=binomial,data=data_input)
fit1f=glmer(position ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit1g=glmer(position ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit1h=glmer(position ~ WL*sex+(1|IID),family=binomial,data=data_input)
fit1i=glmer(position ~ WL+sex+(1|IID),family=binomial,data=data_input)
fit1j=glmer(position ~ EL+(1|IID),family=binomial,data=data_input)
fit1k=glmer(position ~ WL+(1|IID),family=binomial,data=data_input)
fit1l=glmer(position ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g,fit1h,fit1i,fit1j,fit1k,fit1l)
Anova(fit1i,type="III") 
Anova(fit1l,type="III") 
plot(allEffects(fit1i), type='response') 
plot(allEffects(fit1l), type='response') 
plot(allEffects(fit1k), type='response') 

## WEratio ####
fit2a=glmer(position ~ EL*WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2b=glmer(position ~ EL*WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2c=glmer(position ~ EL+WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2d=glmer(position ~ EL*WEratio+(1|IID),family=binomial,data=data_input)
fit2e=glmer(position ~ EL+WEratio+(1|IID),family=binomial,data=data_input)
fit2f=glmer(position ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit2g=glmer(position ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit2h=glmer(position ~ WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2i=glmer(position ~ WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2j=glmer(position ~ EL+(1|IID),family=binomial,data=data_input)
fit2k=glmer(position ~ WEratio+(1|IID),family=binomial,data=data_input)
fit2l=glmer(position ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit2a,fit2b,fit2c,fit2d,fit2e,fit2f,fit2g,fit2h,fit2i,fit2j,fit2k,fit2l)
Anova(fit2l,type="III") # significant effect of sex 0.01972 *

## WEsurfaceratio ####
fit3a=glmer(position ~ EL*WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3b=glmer(position ~ EL*WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3c=glmer(position ~ EL+WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3d=glmer(position ~ EL*WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3e=glmer(position ~ EL+WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3f=glmer(position ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit3g=glmer(position ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit3h=glmer(position ~ WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3i=glmer(position ~ WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3j=glmer(position ~ EL+(1|IID),family=binomial,data=data_input)
fit3k=glmer(position ~ WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3l=glmer(position ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit3a,fit3b,fit3c,fit3d,fit3e,fit3f,fit3g,fit3h,fit3i,fit3j,fit3k,fit3l)
Anova(fit3l,type="III") # significant effect of sex 0.01972 *

## %MRWS ####
fit4a=glmer(position ~ EL*relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4b=glmer(position ~ EL*relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4c=glmer(position ~ EL+relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4d=glmer(position ~ EL*relMRWS+(1|IID),family=binomial,data=data_input)
fit4e=glmer(position ~ EL+relMRWS+(1|IID),family=binomial,data=data_input)
fit4f=glmer(position ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit4g=glmer(position ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit4h=glmer(position ~ relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4i=glmer(position ~ relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4j=glmer(position ~ EL+(1|IID),family=binomial,data=data_input)
fit4k=glmer(position ~ relMRWS+(1|IID),family=binomial,data=data_input)
fit4l=glmer(position ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit4a,fit4b,fit4c,fit4d,fit4e,fit4f,fit4g,fit4h,fit4i,fit4j,fit4k,fit4l)
Anova(fit4f,type="III") 
Anova(fit4i,type="III") 
Anova(fit4l,type="III") 
plot(allEffects(fit4i), type='response')


# TIME ####
## Data filtering ####
# can time (until first breath) be predicted by morphology?

# check distribution of continuous time Variable
data_input <- data_long_trial12345
hist(data_input$time, breaks = 20, xlab = "Time (s)")

# Exclude 0 and 300s to approach a normal distribution
data_input <- subset(data_long_trial12345, time != 0 & time != 300)
hist(data_input$time, breaks = 20)
shapiro.test(data_input$time) #W>0.9

# OPTIONAL: better with only trial 3, 4 and 5?
data_input <- subset(data_long_trial345, time != 0 & time != 300)
hist(data_input$time, breaks = 20)
shapiro.test(data_input$time) #W>0.9

# OPTIONAL: try log transformation to normalize data
hist(log(data_input$time))
shapiro.test(log(data_input$time))
# log transformed is little bit better, but not worth it?

# check trend by eye
xyplot(time ~ WL, data=data_input ,type=c("p","r"))
xyplot(time ~ relMRWS, data=data_input ,type=c("p","r"))
xyplot(time ~ EL, data=data_input ,type=c("p","r"))
xyplot(time ~ sex, data=data_input ,type=c("p","r"))

## WL ####
fit1a <- lmer(time ~ EL*WL*sex + (1|IID),data=data_input)
fit1b <- lmer(time ~ EL*WL+sex + (1|IID),data=data_input)
fit1c <- lmer(time ~ EL+WL*sex + (1|IID),data=data_input)
fit1d <- lmer(time ~ EL+WL+sex + (1|IID),data=data_input)
fit1e <- lmer(time ~ EL*WL + (1|IID),data=data_input)
fit1f <- lmer(time ~ EL+WL + (1|IID),data=data_input)
fit1g <- lmer(time ~ EL*sex + (1|IID),data=data_input)
fit1h <- lmer(time ~ EL+sex + (1|IID),data=data_input)
fit1i <- lmer(time ~ WL*sex + (1|IID),data=data_input)
fit1j <- lmer(time ~ WL+sex + (1|IID),data=data_input)
fit1k <- lmer(time ~ EL + (1|IID),data=data_input)
fit1l <- lmer(time ~ WL + (1|IID),data=data_input)
fit1m <- lmer(time ~ sex + (1|IID),data=data_input)
MuMIn::AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g,fit1h,fit1i,fit1j,fit1k,fit1l,fit1m)
Anova(fit1e, type = 'III')
Anova(fit1a, type = 'III')
plot(allEffects(fit1a), type='response') 

## %MRWS ####
fit4a <- lmer(time ~ relMRWS*EL*sex + (1|IID), data=data_input)
fit4b <- lmer(time ~ relMRWS*EL+sex + (1|IID), data=data_input)
fit4c <- lmer(time ~ EL+relMRWS*sex + (1|IID), data=data_input)
fit4d <- lmer(time ~ EL+relMRWS+sex + (1|IID), data=data_input)
fit4e <- lmer(time ~ relMRWS*EL + (1|IID), data=data_input)
fit4f <- lmer(time ~ EL+relMRWS + (1|IID), data=data_input)
fit4g <- lmer(time ~ EL*sex + (1|IID), data=data_input)
fit4h <- lmer(time ~ EL+sex + (1|IID), data=data_input)
fit4i <- lmer(time ~ relMRWS*sex + (1|IID), data=data_input)
fit4j <- lmer(time ~ relMRWS+sex + (1|IID), data=data_input)
fit4k <- lmer(time ~ EL + (1|IID), data=data_input)
fit4l <- lmer(time ~ relMRWS + (1|IID), data=data_input)
fit4m <- lmer(time ~ sex + (1|IID), data=data_input)
MuMIn::AICc(fit4a,fit4b,fit4c,fit4d,fit4e,fit4f,fit4g,fit4h,fit4i,fit4j,fit4k,fit4l,fit4m)
Anova(fit4a, type = 'III')
Anova(fit4b, type = 'III')
Anova(fit4e, type = 'III')
Anova(fit4l, type = 'III')
plot(allEffects(fit4a), type='response') 
plot(allEffects(fit4e), type='response') 

## old code to make a nice plot ####
# Get effect object
Esize_cuts <- quantile(data_input$Esize, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
data_input$Esize_tertile <- cut(data_input$Esize, 
                                breaks = Esize_cuts, 
                                include.lowest = TRUE, 
                                labels = c("Low", "Medium", "High"))
esize_representatives <- data_input %>%
  group_by(Esize_tertile) %>%
  summarize(Esize_val = median(Esize, na.rm = TRUE)) %>%
  arrange(Esize_tertile) %>%
  pull(Esize_val)

MRWS_range <- range(data_input$relMRWS, na.rm = TRUE)

eff <- Effect(c("relMRWS", "Esize"), fit5a, 
              xlevels = list(
                Esize = esize_representatives,
                relMRWS = seq(MRWS_range[1], MRWS_range[2], length.out = 100)
              ))

eff_df <- as.data.frame(eff)
eff_df$Esize_tertile <- eff_df$Esize

# Replace with your actual data name if different
#data_input$Esize_facet <- sapply(data_input$Esize, function(x) {
#  closest <- which.min(abs(x - c(3.75, 4.0, 4.25)))
#  return(c(3.75, 4.0, 4.25)[closest])
#})
data_input$Esize_facet <- factor(data_input$Esize_tertile,
                                 levels = c("Low", "Medium", "High"),
                                 labels = c("Low Esize", "Medium Esize", "High Esize"))

ggplot() +
  # Effect lines and ribbons
  geom_ribbon(data = eff_df,
              aes(x = relMRWS, ymin = lower, ymax = upper),
              fill = "lightblue", alpha = 0.3) +
  geom_line(data = eff_df,
            aes(x = relMRWS, y = fit),
            color = "blue", size = .1) +
  
  # Raw data
  geom_point(data = data_input,
             aes(x = relMRWS, y = time),
             color = "black", alpha = 0.4) +
  
  # Facet by Esize tertile
  facet_wrap(~ Esize_facet, labeller = label_bquote('Elytra size (mm2)' == .(Esize_facet))) +
  
  
  labs(x = "relMRWS", y = "Time",
       title = "Effect of relMRWS on Time at Esize Tertiles") +
  theme_minimal(base_size = 14)


MuMIn::AICc(fit1a,fit2a,fit3a,fit4a)

EL_levels <- quantile(data_input$EL, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
em_df <- data.frame(emmeans(
  fit4a,
  ~ relMRWS | EL,
  at = list(
    relMRWS = seq(min(data_input$relMRWS, na.rm = TRUE),
                  max(data_input$relMRWS, na.rm = TRUE),
                  length.out = 100),
    EL = EL_levels,
    sex = "F"  # fixed level to isolate relMRWS × EL
  ),
  type = "response"
))
install.packages("ggnewscale")
library(ggnewscale)

ggplot() +
  # Observed data: colored by EL (continuous)
  geom_point(data = data_input,
             aes(x = relMRWS, y = time, color = EL),
             alpha = 0.6, size = 1.7) +
  
  # Continuous color scale for observed points
  scale_color_gradient(low = "skyblue", high = "darkblue", name = "Observed EL") +
  
  # Reset color scale before predicted lines
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Confidence ribbons from model
  geom_ribbon(data = em_df,
              aes(x = relMRWS, ymin = lower.CL, ymax = upper.CL, fill = as.factor(EL)),
              alpha = 0.2, color = NA) +
  
  # Predicted lines: new color scale (discrete)
  geom_line(data = em_df,
            aes(x = relMRWS, y = emmean, color = as.factor(EL)),
            linewidth = 1.2) +
  
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"), name = "Predicted EL") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"), name = "Predicted EL") +
  
  theme_few(base_size = 16) +
  labs(
    x = "relMRWS",
    y = "Time",
    title = "Predicted Time vs relMRWS by EL (Observed & Fitted)"
  ) +
  theme(legend.position = "top")

em_df <- data.frame(emmeans(
  fit4a,
  ~ relMRWS | sex,
  at = list(
    relMRWS = seq(min(data_input$relMRWS, na.rm = TRUE),
                  max(data_input$relMRWS, na.rm = TRUE),
                  length.out = 100),
    EL = mean(data_input$EL, na.rm = TRUE)  # Hold EL constant at its mean
  ),
  type = "response"
))
head(em_df)
colnames(em_df) # use these column names in your ggplot2 code

ggplot(em_df, aes(x = relMRWS, y = emmean, color = sex, fill = sex)) +
  geom_point(data = data_input, 
             aes(x = relMRWS, y = time, color = sex), 
             alpha = 0.4, size = 1.5) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  theme_few(base_size = 16) +
  labs(x = "relMRWS", y = "Predicted Time", title = "Effect of relMRWS on Time by Sex") +
  theme(legend.position = "top")


# TIME binomial trans. ####
## time50 ####
# Find time50: half of the beetles emerged before that time
data_input <- data_long_trial12345 
#data_input <- data_long_trial345 

# OPTIONAL: consider experiment only starting at t=30s to compare with non parametric survival analysis
data_input <- data_long_trial12345 
data_input <- subset(data_input,time>20)

#time50 <- median(data_input$time)
T_cutoff <- 95
data_input$emerged_T_cutoff <- ifelse(data_input$time<T_cutoff,1,0)
count(data_input, emerged_T_cutoff) #check if data is split 50/50

# check trends by eye
xyplot(emerged_T_cutoff ~ WL, data=data_input ,type=c("p","r"))
xyplot(emerged_T_cutoff ~ relMRWS, data=data_input ,type=c("p","r"))
xyplot(emerged_T_cutoff ~ EL, data=data_input ,type=c("p","r"))
xyplot(emerged_T_cutoff ~ sex, data=data_input ,type=c("p","r"))

### WL ####
fit1a=glmer(emerge_time50 ~ EL*WL*sex+(1|IID),family=binomial,data=data_input)
fit1b=glmer(emerge_time50 ~ EL*WL+sex+(1|IID),family=binomial,data=data_input)
fit1c=glmer(emerge_time50 ~ EL+WL+sex+(1|IID),family=binomial,data=data_input)
fit1d=glmer(emerge_time50 ~ EL*WL+(1|IID),family=binomial,data=data_input)
fit1e=glmer(emerge_time50 ~ EL+WL+(1|IID),family=binomial,data=data_input)
fit1f=glmer(emerge_time50 ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit1g=glmer(emerge_time50 ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit1h=glmer(emerge_time50 ~ WL*sex+(1|IID),family=binomial,data=data_input)
fit1i=glmer(emerge_time50 ~ WL+sex+(1|IID),family=binomial,data=data_input)
fit1j=glmer(emerge_time50 ~ EL+(1|IID),family=binomial,data=data_input)
fit1k=glmer(emerge_time50 ~ WL+(1|IID),family=binomial,data=data_input)
fit1l=glmer(emerge_time50 ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g,fit1h,fit1i,fit1j,fit1k,fit1l)
Anova(fit1c,type="III")
Anova(fit1d,type="III")
Anova(fit1e,type="III")
Anova(fit1f,type="III")
Anova(fit1g,type="III")
Anova(fit1h,type="III")
Anova(fit1i,type="III")
Anova(fit1j,type="III")
Anova(fit1k,type="III") 
Anova(fit1l,type="III") 

### %MRWS ####
fit4a=glmer(emerged_T_cutoff ~ EL*relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4b=glmer(emerged_T_cutoff ~ EL*relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4c=glmer(emerged_T_cutoff ~ EL+relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4d=glmer(emerged_T_cutoff ~ EL*relMRWS+(1|IID),family=binomial,data=data_input)
fit4e=glmer(emerged_T_cutoff ~ EL+relMRWS+(1|IID),family=binomial,data=data_input)
fit4f=glmer(emerged_T_cutoff ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit4g=glmer(emerged_T_cutoff ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit4h=glmer(emerged_T_cutoff ~ relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit4i=glmer(emerged_T_cutoff ~ relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit4j=glmer(emerged_T_cutoff ~ EL+(1|IID),family=binomial,data=data_input)
fit4k=glmer(emerged_T_cutoff ~ relMRWS+(1|IID),family=binomial,data=data_input)
fit4l=glmer(emerged_T_cutoff ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit4a,fit4b,fit4c,fit4d,fit4e,fit4f,fit4g,fit4h,fit4i,fit4j,fit4k,fit4l)
Anova(fit4c,type="III")
Anova(fit4d,type="III")
Anova(fit4f,type="III")
Anova(fit4g,type="III")
Anova(fit4h,type="III")
Anova(fit4i,type="III")
Anova(fit4k,type="III") 
Anova(fit4l,type="III") 
plot(allEffects(fit4d), type='response') 
plot(allEffects(fit4h), type='response') 

# testing for linearity (on transformed scale) and testing for outliers/infl obs are not discussed in ABDA

# test for collinearity: check with variance inflation factors
# switch to glm, remove random factors, and remove interaction?
fit4h2=glm(emerge_time50 ~ relMRWS+sex,family=binomial,data=data_input)
vif(fit4h2) # vif < 5, no collinearity

# check for overdispersion
# include random factor for every observation
data_input$obs <- factor(1:nrow(data_input))
fit4h3=glmer(emerge_time50 ~ relMRWS*sex+(1|IID)+(1|obs),family=binomial,data=data_input)
# singular error, so is not needed. no overdispersion


## 'positionStart' ####
# positionStart = time as binary, t=0 -> 0, all others -> 1
# use all trials
data_input <- data_long_trial12345
# use all trials, but exclude beetles that emerge ~50% of the time
data_input <- subset(data_long_trial12345, emergeCount != 2 & emergeCount != 3)
# use all trials, but only include beetles with 'extreme' times 
data_input <- subset(data_long_trial12345, time < 50 | time > 250)
# use only trial 3, 4 and 5
data_input <- data_long_trial345
# use only trial 3, 4 and 5, and only include beetles with 'extreme' times 
data_input <- subset(data_long_trial345, time < 50 | time > 250)

# check trends by eye
xyplot(positionStart ~ WL, data=data_input ,type=c("p","r"))
xyplot(positionStart ~ relMRWS, data=data_input ,type=c("p","r"))
xyplot(positionStart ~ EL, data=data_input ,type=c("p","r"))
xyplot(positionStart ~ sex, data=data_input ,type=c("p","r"))

### WL ####
fit1a=glmer(positionStart ~ EL*WL*sex+(1|IID),family=binomial,data=data_input)
fit1b=glmer(positionStart ~ EL*WL+sex+(1|IID),family=binomial,data=data_input)
fit1c=glmer(positionStart ~ EL+WL+sex+(1|IID),family=binomial,data=data_input)
fit1d=glmer(positionStart ~ EL*WL+(1|IID),family=binomial,data=data_input)
fit1e=glmer(positionStart ~ EL+WL+(1|IID),family=binomial,data=data_input)
fit1f=glmer(positionStart ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit1g=glmer(positionStart ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit1h=glmer(positionStart ~ WL*sex+(1|IID),family=binomial,data=data_input)
fit1i=glmer(positionStart ~ WL+sex+(1|IID),family=binomial,data=data_input)
fit1j=glmer(positionStart ~ EL+(1|IID),family=binomial,data=data_input)
fit1k=glmer(positionStart ~ WL+(1|IID),family=binomial,data=data_input)
fit1l=glmer(positionStart ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit1a,fit1b,fit1c,fit1d,fit1e,fit1f,fit1g,fit1h,fit1i,fit1j,fit1k,fit1l)
Anova(fit1g,type="III")
Anova(fit1i,type="III")
Anova(fit1j,type="III")
Anova(fit1k,type="III") #nothing significant
Anova(fit1l,type="III") #nothing significant
plot(allEffects(fit1l, residuals=TRUE)) 
plot(allEffects(fit1l), type='response') 

### WEratio ####
data_input <- data_long_trial12345
fit2a=glmer(positionStart ~ EL*WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2b=glmer(positionStart ~ EL*WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2c=glmer(positionStart ~ EL+WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2d=glmer(positionStart ~ EL*WEratio+(1|IID),family=binomial,data=data_input)
fit2e=glmer(positionStart ~ EL+WEratio+(1|IID),family=binomial,data=data_input)
fit2f=glmer(positionStart ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit2g=glmer(positionStart ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit2h=glmer(positionStart ~ WEratio*sex+(1|IID),family=binomial,data=data_input)
fit2i=glmer(positionStart ~ WEratio+sex+(1|IID),family=binomial,data=data_input)
fit2j=glmer(positionStart ~ EL+(1|IID),family=binomial,data=data_input)
fit2k=glmer(positionStart ~ WEratio+(1|IID),family=binomial,data=data_input)
fit2l=glmer(positionStart ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit2a,fit2b,fit2c,fit2d,fit2e,fit2f,fit2g,fit2h,fit2i,fit2j,fit2k,fit2l)
Anova(fit2g,type="III")
Anova(fit2i,type="III")
Anova(fit2j,type="III")
Anova(fit2k,type="III") #nothing significant
Anova(fit2l,type="III")

### WEsurfaceratio ####
data_input <- data_long_trial12345
fit3a=glmer(positionStart ~ EL*WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3b=glmer(positionStart ~ EL*WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3c=glmer(positionStart ~ EL+WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3d=glmer(positionStart ~ EL*WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3e=glmer(positionStart ~ EL+WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3f=glmer(positionStart ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit3g=glmer(positionStart ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit3h=glmer(positionStart ~ WEsurfaceratio*sex+(1|IID),family=binomial,data=data_input)
fit3i=glmer(positionStart ~ WEsurfaceratio+sex+(1|IID),family=binomial,data=data_input)
fit3j=glmer(positionStart ~ EL+(1|IID),family=binomial,data=data_input)
fit3k=glmer(positionStart ~ WEsurfaceratio+(1|IID),family=binomial,data=data_input)
fit3l=glmer(positionStart ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit3a,fit3b,fit3c,fit3d,fit3e,fit3f,fit3g,fit3h,fit3i,fit3j,fit3k,fit3l)
Anova(fit3g,type="III")
Anova(fit3i,type="III")
Anova(fit3j,type="III")
Anova(fit3k,type="III") #nothing significant
Anova(fit3l,type="III")

### %MRWS ####
data_input <- data_long_trial345
fit3a=glmer(positionStart ~ EL*relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit3b=glmer(positionStart ~ EL*relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit3c=glmer(positionStart ~ EL+relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit3d=glmer(positionStart ~ EL*relMRWS+(1|IID),family=binomial,data=data_input)
fit3e=glmer(positionStart ~ EL+relMRWS+(1|IID),family=binomial,data=data_input)
fit3f=glmer(positionStart ~ EL*sex+(1|IID),family=binomial,data=data_input)
fit3g=glmer(positionStart ~ EL+sex+(1|IID),family=binomial,data=data_input)
fit3h=glmer(positionStart ~ relMRWS*sex+(1|IID),family=binomial,data=data_input)
fit3i=glmer(positionStart ~ relMRWS+sex+(1|IID),family=binomial,data=data_input)
fit3j=glmer(positionStart ~ EL+(1|IID),family=binomial,data=data_input)
fit3k=glmer(positionStart ~ relMRWS+(1|IID),family=binomial,data=data_input)
fit3l=glmer(positionStart ~ sex+(1|IID),family=binomial,data=data_input)
MuMIn::AICc(fit3a,fit3b,fit3c,fit3d,fit3e,fit3f,fit3g,fit3h,fit3i,fit3j,fit3k,fit3l)
Anova(fit3g,type="III")
Anova(fit3j,type="III")
Anova(fit3l,type="III")



# SURVIVAL ANALYSIS ####
# Tried to normalize time data by excluding 0 and 300s
# Also tried binomial GLM to compare 'emerging within 300s' and 'never emerging'
# Use survival analysis to get more info from data

# 0s data not possible in survival analysis
data_input <- subset(data_long_trial12345, time != 0)
# ALTERNATIVE: make it 0.01s
data_long_0.01 <- data_long_trial12345
data_long_0.01$time[data_long$time <= 0] <- 0.01
data_input <- data_long_0.01
# ALTERNATIVE: consider the experiment only started when the beetle is submerged
# therefor we have to ignore all t=0 data, and shift the start of the experiment to t=30s
data_input <- subset(data_long_trial12345, time >= 30)
data_input$time <- data_input$time - 29.9

# censoring column in data already added (see section 'load data')

## Parametric ####
# first decide on distribution by fitting null models
fit1plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'exp')
plot(fit1plot, type = "survival")
fit2plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'weibull')
plot(fit2plot, type = "survival")
fit3plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'gamma')
plot(fit3plot, type = "survival")
fit4plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'gompertz')
plot(fit4plot, type = "survival")
fit5plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'gengamma')
plot(fit5plot, type = "survival")
plot(fit5plot, type = "hazard")
fit6plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'genf')
plot(fit6plot, type = "survival")
plot(fit6plot, type = "hazard")
fit7plot <- flexsurvreg(Surv(time, censored,type = 'right') ~ 1,
                        data=data_input,
                        dist = 'llogis')
plot(fit7plot, type = "survival")
plot(fit7plot, type = "hazard")
MuMIn::AICc(fit1plot,fit2plot, fit3plot,fit4plot,fit5plot,fit6plot,fit7plot)


#llogis is the only distribution also available in survreg package to calculate Anova table
fit1 <- survreg(Surv(time, censored, type="right") ~  WL+EL+sex, dist = "loglogistic", data=data_input)
fit2 <- survreg(Surv(time, censored, type="right") ~  WL*EL+sex, dist = "loglogistic", data=data_input)
fit3 <- survreg(Surv(time, censored, type="right") ~  WL+EL*sex, dist = "loglogistic", data=data_input)
fit4 <- survreg(Surv(time, censored, type="right") ~  WL*EL, dist = "loglogistic", data=data_input)
fit5 <- survreg(Surv(time, censored, type="right") ~  WL+EL, dist = "loglogistic", data=data_input)
fit6 <- survreg(Surv(time, censored, type="right") ~  WL*sex, dist = "loglogistic", data=data_input)
fit7 <- survreg(Surv(time, censored, type="right") ~  WL+sex, dist = "loglogistic", data=data_input)
fit8 <- survreg(Surv(time, censored, type="right") ~  EL*sex, dist = "loglogistic", data=data_input)
fit9 <- survreg(Surv(time, censored, type="right") ~  EL+sex, dist = "loglogistic", data=data_input)
fit10 <- survreg(Surv(time, censored, type="right") ~  WL, dist = "loglogistic", data=data_input)
fit11 <- survreg(Surv(time, censored, type="right") ~  EL, dist = "loglogistic", data=data_input)
fit12 <- survreg(Surv(time, censored, type="right") ~  sex, dist = "loglogistic", data=data_input)
MuMIn::AICc(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12)
Anova(fit6, type = "III") 

# frailty model (for random effect of IID) not possible in survreg of flexsurvreg

## Non parametric ####
# calculates relative risk/hazard
# Is there a 'difference in risk' after treatment a or b
### WL ####
fit1 <- coxph(Surv(time, censored, type="right") ~ WL * EL * sex, data=data_input)
fit2 <- coxph(Surv(time, censored, type="right") ~ WL * EL + sex, data=data_input)
fit3 <- coxph(Surv(time, censored, type="right") ~ WL + EL * sex, data=data_input)
fit4 <- coxph(Surv(time, censored, type="right") ~ WL + EL + sex, data=data_input)
fit5 <- coxph(Surv(time, censored, type="right") ~ WL * EL , data=data_input)
fit6 <- coxph(Surv(time, censored, type="right") ~ WL + EL , data=data_input)
fit7 <- coxph(Surv(time, censored, type="right") ~ WL * sex, data=data_input)
fit8 <- coxph(Surv(time, censored, type="right") ~ WL + sex, data=data_input)
fit9 <- coxph(Surv(time, censored, type="right") ~ EL * sex, data=data_input)
fit10 <- coxph(Surv(time, censored, type="right") ~ EL + sex, data=data_input)
fit11 <- coxph(Surv(time, censored, type="right") ~ WL , data=data_input)
fit12 <- coxph(Surv(time, censored, type="right") ~ EL, data=data_input)
fit13 <- coxph(Surv(time, censored, type="right") ~ sex, data=data_input)
MuMIn::AICc(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13)
fit7
Anova(fit7,type='III')
autoplot(survfit(fit7), ylab = "submerged", xlab = "time (s)", 
         main = "emerging beetles")

# Manually compare WL below/above median

median_wl <- median(data_input$WL)
SW <- 3
LW <- 4
newdata <- expand.grid(
  WL = c(SW, LW),
  sex = levels(data_input$sex)
)

newdata$group <- paste0(newdata$sex,ifelse(newdata$WL < median_wl, " - 3 mm", " - 4 mm"))

#make separate survival curves per group
sf_df_list <- map(1:nrow(newdata), function(i) {
  sf_i <- survfit(fit7, newdata = newdata[i, , drop = FALSE])
  df_i <- fortify(sf_i)
  df_i$group <- newdata$group[i]
  df_i
})
#combine
sf_df_all <- bind_rows(sf_df_list)


data_plot <- subset(sf_df_all, group== 'F - 3 mm' | group== 'F - 4 mm')
plot_F <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "F - 3 mm" = "blue",  # orange
    "F - 4 mm" = "red"  # purple
    )) +
  scale_fill_manual(values = c(
    "F - 3 mm" = "blue",  # orange
    "F - 4 mm" = "red"  # purple
    )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL  
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8, 0.75))
plot_F

data_plot <- subset(sf_df_all, group== 'M - 3 mm' | group== 'M - 4 mm')
plot_M <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "M - 3 mm" = "blue",  # green
    "M - 4 mm" = "red"   # pinkish magenta (but not typical blue)
  )) +
  scale_fill_manual(values = c(
    "M - 3 mm" = "blue",  # green
    "M - 4 mm" = "red"   # pinkish magenta (but not typical blue)
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL
       ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8, 0.75))
plot_M

plot_MF = grid.arrange(plot_F, 
                     plot_M + rremove("ylab"), 
                     ncol=2) # single column (one figure below the other)

ggsave(
  plot_MF, 
  file="plot_survival_NonPara_WL.svg", width=9, height=4.5)


### %MRWS ####
fit1 <- coxph(Surv(time, censored, type="right") ~ relMRWS * EL * sex, data=data_input)
fit2 <- coxph(Surv(time, censored, type="right") ~ relMRWS * EL + sex, data=data_input)
fit3 <- coxph(Surv(time, censored, type="right") ~ relMRWS + EL * sex, data=data_input)
fit4 <- coxph(Surv(time, censored, type="right") ~ relMRWS + EL + sex, data=data_input)
fit5 <- coxph(Surv(time, censored, type="right") ~ relMRWS * EL , data=data_input)
fit6 <- coxph(Surv(time, censored, type="right") ~ relMRWS + EL , data=data_input)
fit7 <- coxph(Surv(time, censored, type="right") ~ relMRWS * sex, data=data_input)
fit8 <- coxph(Surv(time, censored, type="right") ~ relMRWS + sex, data=data_input)
fit9 <- coxph(Surv(time, censored, type="right") ~ EL * sex, data=data_input)
fit10 <- coxph(Surv(time, censored, type="right") ~ EL + sex, data=data_input)
fit11 <- coxph(Surv(time, censored, type="right") ~ relMRWS , data=data_input)
fit12 <- coxph(Surv(time, censored, type="right") ~ EL, data=data_input)
fit13 <- coxph(Surv(time, censored, type="right") ~ sex, data=data_input)
MuMIn::AICc(fit1,fit2,fit3,fit4,fit5,fit6,fit7,fit8,fit9,fit10,fit11,fit12,fit13)
fit7
Anova(fit7,type='III')
autoplot(survfit(fit7), ylab = "submerged", xlab = "time (s)", 
         main = "emerging beetles")

cox.zph(fit7)  # p>0.05 if proportional hazards assumption not violated
plot(cox.zph(fit7))

# check for repeated measurements with frailty model
fit7a <- coxph(Surv(time, censored, type="right") ~ relMRWS * sex + frailty(IID, dist="gamma"), data=data_input)
fit7b <- coxph(Surv(time, censored, type="right") ~ relMRWS * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit7a,fit7b)
fit7b
summary(fit7b)
summary(contrast(emmeans(fit7b, ~sex), method="pairwise", adjust="tukey"), type="response")


# 'Manually' compare relMRWS below/above median

median_relMRWS <- median(data_input$relMRWS)
SW <- 0.35
LW <- 0.60
newdata <- expand.grid(
  relMRWS = c(SW, LW),
  sex = levels(data_input$sex)
)

newdata$group <- paste0(newdata$sex,ifelse(newdata$relMRWS < median_relMRWS, " - 0.35", " - 0.60"))

#make separate survival curves per group
sf_df_list <- map(1:nrow(newdata), function(i) {
  sf_i <- survfit(fit7, newdata = newdata[i, , drop = FALSE])
  df_i <- fortify(sf_i)
  df_i$group <- newdata$group[i]
  df_i})
#combine
sf_df_all <- bind_rows(sf_df_list)


data_plot <- subset(sf_df_all, group== 'F - 0.35' | group== 'F - 0.60')
plot_F <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "F - 0.35" = "blue",  # orange
    "F - 0.60" = "red"  # purple
  )) +
  scale_fill_manual(values = c(
    "F - 0.35" = "blue",  # orange
    "F - 0.60" = "red"  # purple
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.15, 1)) + 
  labs(x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL  
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8, 0.75))
plot_F

data_plot <- subset(sf_df_all, group== 'M - 0.35' | group== 'M - 0.60')
plot_M <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "M - 0.35" = "blue",  # green
    "M - 0.60" = "red"   # pinkish magenta (but not typical blue)
  )) +
  scale_fill_manual(values = c(
    "M - 0.35" = "blue",  # green
    "M - 0.60" = "red"   # pinkish magenta (but not typical blue)
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.15, 1)) + 
  labs(x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.8, 0.75))
plot_M

plot_MF = grid.arrange(plot_F, 
                       plot_M + rremove("ylab"), 
                       ncol=2) # single column (one figure below the other)

ggsave(
  plot_MF, 
  file="plot_survival_NonPara_relMRWS.svg", width=9, height=4.5)


