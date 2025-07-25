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
packageVersion("survival")
library(flexsurv)    # Flexible parametric survival models
library(scales)      # for percent_format()
library(ggfortify)   # Plot survival analysis
library(ggpubr)      # in combined plots, ex rremove("ylab")
library(readxl)


# Set working directory
dir <- "/Users/Sander/Documents/Universiteit/KU Leuven/Masterproef/"
setwd(dir)

# to compare levels (when using factors) to the grand mean and coefficients sum to zero
set_sum_contrasts()

# Load data ####
#data <- read.xlsx("CLEANED.xlsx")
data <- read.xlsx("NP25_phenotypes/Pogonus_Nieuwpoort_2025_local.xlsx")

for (i in 1:5) {
  action_col <- paste0("action", i)
  time_col <- paste0("time", i)
  new_col <- paste0("actionInduced", i)
  
  data[[new_col]] <- ifelse(
    data[[time_col]] == 0, "emerge",
    ifelse(data[[time_col]] == 300, "submerge", data[[action_col]])
  )
}

# calculate ratio of emerge/submerge ation 
data$action5trials <- apply(data[, c("actionInduced1", "actionInduced2", "actionInduced3", "actionInduced4", "actionInduced5")], 1, function(row) {
  count_emerge <- sum(row == "emerge", na.rm = TRUE)
  non_na_count <- sum(!is.na(row))
  divisor <- if (non_na_count == 5) 5 else 4
  count_emerge / divisor
})

data$action3trials <- apply(data[, c("actionInduced3", "actionInduced4", "actionInduced5")], 1, function(row) {
  count_emerge <- sum(row == "emerge", na.rm = TRUE)
  non_na_count <- sum(!is.na(row))
  divisor <- if (non_na_count == 3) 3 else 2
  count_emerge / divisor
})

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

data$sex_numeric <- data$sex
data$sex_numeric <- ifelse(data$sex=='F',0,1)

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

# SURVIVAL ANALYSIS ####
# Tried to normalize time data by excluding 0 and 300s
# Also tried binomial GLM to compare 'emerging within 300s' and 'never emerging'
# Use survival analysis to get more info from data
data_input <- data_long_trial345
data_input <- data_long_trial12345


# 0s data not possible in paramteric survival analysis, but fine in non-para?
data_input <- subset(data_long_trial12345, time != 0)
# ALTERNATIVE: make it 0.01s
data_long_0.01 <- data_long_trial12345
data_long_0.01$time[data_long$time <= 0] <- 0.01
data_input <- data_long_0.01 # 945 observations
# ALTERNATIVE 2: consider the experiment only started when the beetle is submerged
# therefor we have to ignore all t=0 data, and shift the start of the experiment to t=30s
data_input <- subset(data_long_trial12345, time >= 30) #740 observations
data_input$time <- data_input$time - 29.9

# censoring column in data already added (see section 'load data')

## PARAMETRIC ####
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

## NON PARAMETRIC ####
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
fit1
Anova(fit7,type='III')
autoplot(survfit(fit7), ylab = "submerged", xlab = "time (s)", 
         main = "emerging beetles")

# frailty
fit1a <- coxph(Surv(time, censored, type="right") ~ WL * EL * sex + frailty(IID, dist="gamma"), data=data_input)
fit1b <- coxph(Surv(time, censored, type="right") ~ WL * EL * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit1a,fit1b)
fit1b

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
fit13
Anova(fit1,type='III')
Anova(fit3,type='III')
Anova(fit8,type='III')
Anova(fit9,type='III')
Anova(fit10,type='III')
Anova(fit13,type='III')
summary(fit9)
summary(fit13)
autoplot(survfit(fit9), ylab = "submerged", xlab = "time (s)", 
         main = "emerging beetles")

cox.zph(fit13)  # p>0.05 if proportional hazards assumption not violated
plot(cox.zph(fit13))

# check for repeated measurements with frailty model
fit1a <- coxph(Surv(time, censored, type="right") ~ relMRWS * EL * sex + frailty(IID, dist="gamma"), data=data_input)
fit1b <- coxph(Surv(time, censored, type="right") ~ relMRWS * EL * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit1a,fit1b)
fit1b

fit3a <- coxph(Surv(time, censored, type="right") ~ relMRWS + EL * sex + frailty(IID, dist="gamma"), data=data_input)
fit3b <- coxph(Surv(time, censored, type="right") ~ relMRWS + EL * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit9a,fit9b)
fit3b
summary(fit3b)
summary(contrast(emmeans(fit3b, ~sex), method="pairwise", adjust="tukey"), type="response")

fit7a <- coxph(Surv(time, censored, type="right") ~ relMRWS * sex + frailty(IID, dist="gamma"), data=data_input)
fit7b <- coxph(Surv(time, censored, type="right") ~ relMRWS * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit7a,fit7b)
fit7b

fit9a <- coxph(Surv(time, censored, type="right") ~ EL * sex + frailty(IID, dist="gamma"), data=data_input)
fit9b <- coxph(Surv(time, censored, type="right") ~ EL * sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit9a,fit9b)
fit9b
summary(fit9b)
summary(contrast(emmeans(fit9b, ~sex), method="pairwise", adjust="tukey"), type="response")

fit13a <- coxph(Surv(time, censored, type="right") ~ sex + frailty(IID, dist="gamma"), data=data_input)
fit13b <- coxph(Surv(time, censored, type="right") ~ sex + frailty(IID, dist="gaussian"), data=data_input)
MuMIn::AICc(fit13a,fit13b)
fit13b
summary(contrast(emmeans(fit13b, ~sex), method="pairwise", adjust="tukey"), type="response")
anova(fit13, fit13b)

# plot female male curves. cannot plot frailty model, so use previous model for visualisation
# Fix EL at median
median_EL <- median(data_input$EL, na.rm = TRUE)

# New data for predictions with median EL, separately for sex
newdata_F <- data.frame(EL = median_EL, sex = "F")
newdata_M <- data.frame(EL = median_EL, sex = "M")

# Get survival fits for both sexes at median EL
sf_F <- survfit(fit13, newdata = newdata_F)
sf_M <- survfit(fit13, newdata = newdata_M)

# Fortify survival objects
df_F <- fortify(sf_F)
df_F$sex <- "Female"

df_M <- fortify(sf_M)
df_M$sex <- "Male"

# Combine data frames
sf_df <- rbind(df_F, df_M)

# Plot survival curves
plot_survival <- ggplot(sf_df, aes(x = time, y = surv, color = sex, fill = sex)) +
  geom_step(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Female" = "orange", "Male" = "turquoise")) +
  scale_fill_manual(values = c("Female" = "orange", "Male" = "turquoise")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = NULL,
    x = "Time (s)",
    y = "Submerged beetles (%)",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = c(0.75, 0.85),plot.title = element_text(hjust = 0.5,size=16))
plot_survival
ggsave(
  plot_survival, 
  file="NP25_phenotypes/plot_survival_NonPara_AllData_sex.svg", width=4.5, height=3.5)

# 'Manually' compare relMRWS below/above median

median_relMRWS <- median(data_input$relMRWS)
SW <- 0.25
LW <- 0.50
newdata <- expand.grid(
  relMRWS = c(SW, LW),
  sex = levels(data_input$sex)
)

newdata$group <- ifelse(newdata$relMRWS < median_relMRWS, "25% MRWS", "50% MRWS")

data_long_trial345 %>%
  group_by(sex) %>%
  summarise(mean_time = mean(time, na.rm = TRUE))

#make separate survival curves per group
sf_df_list <- map(1:nrow(newdata), function(i) {
  sf_i <- survfit(fit7, newdata = newdata[i, , drop = FALSE])
  df_i <- fortify(sf_i)
  df_i$sex <- newdata$sex[i]  
  df_i$group <- newdata$group[i]
  df_i})
#combine
sf_df_all <- bind_rows(sf_df_list)


data_plot <- subset(sf_df_all, sex == "F")
plot_F <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "25% MRWS" = "blue",  
    "50% MRWS" = "red"
  )) +
  scale_fill_manual(values = c(
    "25% MRWS" = "blue",
    "50% MRWS" = "red"
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.1, 1)) + 
  labs(title='Females',
       x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL  
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.75, 0.85),plot.title = element_text(hjust = 0.5,size=18))
plot_F


data_plot <- subset(sf_df_all, sex == "M")
plot_M <- ggplot(data_plot, aes(x = time, y = surv, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "25% MRWS" = "blue",  
    "50% MRWS" = "red"
    )) +
  scale_fill_manual(values = c(
    "25% MRWS" = "blue",  
    "50% MRWS" = "red"
    )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0.1, 1)) + 
  labs(title='Males',
       x = "Time (s)",
       y = "Submerged beetles (%)",
       color = NULL,
       fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.75, 0.85),plot.title = element_text(hjust = 0.5,size=18))
plot_M

plot_MF = grid.arrange(plot_F, 
                       plot_M + rremove("ylab"), 
                       ncol=2) # single column (one figure below the other)

ggsave(
  plot_MF, 
  file="results/plot_survival_NonPara_relMRWS.svg", width=9, height=4.5)


