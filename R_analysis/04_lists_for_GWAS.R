# collapse all: alt+cmd+O (mac)
# Initialise r ####
# loading packages, working directory and set sum contrast 

library(openxlsx)    # Read/write Excel files
library(dplyr)       # Data manipulation and transformation (in tidyverse)
library(afex)        # ANOVA and mixed models for factorial experiments
library(reshape2)    # Data reshaping (wide <-> long format)



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

# phenotypes ####
write.xlsx(data[, c("IID","WL","EL","relMRWS","action5trials",'action3trials')], file = "phenotype.xlsx", rowNames = FALSE)


# average time and divide beetles in two groups (binary)
data$time_average <- rowMeans(data[, c("time1", "time2", "time3", "time4", "time5")], na.rm = TRUE)
hist(data$time_average)

data$time_binary <- NA  # initialize with NA
data$time_binary[data$time_average < 100] <- 0
data$time_binary[data$time_average > 200] <- 1
hist(data$time_binary)
length(na.omit(data$time_binary))

# make phenotype file
# Filter data to include only IID with time_binary 0 or 1
filtered_data <- data[data$time_binary %in% c(0,1), ]

# Function to convert IID like Pc25Np002 to filename like S-02.BarSW.filtered.sorted.bam
convert_iid_to_filename <- function(iid) {
  # Extract the number part
  num <- sub("Pc25Np0*", "", iid)  # remove prefix and leading zeros
  num <- as.numeric(num)
  # Format number with leading zero if less than 10
  num_str <- sprintf("%02d", num)
  paste0("S-", num_str, ".BarSW.filtered.sorted.bam")
}

# Apply conversion
filtered_data$filename <- sapply(filtered_data$IID, convert_iid_to_filename)

# Prepare output dataframe with 3 columns
output_df <- data.frame(
  filename1 = filtered_data$filename,
  filename2 = filtered_data$filename,
  time_binary = filtered_data$time_binary
)
sum(output_df$time_binary)
# Write to a tab-delimited file
write.table(output_df, "phenotype_time_binary.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# sample list ####
# export txt file with sample names of filtered beetles
filtered_iids <- NA
filtered_iids <- data$IID[data$time_binary %in% c(0,1)]
length(filtered_iids)
# Extract the numeric part from each IID
numbers <- sub("Pc25Np0*", "", filtered_iids)  # remove "Pc25Np" and leading zeros
# Pad with leading zero if needed (to always have 2 digits)
numbers_padded <- sprintf("%02d", as.numeric(numbers))
# Build filenames
filenames <- paste0("S-", numbers_padded, ".BarSW.filtered.sorted.bam")
filenames
# Export only those IIDs
write.table(filenames, "IID_before100after200.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)



# covariates ####
#export covariates file for gemma as df(1, sex)
sex <- as.numeric(data$sex) - 1
df <- data.frame(1,sex)[0:192,]
write.table(df,"covariates.txt",sep = '\t',
            row.names = FALSE,
            col.names = FALSE)