# Sawangjit et al. (2021): 1 week remote novel object recognition (NOR) memory task
# Mixed model statistical analyses for primary outcome variables (discrimination ratio, change in mean rearing duration)
# By: Max Harkotte
# Contact: maximilian.harkotte@gmail.com
# Last update: June 2021

rm(list = ls()) # clear workspace
cat("\014") # clear console

# 0 - Load packages -------------------------------------------------------
library(dplyr)
library(ggpubr)

# 1 - Source file ---------------------------------------------------------
dataPath <- "Z:/Max/03_Sleep_vs_wake_consolidation/1wk_NOR_new/Behavior/Data"
setwd(dataPath)


# 2 - Read in data --------------------------------------------------------
Test_wide <-
  read.csv2(
    "2021-06-Behavior_data.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )

SleepArch <-
  read.csv2(
    "SleepArch.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = FALSE
  )

Spindles <-
  read.csv2(
    "Spindles_Behavior_EEGLeft.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )

SOs <-
  read.csv2(
    "SO_new.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )

Theta <-
  read.csv2(
    "Theta_Behavior_EEGLeft.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )


# 3 - Wrangling -----------------------------------------------------------
for (i in 2:(length(SleepArch))) {
  SleepArch[, i]  <- as.numeric(as.character(SleepArch[, i]))
}

for (i in 5:(length(Spindles))) {
  Spindles[, i]  <- as.numeric(as.character(Spindles[, i]))
}

Spindles_Mus <- subset(Spindles, Spindles$Drug == "MUS")
Spindles_Sal <- subset(Spindles, Spindles$Drug == "SAL")
Spindles_Mus <- subset(Spindles_Mus, Spindles_Mus$Animal != "05-P10-S1")
Spindles = merge(Spindles_Mus, Spindles_Sal, by.x = "Animal", by.y = "Animal") # x = mus, y = sal

for (i in 4:(length(SOs))) {
  SOs[, i]  <- as.numeric(as.character(SOs[, i]))
}

SOs_Mus <- subset(SOs, SOs$Drug == "MUS")
SOs_Sal <- subset(SOs, SOs$Drug == "SAL")
SOs_Mus <- subset(SOs_Mus, SOs_Mus$Animal != "05-P10-S1")
SOs = merge(SOs_Mus, SOs_Sal, by.x = "Animal", by.y = "Animal") # x = mus, y = sal

Theta$Animal <- Theta$Animal.y
Theta <- select(Theta, -Animal.y)
Theta <- select(Theta, Drug, everything())

for (i in 4:(length(Theta))) {
  Theta[, i]  <- as.numeric(as.character(Theta[, i]))
}

Theta_Mus <- subset(Theta, Theta$Drug == "MUS")
Theta_Sal <- subset(Theta, Theta$Drug == "SAL")
Theta_Mus <- subset(Theta_Mus, Theta_Mus$Animal != "05-P10-S1")
Theta = merge(Theta_Mus, Theta_Sal, by.x = "Animal", by.y = "Animal") # x = mus, y = sal

Memory_performance <- subset(Theta, select = c(Animal, DR_1.x, DR_1.y, dRear_1.x, dRear_1.y))
SleepArch <- SleepArch[!(SleepArch$ï..animal=="05-P10-BL1"),]

SleepArch <- merge(Memory_performance, SleepArch, by.x = "Animal", by.y = "ï..animal")


# 4 - Comparisons ---------------------------------------------------------
t.test(
  SOs$mean_Pk2Pk_Amp.y,
  SOs$mean_Pk2Pk_Amp.x,
  alternative = "two.sided",
  paired = TRUE
)

mean(SOs$mean_Pk2Pk_Amp.x, na.rm = TRUE)
sd(SOs$mean_Pk2Pk_Amp.x, na.rm = TRUE)

cor.test(SleepArch$dRear_1.x, SleepArch$Rtimemus, method="spearman")


