# Sawangjit et al. (2021): 1 week remote novel object recognition (NOR) memory task
# Mixed model statistical analyses for primary outcome variables (discrimination ratio, change in mean rearing duration)
# By: Max Harkotte
# Contact: maximilian.harkotte@gmail.com
# Last update: June 2021

rm(list = ls()) # clear workspace
cat("\014") # clear console

# 0 - Load packages -------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(psych)
library(lme4)
library(effects)
library(sjPlot)
library(ggpubr)
library(rstatix)
library(lmerTest)
library(stats)

# 1 - Source file ---------------------------------------------------------
dataPath <- "Z:/Max/1wk_NOR_new/Behavior/Data"
setwd(dataPath)

# Add aestetics for plotting
limits = aes(ymax = mean + (se), ymin = mean - (se)) # (1.96*se) for confidence intervals
dodge = position_dodge(width = 0.8)

# 2 - Read in data --------------------------------------------------------

# Test wide contains Discrimination Ratios (dRear) for three time intervals
# of the memory test: 0-1 min (DR1), 0-3 min (DR3), 0-5 min (DR5) (Minute)
# Between Factor: Sleep or Wake during Post-Encoding Interval (PostInterval)
# Mixed Factor: Muscimol or Saline infusion in the dorsal hippocampus during
# Post-Encoding Interval (Drug)

Test_wide <-
  read.csv2(
    "context_change.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = FALSE
  )


Test_wide$Context[Test_wide$Context == 1] = "OPTO"
Test_wide$Context[Test_wide$Context == 2] = "SURG"
Test_wide$Context = as.factor(Test_wide$Context)

Test_wide$Group[Test_wide$Group == 1] = "SLEEP"
Test_wide$Group[Test_wide$Group == 2] = "WAKE"
Test_wide$Group = as.factor(Test_wide$Group)

for (i in 4:(length(Test_wide))) {
  Test_wide[, i]  <- as.numeric(as.character(Test_wide[, i]))
}

Test_wide <- Test_wide %>%
  rename(
    PostInterval = Group,
    animal =  ï..animal,
    dRear_mean_dur_1 = perchangerearmeandur1,
    dRear_mean_dur_3 = perchangerearmeandur3
  )


## DR 1
DR1_Sleep <- DR_long[which(DR_long$Minute == 'DR1' &
                                 DR_long$PostInterval == 'SLEEP'), ]

DR3_Sleep <- DR_long[which(DR_long$Minute == 'DR3' &
                             DR_long$PostInterval == 'SLEEP'), ]

DR1_Wake <- DR_long[which(DR_long$Minute == 'DR1' &
                             DR_long$PostInterval == 'WAKE'), ]

DR3_Wake <- DR_long[which(DR_long$Minute == 'DR3' &
                             DR_long$PostInterval == 'WAKE'), ]

t.test(DR1_Sleep$DR,
       DR1_Wake$DR,
       alternative = "two.sided",
       paired = FALSE)

t.test(DR1_Sleep$DR, mu = 0)
t.test(DR1_Sleep$DR, mu = 0, alternative = "greater")
describe(DR1_Sleep$DR)

t.test(DR3_Sleep$DR, mu = 0)
t.test(DR3_Sleep$DR, mu = 0, alternative = "greater") 
describe(DR3_Sleep$DR)

t.test(DR1_Wake$DR, mu = 0)
t.test(DR1_Wake$DR, mu = 0, alternative = "greater")
describe(DR1_Wake$DR)

t.test(DR3_Wake$DR, mu = 0)
t.test(DR3_Wake$DR, mu = 0, alternative = "greater")
describe(DR3_Wake$DR)



# Discrimination Ratio  ---------------------------------------------------


test_DR <-
  Test_wide[c("animal", "PostInterval", "Context" , "DR1", "DR3")]
DR_long <-
  pivot_longer(test_DR,
               cols = 4:5 ,
               names_to = "Minute",
               values_to = "DR")

DR_long$Minute <- as.factor(DR_long$Minute)

## Data exploration

# Data summary
DR_long %>%
  group_by(PostInterval, Minute) %>%
  get_summary_stats(DR, type = "mean_sd")

bxp_interval <- ggboxplot(
  DR_long,
  x = "Minute",
  y = "DR",
  palette = "jco",
  add = "jitter",
  facet.by = "PostInterval"
)
bxp_interval

# Basic model without interactions
basic_hlm <-
  lmer(DR ~ Minute + PostInterval + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Add PostInterval*Minute
Minute_Sleep_interaction_hlm <-
  lmer(DR ~ Minute * PostInterval + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Remove Postinterval
rm_PostInterval_hlm <-
  lmer(DR ~ Minute + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Remote Minute
rm_Minute_hlm <-
  lmer(DR ~ PostInterval + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Test nested interactions
anova(basic_hlm, Minute_Sleep_interaction_hlm) # Test for PostInterval*Minute interaction
anova(rm_PostInterval_hlm, basic_hlm) # Test for Post Interval
anova(rm_Minute_hlm, basic_hlm)# Test for overall interaction


# Rearing -----------------------------------------------------------------

dRear_wide <-
  Test_wide[c("animal",
              "PostInterval",
              "Context",
              "dRear_mean_dur_1",
              "dRear_mean_dur_3")]

# Wide to long data format
dRear_long <-
  pivot_longer(
    dRear_wide,
    cols = 4:5 ,
    names_to = "Minute",
    values_to = "dRear"
  )

dRear_long$Minute <- as.factor(dRear_long$Minute)


## Data exploration

# Data summary
dRear_long %>%
  group_by(PostInterval, Drug, Minute) %>%
  get_summary_stats(dRear, type = "mean_se")


# Visualization
bxp_interval <- ggboxplot(
  dRear_long,
  x = "Minute",
  y = "dRear",
  palette = "jco",
  add = "jitter",
  facet.by = "PostInterval"
)
bxp_interval


# 6 - Statistics ----------------------------------------------------------
# Omnibus mixed models ----------------------------------------------------

# Basic model without interactions
basic_hlm <-
  lmer(dRear ~ Minute + PostInterval  + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Add PostInterval*Minute
Minute_Sleep_interaction_hlm <-
  lmer(dRear ~ Minute * PostInterval + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# rm PI
rm_PostInterval_hlm <-
  lmer(dRear ~ Minute  + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# rm PI
rm_Minute_hlm <-
  lmer(dRear ~ PostInterval  + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Test nested interactions
anova(basic_hlm, Minute_Sleep_interaction_hlm) # Test for PostInterval*Minute interaction
anova(rm_PostInterval_hlm, basic_hlm) # Test for PostInterval
anova(rm_Minute_hlm, basic_hlm) # Test for Minute


