# Sawangjit et al. (2022)
# Mixed model statistical analyses for primary outcome variables (discrimination ratio, change in mean rearing duration)
# By: Max Harkotte
# Contact: maximilian.harkotte@gmail.com
# Last update: July 2022

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
dataPath <- "Z:/Max/03_Sleep_vs_wake_consolidation/1wk_NOR_new/Behavior/Data"
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
    "2021-06-Behavior_data.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = TRUE
  )

Test_wide <-
  subset(Test_wide,!(animal %in% c("07-CAN-G1"))) # exclude animals from analysis

# 3 - Discrimination ratio ------------------------------------------------

test_DR <-
  Test_wide[c("animal", "PostInterval", "Drug", "DR1", "DR3")]

for (i in 4:(length(test_DR))) {
  test_DR[, i]  <- as.numeric(as.character(test_DR[, i]))
}

# Wide to long data format
DR_long <-
  pivot_longer(test_DR,
               cols = 4:5 ,
               names_to = "Minute",
               values_to = "DR")

DR_long$Minute <- as.factor(DR_long$Minute)

## Data exploration

# Data summary
DR_long %>%
  group_by(PostInterval, Drug, Minute) %>%
  get_summary_stats(DR, type = "mean_sd")


# Visualization
bxp_interval <- ggboxplot(
  DR_long,
  x = "Minute",
  y = "DR",
  color = "Drug",
  palette = "jco",
  add = "jitter",
  shape = "Drug",
  facet.by = "PostInterval"
)
bxp_interval

bxp_min <- ggboxplot(
  DR_long,
  x = "PostInterval",
  y = "DR",
  color = "Drug",
  palette = "jco",
  add = "jitter",
  shape = "Drug",
  facet.by = "Minute"
)
bxp_min

bxp_drug <- ggboxplot(
  DR_long,
  x = "Minute",
  y = "DR",
  color = "Drug",
  palette = "jco",
  add = "jitter",
  shape = "Drug",
  facet.by = "PostInterval"
)
bxp_drug

# 4 - Statistics ----------------------------------------------------------
# Omnibus mixed models ----------------------------------------------------

# Basic model without interactions
basic_hlm <-
  lmer(DR ~ Minute + PostInterval + Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Add Drug*PostInterval interaction
Drug_Sleep_interaction_hlm <-
  lmer(DR ~ Minute + PostInterval * Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Add PostInterval*Minute
Minute_Sleep_interaction_hlm <-
  lmer(DR ~ Minute * PostInterval + Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Add Drug*Minute
Minute_Drug_interaction_hlm <-
  lmer(DR ~ Minute * Drug + PostInterval + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Three way mixed model with all possible interactions
complete_hlm <-
  lmer(DR ~ Minute * PostInterval * Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Test nested interactions
anova(basic_hlm, Drug_Sleep_interaction_hlm) # Test for PostInterval*Drug interaction
anova(basic_hlm, Minute_Sleep_interaction_hlm) # Test for PostInterval*Minute interaction
anova(basic_hlm, Minute_Drug_interaction_hlm) # Test for Drug*Minute interaction
anova(Drug_Sleep_interaction_hlm, complete_hlm) # Test for overall interaction

# Summary for selected model after interactions were tested
summary(Drug_Sleep_interaction_hlm)

# Remove factor Minute
rm_minute_hlm <-
  lmer(DR ~ PostInterval * Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Remove factor PostInterval
rm_interval_hlm <-
  lmer(DR ~ Minute + Drug + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Remove factor Drug
rm_drug_hlm <-
  lmer(DR ~ Minute + PostInterval + (1 | animal),
       data = DR_long,
       REML = FALSE)

# Test nested models for single factors
anova(rm_minute_hlm, Drug_Sleep_interaction_hlm) # Test for Factor Minute
anova(rm_interval_hlm, Drug_Sleep_interaction_hlm) # Test for Factor PostInterval
anova(rm_drug_hlm, Drug_Sleep_interaction_hlm) # Test for Factor Drug


# Check assumptions of selected model
plot_model(Drug_Sleep_interaction_hlm, type = 'diag')
tab_model(Drug_Sleep_interaction_hlm, show.intercept = FALSE)

# Post-hoc tests for Post-Interval*Drug Interaction ------------------------
# Post hoc analyses include two-way mixed models with animals as random
# effect and factors PostInterval/Drug and Minute for four different groups:
# 1. Sleep group
# 2. Wake groups
# 3. Muscimol infusions
# 4. Saline infusions


# Post hoc Sleep condition ------------------------------------------------
df_sleep <- subset(DR_long, DR_long$PostInterval == "Sleep")

basic_sleep_hlm <-
  lmer(DR ~ Minute + Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

rm_drug_sleep_hlm <-
  lmer(DR ~ Minute + (1 | animal),
       data = df_sleep,
       REML = FALSE)

rm_min_sleep_hlm <-
  lmer(DR ~ Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

complete_sleep_hlm <-
  lmer(DR ~ Minute * Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

# Post-hoc sleep group
anova(rm_drug_sleep_hlm, basic_sleep_hlm) # test factor drug
anova(rm_min_sleep_hlm, basic_sleep_hlm) # test factor minute
anova(basic_sleep_hlm, complete_sleep_hlm) # test Drug*Minute interaction

tab_model(complete_sleep_hlm, show.intercept = FALSE)

# Post hoc Wake condition -------------------------------------------------
df_wake <- subset(DR_long, DR_long$PostInterval == "Wake")

basic_wake_hlm <-
  lmer(DR ~ Minute + Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

rm_drug_wake_hlm <-
  lmer(DR ~ Minute + (1 | animal),
       data = df_wake,
       REML = FALSE)

rm_min_wake_hlm <-
  lmer(DR ~ Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

complete_wake_hlm <-
  lmer(DR ~ Minute * Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

# Post-hoc wake group
anova(rm_drug_wake_hlm, basic_wake_hlm) # test factor drug
anova(rm_min_wake_hlm, basic_wake_hlm) # test factor minute
anova(basic_wake_hlm, complete_wake_hlm) # test Drug*Minute interaction

tab_model(complete_wake_hlm, show.intercept = FALSE)

# Post hoc Muscimol condition ---------------------------------------------
df_mus <- subset(DR_long, DR_long$Drug == "Muscimol")

basic_mus_hlm <-
  lmer(DR ~ Minute + PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

rm_interval_mus_hlm <-
  lmer(DR ~ Minute + (1 | animal),
       data = df_mus,
       REML = FALSE)

rm_min_mus_hlm <-
  lmer(DR ~ PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

complete_mus_hlm <-
  lmer(DR ~ Minute * PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

# Post-hoc mus group
anova(rm_interval_mus_hlm, basic_mus_hlm) # test factor PostInterval
anova(rm_min_mus_hlm, basic_mus_hlm) # test factor minute
anova(basic_mus_hlm, complete_mus_hlm) # test PostInterval*Minute interaction

tab_model(complete_mus_hlm, show.intercept = FALSE)


# Post hoc Saline condition -----------------------------------------------
df_sal <- subset(DR_long, DR_long$Drug == "Saline")

basic_sal_hlm <-
  lmer(DR ~ Minute + PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

rm_interval_sal_hlm <-
  lmer(DR ~ Minute + (1 | animal),
       data = df_sal,
       REML = FALSE)

rm_min_sal_hlm <-
  lmer(DR ~ PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

complete_sal_hlm <-
  lmer(DR ~ Minute * PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

# Post-hoc sal group
anova(rm_interval_sal_hlm, basic_sal_hlm) # test factor PostInterval
anova(rm_min_sal_hlm, basic_sal_hlm) # test factor minute
anova(basic_sal_hlm, complete_sal_hlm) # test PostInterval*Minute interaction

tab_model(complete_sal_hlm, show.intercept = FALSE)

# Multiple pairwise comparisons -------------------------------------------

## DR 1
DR1_Sleep_Mus <- DR_long[which(DR_long$Minute == 'DR1' &
                                 DR_long$PostInterval == 'Sleep' &
                                 DR_long$Drug == "Muscimol"), ]

DR1_Sleep_Sal <- DR_long[which(DR_long$Minute == 'DR1'
                               &
                                 DR_long$PostInterval == 'Sleep' &
                                 DR_long$Drug == "Saline"), ]

DR1_Wake_Mus <- DR_long[which(DR_long$Minute == 'DR1'
                              &
                                DR_long$PostInterval == 'Wake' &
                                DR_long$Drug == "Muscimol"), ]

DR1_Wake_Sal <- DR_long[which(DR_long$Minute == 'DR1'
                              &
                                DR_long$PostInterval == 'Wake' &
                                DR_long$Drug == "Saline"), ]

DR3_Sleep_Mus <- DR_long[which(DR_long$Minute == 'DR3' &
                                 DR_long$PostInterval == 'Sleep' &
                                 DR_long$Drug == "Muscimol"), ]

DR3_Sleep_Sal <- DR_long[which(DR_long$Minute == 'DR3'
                               &
                                 DR_long$PostInterval == 'Sleep' &
                                 DR_long$Drug == "Saline"), ]

DR3_Wake_Mus <- DR_long[which(DR_long$Minute == 'DR3'
                              &
                                DR_long$PostInterval == 'Wake' &
                                DR_long$Drug == "Muscimol"), ]

DR3_Wake_Sal <- DR_long[which(DR_long$Minute == 'DR3'
                              &
                                DR_long$PostInterval == 'Wake' &
                                DR_long$Drug == "Saline"), ]


t.test(DR1_Sleep_Mus$DR,
       DR1_Sleep_Sal$DR,
       alternative = "two.sided",
       paired = TRUE)

t.test(DR1_Wake_Mus$DR,
       DR1_Sleep_Sal$DR,
       alternative = "two.sided",
       paired = FALSE)

t.test(DR3_Sleep_Mus$DR,
       DR3_Sleep_Sal$DR,
       alternative = "two.sided",
       paired = TRUE)

var.test(DR1_Sleep_Sal$DR, DR1_Wake_Sal$DR)

t.test(
  DR1_Sleep_Sal$DR,
  DR1_Wake_Sal$DR,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)

var.test(DR3_Sleep_Sal$DR, DR3_Wake_Sal$DR)

t.test(
  DR3_Sleep_Sal$DR,
  DR3_Wake_Sal$DR,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)

t.test(DR1_Sleep_Mus$DR,
       DR1_Wake_Mus$DR,
       alternative = "two.sided",
       paired = FALSE)

t.test(DR1_Sleep_Mus$DR,
       DR1_Wake_Sal$DR,
       alternative = "two.sided",
       paired = FALSE)

var.test(DR1_Wake_Mus$DR, DR1_Wake_Sal$DR)

t.test(
  DR1_Wake_Mus$DR,
  DR1_Wake_Sal$DR,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)

var.test(DR3_Wake_Mus$DR, DR3_Wake_Sal$DR)

t.test(
  DR3_Wake_Mus$DR,
  DR3_Wake_Sal$DR,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)


t.test(DR1_Sleep_Mus$DR, mu = 0)
t.test(DR1_Sleep_Mus$DR, mu = 0, alternative = "greater")
describe(DR1_Sleep_Mus$DR)
t.test(DR1_Sleep_Sal$DR, mu = 0)
t.test(DR1_Sleep_Sal$DR, mu = 0, alternative = "greater")
describe(DR1_Sleep_Sal$DR)

t.test(DR3_Sleep_Mus$DR, mu = 0)
t.test(DR3_Sleep_Mus$DR, mu = 0, alternative = "greater") # gets signif. 
t.test(DR3_Sleep_Sal$DR, mu = 0)
t.test(DR3_Sleep_Sal$DR, mu = 0, alternative = "greater")
describe(DR3_Sleep_Sal$DR)

t.test(DR1_Wake_Mus$DR, mu = 0)
t.test(DR1_Wake_Mus$DR, mu = 0, alternative = "greater")
t.test(DR1_Wake_Sal$DR, mu = 0)
t.test(DR1_Wake_Sal$DR, mu = 0, alternative = "greater")
describe(DR1_Wake_Sal$DR)

t.test(DR3_Wake_Mus$DR, mu = 0)
t.test(DR3_Wake_Mus$DR, mu = 0, alternative = "greater")
t.test(DR3_Wake_Sal$DR, mu = 0)
t.test(DR3_Wake_Sal$DR, mu = 0, alternative = "greater") # gets signif. 
describe(DR3_Wake_Sal$DR)


# 5 - Change in Rearing mean duration  ------------------------------------

dRear_wide <-
  Test_wide[c("animal",
              "PostInterval",
              "Drug",
              "dRear_mean_dur_1",
              "dRear_mean_dur_3")]

for (i in 4:(length(dRear_wide))) {
  dRear_wide[, i]  <- as.numeric(as.character(dRear_wide[, i]))
}

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
  color = "Drug",
  palette = "jco",
  add = "jitter",
  shape = "Drug",
  facet.by = "PostInterval"
)
bxp_interval

bxp_min <- ggboxplot(
  dRear_long,
  x = "PostInterval",
  y = "dRear",
  color = "Drug",
  palette = "jco",
  add = "jitter",
  shape = "Drug",
  facet.by = "Minute"
)
bxp_min

# 6 - Statistics ----------------------------------------------------------
# Omnibus mixed models ----------------------------------------------------

# Basic model without interactions
basic_hlm <-
  lmer(dRear ~ Minute + PostInterval + Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Add Drug*PostInterval interaction
Drug_Sleep_interaction_hlm <-
  lmer(dRear ~ Minute + PostInterval * Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Add PostInterval*Minute
Minute_Sleep_interaction_hlm <-
  lmer(dRear ~ Minute * PostInterval + Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Three way mixed model with all possible interactions
complete_hlm <-
  lmer(dRear ~ Minute * PostInterval * Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Test nested interactions
anova(basic_hlm, Drug_Sleep_interaction_hlm) # Test for PostInterval*Drug interaction
anova(basic_hlm, Minute_Sleep_interaction_hlm) # Test for PostInterval*Minute interaction

# Summary for selected model after interactions were tested
summary(Drug_Sleep_interaction_hlm)

# Remove factor Minute
rm_minute_hlm <-
  lmer(dRear ~ PostInterval * Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Remove factor PostInterval
rm_interval_hlm <-
  lmer(dRear ~ Minute + Drug + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Remove factor Drug
rm_drug_hlm <-
  lmer(dRear ~ Minute + PostInterval + (1 | animal),
       data = dRear_long,
       REML = FALSE)

# Test nested models for single factors
anova(rm_minute_hlm, Drug_Sleep_interaction_hlm) # Test for Factor Minute
anova(rm_interval_hlm, Drug_Sleep_interaction_hlm) # Test for Factor PostInterval
anova(rm_drug_hlm, Drug_Sleep_interaction_hlm) # Test for Factor Drug

# Check assumptions of selected model
plot_model(Drug_Sleep_interaction_hlm, type = 'diag')
tab_model(Drug_Sleep_interaction_hlm, show.intercept = FALSE)

# Post-hoc tests for Post-Interval*Drug Interaction ------------------------
# Post hoc analyses include two-way mixed models with animals as random
# effect and factors PostInterval/Drug and Minute for four different groups:
# 1. Sleep group
# 2. Wake groups
# 3. Muscimol infusions
# 4. Saline infusions


# Post hoc Sleep condition ------------------------------------------------
df_sleep <- subset(dRear_long, dRear_long$PostInterval == "Sleep")

basic_sleep_hlm <-
  lmer(dRear ~ Minute + Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

rm_drug_sleep_hlm <-
  lmer(dRear ~ Minute + (1 | animal),
       data = df_sleep,
       REML = FALSE)

rm_min_sleep_hlm <-
  lmer(dRear ~ Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

complete_sleep_hlm <-
  lmer(dRear ~ Minute * Drug + (1 | animal),
       data = df_sleep,
       REML = FALSE)

# Post-hoc sleep group
anova(rm_drug_sleep_hlm, basic_sleep_hlm) # test factor drug
anova(rm_min_sleep_hlm, basic_sleep_hlm) # test factor minute
anova(basic_sleep_hlm, complete_sleep_hlm) # test Drug*Minute interaction

tab_model(complete_sleep_hlm, show.intercept = FALSE)

# Post hoc Wake condition -------------------------------------------------
df_wake <- subset(dRear_long, dRear_long$PostInterval == "Wake")

basic_wake_hlm <-
  lmer(dRear ~ Minute + Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

rm_drug_wake_hlm <-
  lmer(dRear ~ Minute + (1 | animal),
       data = df_wake,
       REML = FALSE)

rm_min_wake_hlm <-
  lmer(dRear ~ Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

complete_wake_hlm <-
  lmer(dRear ~ Minute * Drug + (1 | animal),
       data = df_wake,
       REML = FALSE)

# Post-hoc wake group
anova(rm_drug_wake_hlm, basic_wake_hlm) # test factor drug
anova(rm_min_wake_hlm, basic_wake_hlm) # test factor minute
anova(basic_wake_hlm, complete_wake_hlm) # test Drug*Minute interaction

tab_model(complete_wake_hlm, show.intercept = FALSE)


# Post hoc Muscimol condition ---------------------------------------------
df_mus <- subset(dRear_long, dRear_long$Drug == "Muscimol")

basic_mus_hlm <-
  lmer(dRear ~ Minute + PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

rm_interval_mus_hlm <-
  lmer(dRear ~ Minute + (1 | animal),
       data = df_mus,
       REML = FALSE)

rm_min_mus_hlm <-
  lmer(dRear ~ PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

complete_mus_hlm <-
  lmer(dRear ~ Minute * PostInterval + (1 | animal),
       data = df_mus,
       REML = FALSE)

# Post-hoc mus group
anova(rm_interval_mus_hlm, basic_mus_hlm) # test factor PostInterval
anova(rm_min_mus_hlm, basic_mus_hlm) # test factor minute
anova(basic_mus_hlm, complete_mus_hlm) # test PostInterval*Minute interaction

tab_model(complete_mus_hlm, show.intercept = FALSE)

# Post hoc Saline condition -----------------------------------------------
df_sal <- subset(dRear_long, dRear_long$Drug == "Saline")

basic_sal_hlm <-
  lmer(dRear ~ Minute + PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

rm_interval_sal_hlm <-
  lmer(dRear ~ Minute + (1 | animal),
       data = df_sal,
       REML = FALSE)

rm_min_sal_hlm <-
  lmer(dRear ~ PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

complete_sal_hlm <-
  lmer(dRear ~ Minute * PostInterval + (1 | animal),
       data = df_sal,
       REML = FALSE)

# Post-hoc sal group
anova(rm_interval_sal_hlm, basic_sal_hlm) # test factor PostInterval
anova(rm_min_sal_hlm, basic_sal_hlm) # test factor minute
anova(basic_sal_hlm, complete_sal_hlm) # test PostInterval*Minute interaction

tab_model(complete_sal_hlm, show.intercept = FALSE)

# Multiple pairwise comparisons -------------------------------------------

## DR 1
dRear1_Sleep_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_1' &
    dRear_long$PostInterval == 'Sleep' &
    dRear_long$Drug == "Muscimol"
), ]

dRear3_Sleep_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3' &
    dRear_long$PostInterval == 'Sleep' &
    dRear_long$Drug == "Muscimol"
), ]

dRear1_Sleep_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_1'
  &
    dRear_long$PostInterval == 'Sleep' & dRear_long$Drug == "Saline"
), ]

dRear3_Sleep_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Sleep' & dRear_long$Drug == "Saline"
), ]

dRear1_Wake_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_1'
  &
    dRear_long$PostInterval == 'Wake' &
    dRear_long$Drug == "Muscimol"
), ]

dRear3_Wake_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Wake' &
    dRear_long$Drug == "Muscimol"
), ]

dRear1_Wake_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_1'
  &
    dRear_long$PostInterval == 'Wake' & dRear_long$Drug == "Saline"
), ]

dRear3_Wake_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Wake' & dRear_long$Drug == "Saline"
), ]


dRear3_Sleep_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3' &
    dRear_long$PostInterval == 'Sleep' &
    dRear_long$Drug == "Muscimol"
), ]

dRear3_Sleep_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Sleep' & dRear_long$Drug == "Saline"
), ]

dRear3_Wake_Mus <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Wake' &
    dRear_long$Drug == "Muscimol"
), ]

dRear3_Wake_Sal <- dRear_long[which(
  dRear_long$Minute == 'dRear_mean_dur_3'
  &
    dRear_long$PostInterval == 'Wake' & dRear_long$Drug == "Saline"
), ]


t.test(
  dRear1_Sleep_Mus$dRear,
  dRear1_Sleep_Sal$dRear,
  alternative = "two.sided",
  paired = TRUE
)

t.test(
  dRear3_Sleep_Mus$dRear,
  dRear3_Sleep_Sal$dRear,
  alternative = "two.sided",
  paired = TRUE
)

t.test(
  dRear1_Sleep_Sal$dRear,
  dRear1_Wake_Sal$dRear,
  alternative = "two.sided",
  paired = FALSE
)

t.test(
  dRear1_Sleep_Mus$dRear,
  dRear1_Wake_Mus$dRear,
  alternative = "two.sided",
  paired = FALSE
)

t.test(
  dRear1_Sleep_Mus$dRear,
  dRear1_Wake_Sal$dRear,
  alternative = "two.sided",
  paired = FALSE
)

var.test(dRear1_Wake_Mus$dRear, dRear1_Wake_Sal$dRear)

t.test(
  dRear1_Wake_Mus$dRear,
  dRear1_Wake_Sal$dRear,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)

var.test(dRear3_Wake_Mus$dRear, dRear3_Wake_Sal$dRear)

t.test(
  dRear3_Wake_Mus$dRear,
  dRear3_Wake_Sal$dRear,
  alternative = "two.sided",
  paired = FALSE,
  var.equal = TRUE
)

t.test(dRear1_Sleep_Mus$dRear, mu = 0)
describe(dRear1_Sleep_Mus$dRear)
t.test(dRear1_Sleep_Sal$dRear, mu = 0)
describe(dRear1_Sleep_Sal$dRear)

t.test(dRear3_Sleep_Mus$dRear, mu = 0)
t.test(dRear3_Sleep_Sal$dRear, mu = 0)
describe(dRear3_Sleep_Sal$dRear)

t.test(dRear1_Wake_Mus$dRear, mu = 0)
describe(dRear1_Wake_Mus$dRear)
t.test(dRear1_Wake_Sal$dRear, mu = 0)
describe(dRear1_Wake_Sal$dRear)

t.test(dRear3_Wake_Mus$dRear, mu = 0)
describe(dRear3_Wake_Mus$dRear)
t.test(dRear3_Wake_Sal$dRear, mu = 0)
describe(dRear3_Wake_Sal$dRear)

rear_DR_cor_1 <- merge(
  dRear1_Sleep_Sal,
  DR1_Sleep_Sal,
  by = c("animal", "animal")
)

ggscatter(
  rear_DR_cor_1,
  x = "DR",
  y = "dRear",
  title = "Correlation context and object exploration (1st min)",
  color = "black",
  shape = 21,
  size = 3,
  # Points color, shape and size
  add = "reg.line",
  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"),
  # Customize reg. line
  conf.int = TRUE,
  # Add confidence interval
  cor.coef = TRUE,
  # Add correlation coefficient. see ?stat_cor
  cor.coeff.args = list(
    method = "spearman",
    label.y = 1,
    label.sep = "\n"
  )
)

rear_DR_cor_3 <- merge(
  dRear3_Sleep_Sal,
  DR3_Sleep_Sal,
  by = c("animal", "animal")
)

ggscatter(
  rear_DR_cor_3,
  x = "dRear",
  y = "DR",
  title = "Correlation context and object exploration (1st min)",
  color = "black",
  shape = 21,
  size = 3,
  # Points color, shape and size
  add = "reg.line",
  # Add regressin line
  add.params = list(color = "blue", fill = "lightgray"),
  # Customize reg. line
  conf.int = TRUE,
  # Add confidence interval
  cor.coef = TRUE,
  # Add correlation coefficient. see ?stat_cor
  cor.coeff.args = list(
    method = "spearman",
    label.y = 1,
    label.sep = "\n"
  )
)



# 7 - Test Control parameters ---------------------------------------------
# Distance traveled
test_dist <-
  Test_wide[c("animal", "PostInterval", "Drug", "test_dist_1", "test_dist_3")]

test_dist <-
  pivot_longer(test_dist,
               cols = 4:5 ,
               names_to = "Minute",
               values_to = "dist")

test_dist$dist <- as.numeric(as.character(test_dist$dist))
test_dist$Minute <- as.factor(test_dist$Minute)
test_dist$Condition <- paste(test_dist$PostInterval, test_dist$Drug)

test_dist_summary = describeBy(
  test_dist$dist,
  list(test_dist$Condition, test_dist$Minute),
  mat = TRUE,
  digits = 2
)

test_dist_plot <-
  ggplot(test_dist_summary,
         aes(x = group2, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Distance travelled [m]",
                     breaks = seq(0, 30, 5),
                     limits = c(0, 30)) +
  scale_x_discrete(
    name = "Minute",
    labels = c("1", "3"),
    limits = c("test_dist_1",
               "test_dist_3")
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Test | Distance travelled (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

test_dist_plot

Test_Dis_sal <- subset(test_dist,test_dist$Drug == "Saline")

Test_Dis_basic_sal_hlm <-
  lmer(dist ~ Minute + PostInterval + (1 | animal),
       data = Test_Dis_sal,
       REML = FALSE)

Test_Dis_rm_interval_sal_hlm <-
  lmer(dist ~ Minute + (1 | animal),
       data = Test_Dis_sal,
       REML = FALSE)

Test_Dis_rm_min_sal_hlm <-
  lmer(dist ~ PostInterval + (1 | animal),
       data = Test_Dis_sal,
       REML = FALSE)

Test_Dis_complete_sal_hlm <-
  lmer(dist ~ Minute * PostInterval + (1 | animal),
       data = Test_Dis_sal,
       REML = FALSE)

# Post-hoc sal group
anova(Test_Dis_rm_interval_sal_hlm, Test_Dis_basic_sal_hlm) # test factor PostInterval
anova(Test_Dis_rm_min_sal_hlm, Test_Dis_basic_sal_hlm) # test factor minute
anova(Test_Dis_basic_sal_hlm, Test_Dis_complete_sal_hlm) # test PostInterval*Minute interaction

Test_dist3_SAL <- test_dist[which(
  test_dist$Minute == 'test_dist_3' & test_dist$Drug == "Saline"
), ]

Test_dist3_MUS <- test_dist[which(
  test_dist$Minute == 'test_dist_3' & test_dist$Drug == "Muscimol"
), ]


t.test(
  Test_dist3_SAL$dist,
  Test_dist3_MUS$dist,
  alternative = "two.sided",
  paired = FALSE
)


# Further 
basic_hlm <-
  lmer(dist ~ Minute + PostInterval + Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

Drug_Sleep_interaction_hlm <-
  lmer(dist ~ Minute + PostInterval * Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Add PostInterval*Minute
Minute_Sleep_interaction_hlm <-
  lmer(dist ~ Minute * PostInterval + Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Three way mixed model with all possible interactions
complete_hlm <-
  lmer(dist ~ Minute * PostInterval * Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Test nested interactions
anova(basic_hlm, Drug_Sleep_interaction_hlm) # Test for PostInterval*Drug interaction
anova(basic_hlm, Minute_Sleep_interaction_hlm) # Test for PostInterval*Minute interaction

# Summary for selected model after interactions were tested
summary(basic_hlm)

# Remove factor Minute
rm_minute_hlm <-
  lmer(dist ~ PostInterval + Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Remove factor PostInterval
rm_interval_hlm <-
  lmer(dist ~ Minute + Drug + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Remove factor Drug
rm_drug_hlm <-
  lmer(dist ~ Minute + PostInterval + (1 | animal),
       data = test_dist,
       REML = FALSE)

# Test nested models for single factors
anova(rm_minute_hlm, basic_hlm) # Test for Factor Minute
anova(rm_interval_hlm, basic_hlm) # Test for Factor PostInterval
anova(rm_drug_hlm, basic_hlm) # Test for Factor Drug

# Total exploration time
test_total_expl <-
  Test_wide[c(
    "animal",
    "PostInterval",
    "Drug",
    "test_total_expl_time_1",
    "test_total_expl_time_3"
  )]

test_total_expl <-
  pivot_longer(
    test_total_expl,
    cols = 4:5 ,
    names_to = "Minute",
    values_to = "expl_time"
  )

test_total_expl$expl_time <-
  as.numeric(as.character(test_total_expl$expl_time))
test_total_expl$Minute <- as.factor(test_total_expl$Minute)
test_total_expl$Condition <-
  paste(test_total_expl$PostInterval, test_total_expl$Drug)

test_total_expl_summary = describeBy(
  test_total_expl$expl_time,
  list(test_total_expl$Condition, test_total_expl$Minute),
  mat = TRUE,
  digits = 2
)

test_expl_time_plot <-
  ggplot(test_total_expl_summary,
         aes(x = group2, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Total exploration time [s]",
                     breaks = seq(0, 7, 1),
                     limits = c(0, 7)) +
  scale_x_discrete(
    name = "Minute",
    labels = c("1", "3"),
    limits = c("test_total_expl_time_1",
               "test_total_expl_time_3")
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Test | Total exploration time (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

test_expl_time_plot

Test_expl_sal <- subset(test_total_expl,test_total_expl$Drug == "Saline")

Test_expl_basic_sal_hlm <-
  lmer(expl_time ~ Minute + PostInterval + (1 | animal),
       data = Test_expl_sal,
       REML = FALSE)

Test_expl_rm_interval_sal_hlm <-
  lmer(expl_time ~ Minute + (1 | animal),
       data = Test_expl_sal,
       REML = FALSE)

Test_expl_rm_min_sal_hlm <-
  lmer(expl_time ~ PostInterval + (1 | animal),
       data = Test_expl_sal,
       REML = FALSE)

Test_expl_complete_sal_hlm <-
  lmer(expl_time ~ Minute * PostInterval + (1 | animal),
       data = Test_expl_sal,
       REML = FALSE)

# Post-hoc sal group
anova(Test_expl_rm_interval_sal_hlm, Test_expl_basic_sal_hlm) # test factor PostInterval
anova(Test_expl_rm_min_sal_hlm, Test_expl_basic_sal_hlm) # test factor minute
anova(Test_expl_basic_sal_hlm, Test_expl_complete_sal_hlm) # test PostInterval*Minute interaction

Test_expl3_SAL <- test_total_expl[which(
  test_total_expl$Minute == 'test_total_expl_time_3' & test_total_expl$Drug == "Saline"
), ]

Test_expl3_MUS <- test_total_expl[which(
  test_total_expl$Minute == 'test_total_expl_time_3' & test_total_expl$Drug == "Muscimol"
), ]


t.test(
  Test_expl3_SAL$expl_time,
  Test_expl3_MUS$expl_time,
  alternative = "two.sided",
  paired = FALSE
)


# Number of rearing events
test_rear_count <-
  Test_wide[c("animal", "PostInterval", "Drug", "test_rear_count_1", "test_rear_count_3")]

test_rear_count <-
  pivot_longer(test_rear_count,
               cols = 4:5 ,
               names_to = "Minute",
               values_to = "rear_count")

test_rear_count$rear_count <- as.numeric(as.character(test_rear_count$rear_count))
test_rear_count$Minute <- as.factor(test_rear_count$Minute)
test_rear_count$Condition <- paste(test_rear_count$PostInterval, test_rear_count$Drug)

test_rear_count_summary = describeBy(
  test_rear_count$rear_count,
  list(test_rear_count$Condition, test_rear_count$Minute),
  mat = TRUE,
  digits = 2
)

test_rear_count_plot <-
  ggplot(test_rear_count_summary,
         aes(x = group2, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Number of rearing events",
                     breaks = seq(0, 50, 5),
                     limits = c(0, 50)) +
  scale_x_discrete(
    name = "Minute",
    labels = c("1", "3"),
    limits = c("test_rear_count_1",
               "test_rear_count_3")
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Test | Number of rearing events (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

test_rear_count_plot

test_rear_count_sal <- subset(test_rear_count,test_rear_count$Drug == "Saline")

test_rear_count_basic_sal_hlm <-
  lmer(rear_count ~ Minute + PostInterval + (1 | animal),
       data = test_rear_count_sal,
       REML = FALSE)

test_rear_count_rm_interval_sal_hlm <-
  lmer(rear_count ~ Minute + (1 | animal),
       data = test_rear_count_sal,
       REML = FALSE)

test_rear_count_rm_min_sal_hlm <-
  lmer(rear_count ~ PostInterval + (1 | animal),
       data = test_rear_count_sal,
       REML = FALSE)

test_rear_count_complete_sal_hlm <-
  lmer(rear_count ~ Minute * PostInterval + (1 | animal),
       data = test_rear_count_sal,
       REML = FALSE)

# Post-hoc sal group
anova(test_rear_count_rm_interval_sal_hlm, test_rear_count_basic_sal_hlm) # test factor PostInterval
anova(test_rear_count_rm_min_sal_hlm, test_rear_count_basic_sal_hlm) # test factor minute
anova(test_rear_count_basic_sal_hlm, test_rear_count_complete_sal_hlm) # test PostInterval*Minute interaction

Test_rear_count3_SAL <- test_rear_count[which(
  test_rear_count$Minute == 'test_rear_count_3' & test_rear_count$Drug == "Saline"
), ]

Test_rear_count3_MUS <- test_rear_count[which(
  test_rear_count$Minute == 'test_rear_count_3' & test_rear_count$Drug == "Muscimol"
), ]


t.test(
  Test_rear_count3_SAL$rear_count,
  Test_rear_count3_MUS$rear_count,
  alternative = "two.sided",
  paired = FALSE
)


# Mean duration of rearing events 
test_rear_mean_dur <-
  Test_wide[c("animal", "PostInterval", "Drug", "test_rear_mean_dur_1", "test_rear_mean_dur_3")]

test_rear_mean_dur <-
  pivot_longer(test_rear_mean_dur,
               cols = 4:5 ,
               names_to = "Minute",
               values_to = "rear_mean_dur")

test_rear_mean_dur$rear_mean_dur <- as.numeric(as.character(test_rear_mean_dur$rear_mean_dur))
test_rear_mean_dur$Minute <- as.factor(test_rear_mean_dur$Minute)
test_rear_mean_dur$Condition <- paste(test_rear_mean_dur$PostInterval, test_rear_mean_dur$Drug)

test_rear_mean_dur_summary = describeBy(
  test_rear_mean_dur$rear_mean_dur,
  list(test_rear_mean_dur$Condition, test_rear_mean_dur$Minute),
  mat = TRUE,
  digits = 2
)

test_rear_mean_dur_plot <-
  ggplot(test_rear_mean_dur_summary,
         aes(x = group2, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Mean duration of rearing events [s]",
                     breaks = seq(0, 1.5, 0.1),
                     limits = c(0, 1.5)) +
  scale_x_discrete(
    name = "Minute",
    labels = c("1", "3"),
    limits = c("test_rear_mean_dur_1",
               "test_rear_mean_dur_3")
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Test | Mean duration of rearing events (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

test_rear_mean_dur_plot

Test_rear_mean_dur3_SAL <- test_rear_mean_dur[which(
  test_rear_mean_dur$Minute == 'test_rear_mean_dur_3' & test_rear_mean_dur$Drug == "Saline"
), ]

Test_rear_mean_dur3_MUS <- test_rear_mean_dur[which(
  test_rear_mean_dur$Minute == 'test_rear_mean_dur_3' & test_rear_mean_dur$Drug == "Muscimol"
), ]


t.test(
  Test_rear_mean_dur3_SAL$rear_mean_dur,
  Test_rear_mean_dur3_MUS$rear_mean_dur,
  alternative = "two.sided",
  paired = FALSE
)



# 8 - Sampling Control parameters -----------------------------------------
# Distance traveled
enc_dist <-
  Test_wide[c("animal", "PostInterval", "Drug", "enc_dist_10")]

enc_dist$enc_dist_10 <-
  as.numeric(as.character(enc_dist$enc_dist_10))
enc_dist$Condition <- paste(enc_dist$PostInterval, enc_dist$Drug)

enc_dist_summary = describeBy(enc_dist$enc_dist_10,
                              list(enc_dist$Condition),
                              mat = TRUE,
                              digits = 2)

enc_dist_plot <-
  ggplot(enc_dist_summary,
         aes(x = group1, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Distance travelled [m]",
                     breaks = seq(0, 70, 10),
                     limits = c(0, 70)) +
  scale_x_discrete(name = "Total distance traveled (10 minutes)") +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Sampling | Distance travelled (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)


enc_dist_plot

enc_dist_Sleep_SAL <- enc_dist[which(
  enc_dist$PostInterval == 'Sleep' & enc_dist$Drug == "Saline"
), ]

enc_dist_Wake_SAL <- enc_dist[which(
  enc_dist$PostInterval == 'Wake' & enc_dist$Drug == "Saline"
), ]


t.test(
  enc_dist_Sleep_SAL$enc_dist_10,
  enc_dist_Wake_SAL$enc_dist_10,
  alternative = "two.sided",
  paired = FALSE
)

enc_dist_SAL <- enc_dist[which(
 enc_dist$Drug == "Saline"
), ]

enc_dist_MUS <- enc_dist[which(
 enc_dist$Drug == "Muscimol"
), ]


t.test(
  enc_dist_SAL$enc_dist_10,
  enc_dist_MUS$enc_dist_10,
  alternative = "two.sided",
  paired = FALSE
)


# Total exploration time
enc_total_expl <-
  Test_wide[c("animal", "PostInterval", "Drug", "enc_total_expl_time_10")]

enc_total_expl$enc_total_expl_time_10 <-
  as.numeric(as.character(enc_total_expl$enc_total_expl_time_10))
enc_total_expl$Condition <-
  paste(enc_total_expl$PostInterval, enc_total_expl$Drug)

enc_total_expl_summary = describeBy(
  enc_total_expl$enc_total_expl_time_10,
  list(enc_total_expl$Condition),
  mat = TRUE,
  digits = 2
)

enc_expl_time_plot <-
  ggplot(enc_total_expl_summary,
         aes(x = group1, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Total exploration time [s]",
                     breaks = seq(0, 7, 1),
                     limits = c(0, 7)) +
  scale_x_discrete(name = "Total exploration time (10 minutes)") +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Sampling | Total exploration time (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

enc_expl_time_plot

enc_expl_Sleep_SAL <- enc_total_expl[which(
  enc_total_expl$PostInterval == 'Sleep' & enc_total_expl$Drug == "Saline"
), ]

enc_expl_Wake_SAL <- enc_total_expl[which(
  enc_total_expl$PostInterval == 'Wake' & enc_total_expl$Drug == "Saline"
), ]


t.test(
  enc_expl_Sleep_SAL$enc_total_expl_time_10,
  enc_expl_Wake_SAL$enc_total_expl_time_10,
  alternative = "two.sided",
  paired = FALSE
)

enc_expl_SAL <- enc_total_expl[which(
 enc_total_expl$Drug == "Saline"
), ]

enc_expl_MUS <- enc_total_expl[which(
 enc_total_expl$Drug == "Muscimol"
), ]


t.test(
  enc_expl_SAL$enc_total_expl_time_10,
  enc_expl_MUS$enc_total_expl_time_10,
  alternative = "two.sided",
  paired = FALSE
)

# Number of rearing events
enc_rear_count <-
  Test_wide[c("animal", "PostInterval", "Drug", "enc_rear_count_10")]


enc_rear_count$enc_rear_count_10 <- as.numeric(as.character(enc_rear_count$enc_rear_count_10))
enc_rear_count$Condition <- paste(enc_rear_count$PostInterval, enc_rear_count$Drug)

enc_rear_count_summary = describeBy(
  enc_rear_count$enc_rear_count_10,
  list(enc_rear_count$Condition),
  mat = TRUE,
  digits = 2
)

enc_rear_count_plot <-
  ggplot(enc_rear_count_summary,
         aes(x = group1, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Number of rearing events",
                     breaks = seq(0, 140, 20),
                     limits = c(0, 140)) +
  scale_x_discrete(
    name = "Total number of rearing events (10 minutes)"
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Sampling | Total number of rearing events (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

enc_rear_count_plot

enc_rear_count_Sleep_SAL <- enc_rear_count[which(
  enc_rear_count$PostInterval == 'Sleep' & enc_rear_count$Drug == "Saline"
), ]

enc_rear_count_Wake_SAL <- enc_rear_count[which(
  enc_rear_count$PostInterval == 'Wake' & enc_rear_count$Drug == "Saline"
), ]


t.test(
  enc_rear_count_Sleep_SAL$enc_rear_count_10,
  enc_rear_count_Wake_SAL$enc_rear_count_10,
  alternative = "two.sided",
  paired = FALSE
)

enc_rear_count_SAL <- enc_rear_count[which(
 enc_rear_count$Drug == "Saline"
), ]

enc_rear_count_MUS <- enc_rear_count[which(
enc_rear_count$Drug == "Muscimol"
), ]


t.test(
  enc_rear_count_SAL$enc_rear_count_10,
  enc_rear_count_MUS$enc_rear_count_10,
  alternative = "two.sided",
  paired = FALSE
)

# Mean duration of rearing events 
enc_rear_mean_dur <-
  Test_wide[c("animal", "PostInterval", "Drug", "enc_rear_mean_dur_10")]

enc_rear_mean_dur$enc_rear_mean_dur_10 <- as.numeric(as.character(enc_rear_mean_dur$enc_rear_mean_dur_10))
enc_rear_mean_dur$Condition <- paste(enc_rear_mean_dur$PostInterval, enc_rear_mean_dur$Drug)

enc_rear_mean_dur_summary = describeBy(
  enc_rear_mean_dur$enc_rear_mean_dur_10,
  list(enc_rear_mean_dur$Condition),
  mat = TRUE,
  digits = 2
)

enc_rear_mean_dur_plot <-
  ggplot(enc_rear_mean_dur_summary,
         aes(x = group1, y = mean, fill = group1)) +
  geom_bar(stat = 'identity',
           position = dodge,
           width = .8) +
  geom_errorbar(limits, position = dodge, width = 0.3) +
  theme_minimal() +
  scale_y_continuous(name = "Mean duration of rearing events [s]",
                     breaks = seq(0, 1.5, 0.1),
                     limits = c(0, 1.5)) +
  scale_x_discrete(
    name = "Mean duration after 10 minutes sampling"
  ) +
  scale_fill_grey(
    name = "Condition",
    labels = c(
      "Sleep Muscimol (N = 11)",
      "Sleep Saline (N = 11)",
      "Wake Muscimol (N = 8)",
      "Wake Saline (N = 10)"
    ),
    limits = c("Sleep Muscimol", "Sleep Saline", "Wake Muscimol", "Wake Saline")
  ) +
  ggtitle("Sampling | Mean duration of rearing events (cumulative)") +
  geom_hline(colour = "black",
             yintercept = 0,
             size = .1)

enc_rear_mean_dur_plot

enc_rear_mean_dur_Sleep_SAL <- enc_rear_mean_dur[which(
  enc_rear_mean_dur$PostInterval == 'Sleep' & enc_rear_mean_dur$Drug == "Saline"
), ]

enc_rear_mean_dur_Wake_SAL <- enc_rear_mean_dur[which(
  enc_rear_mean_dur$PostInterval == 'Wake' & enc_rear_mean_dur$Drug == "Saline"
), ]


t.test(
  enc_rear_mean_dur_Sleep_SAL$enc_rear_mean_dur_10,
  enc_rear_mean_dur_Wake_SAL$enc_rear_mean_dur_10,
  alternative = "two.sided",
  paired = FALSE
)

enc_rear_mean_dur_SAL <- enc_rear_mean_dur[which(
   enc_rear_mean_dur$Drug == "Saline"
), ]

enc_rear_mean_dur_MUS <- enc_rear_mean_dur[which(
  enc_rear_mean_dur$Drug == "Muscimol"
), ]


t.test(
  enc_rear_mean_dur_SAL$enc_rear_mean_dur_10,
  enc_rear_mean_dur_MUS$enc_rear_mean_dur_10,
  alternative = "two.sided",
  paired = FALSE
)

