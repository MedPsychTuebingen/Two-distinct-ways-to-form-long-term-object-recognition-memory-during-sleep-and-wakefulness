# Sawangjit et al. (2021): 1 week remote novel object recognition (NOR) memory task
# Correlation analysis for slow oscillations/spindle properties (Mean duration, Mean power, Number)
# By: Max Harkotte
# Contact: maximilian.harkotte@gmail.com
# Last update: January 2021

rm(list = ls()) # clear workspace
cat("\014") # clear console

# 0 - Load packages -------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(psych)
library(cocor)

# 1 - Source file ---------------------------------------------------------
dataPath <- "Z:/Max/03_Sleep_vs_wake_consolidation/1wk_NOR_new/Behavior/"
setwd(dataPath)

# 2 - Read in data --------------------------------------------------------
coupling_SO_left <-
  read.csv2(paste0(dataPath, "Data/Coupling_SO_EEGLeft.csv"), sep = ",")
coupling_SO_right <-
  read.csv2(paste0(dataPath, "Data/Coupling_SO_EEGRight.csv"), sep = ",")

spindles_left <-
  read.csv2(paste0(dataPath, "Data/For Anuck/Spindles_Behavior_EEGLeft.csv"),
            sep = ",")
spindles_right <-
  read.csv2(paste0(dataPath, "Data/For Anuck/Spindles_Behavior_EEGRight.csv"),
            sep = ",")

# Data prep ---------------------------------------------------------------
spindles_left <-
  spindles_left %>% select(Animal,
                           Sampling,
                           Number_Spi,
                           Densitiy_permin,
                           Mean_Duration,
                           Mean_Power)
spindles_right <-
  spindles_right %>% select(Animal,
                            Sampling,
                            Number_Spi,
                            Densitiy_permin,
                            Mean_Duration,
                            Mean_Power)

coupling_SO_left <-
  subset(coupling_SO_left, coupling_SO_left$Animal != "05-P10-S1")
coupling_SO_left$Percent_SO_Coupled <-
  round(as.numeric(as.character(coupling_SO_left$Percent_SO_Coupled)), digits = 2)
coupling_SO_left$Percent_Spi_Coupled <-
  round(as.numeric(as.character(coupling_SO_left$Percent_Spi_Coupled)), digits = 2)
coupling_SO_left$Mean_Angle <-
  round(as.numeric(as.character(coupling_SO_left$Mean_Angle)), digits = 2)
coupling_SO_left$Median_Angle <-
  round(as.numeric(as.character(coupling_SO_left$Median_Angle)), digits = 2)
coupling_SO_left$DR_1 <-
  as.numeric(as.character(coupling_SO_left$DR_1))
coupling_SO_left$DR_3 <-
  as.numeric(as.character(coupling_SO_left$DR_3))
coupling_SO_left$DR_5 <-
  as.numeric(as.character(coupling_SO_left$DR_5))
coupling_SO_left$dRear_1 <-
  as.numeric(as.character(coupling_SO_left$dRear_1))
coupling_SO_left$dRear_3 <-
  as.numeric(as.character(coupling_SO_left$dRear_3))
coupling_SO_left$dRear_5 <-
  as.numeric(as.character(coupling_SO_left$dRear_5))
coupling_SO_left <-
  merge(
    coupling_SO_left,
    spindles_left,
    by.x = c("Animal", "Sampling"),
    by.y = c("Animal", "Sampling")
  )
coupling_SO_left$Non_coupled_Spi <-
  coupling_SO_left$Number_Spi - coupling_SO_left$Number_Spi_Coupled
coupling_mus_left <-
  subset(coupling_SO_left, coupling_SO_left$Drug == "MUS")
coupling_sal_left <-
  subset(coupling_SO_left, coupling_SO_left$Drug == "SAL")


coupling_SO_right <-
  subset(coupling_SO_right, coupling_SO_right$Animal != "05-P10-S1")
coupling_SO_right$Percent_SO_Coupled <-
  round(as.numeric(as.character(coupling_SO_right$Percent_SO_Coupled)), digits = 2)
coupling_SO_right$Percent_Spi_Coupled <-
  round(as.numeric(as.character(coupling_SO_right$Percent_Spi_Coupled)), digits = 2)
coupling_SO_right$Mean_Angle <-
  round(as.numeric(as.character(coupling_SO_right$Mean_Angle)), digits = 2)
coupling_SO_right$Median_Angle <-
  round(as.numeric(as.character(coupling_SO_right$Median_Angle)), digits = 2)
coupling_SO_right$DR_1 <-
  as.numeric(as.character(coupling_SO_right$DR_1))
coupling_SO_right$DR_3 <-
  as.numeric(as.character(coupling_SO_right$DR_3))
coupling_SO_right$DR_5 <-
  as.numeric(as.character(coupling_SO_right$DR_5))
coupling_SO_right$dRear_1 <-
  as.numeric(as.character(coupling_SO_right$dRear_1))
coupling_SO_right$dRear_3 <-
  as.numeric(as.character(coupling_SO_right$dRear_3))
coupling_SO_right$dRear_5 <-
  as.numeric(as.character(coupling_SO_right$dRear_5))
coupling_SO_right <-
  merge(
    coupling_SO_right,
    spindles_right,
    by.x = c("Animal", "Sampling"),
    by.y = c("Animal", "Sampling")
  )
coupling_SO_right$Non_coupled_Spi <-
  coupling_SO_right$Number_Spi - coupling_SO_right$Number_Spi_Coupled
coupling_mus_right <-
  subset(coupling_SO_right, coupling_SO_right$Drug == "MUS")
coupling_sal_right <-
  subset(coupling_SO_right, coupling_SO_right$Drug == "SAL")

# 3 - Correlations --------------------------------------------------------
# EEG left ----------------------------------------------------------------
ggscatter(
  coupling_SO_left,
  x = "Number_Events",
  y = "DR_1",
  add = "reg.line",
  conf.int = TRUE,
  color = "Drug",
  palette = "jco",
) + stat_cor(
  aes(color = Drug),
  method = "spearman",
  label.x = 3,
  label.y = c(1.4, 1.5)
)

library("RVAideMemoire")

spearman.ci(coupling_mus_left$Number_Events, coupling_mus_left$DR_1, nrep = 1000, conf.level = 0.95)
spearman.ci(coupling_sal_left$Number_Events, coupling_sal_left$DR_1, nrep = 1000, conf.level = 0.95)

ggscatter(
  coupling_mus_left,
  x = "DR_1",
  y = "Number_Events",
  title = "Muscimol - number of coupled events  (EEG left)",
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
    label.x = -0.1,
    label.sep = "\n"
  )
)

coupling_sal_left$Densitiy_permin <-
  as.numeric(as.character(coupling_sal_left$Densitiy_permin))

ggscatter(
  coupling_sal_left,
  x = "DR_1",
  y = "Number_Events",
  title = "Saline - number of coupled events (EEG left)",
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
    label.x = -0.1,
    label.sep = "\n"
  )
)

# EEG right ---------------------------------------------------------------
ggscatter(
  coupling_SO_right,
  x = "DR_1",
  y = "Number_Events",
  title = "Overall number of coupled events  (EEG right)",
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
    label.x = -0.1,
    label.sep = "\n"
  )
)

ggscatter(
  coupling_mus_right,
  x = "DR_1",
  y = "Number_Events",
  title = "Muscimol - number of coupled events  (EEG right)",
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
    label.x = -0.1,
    label.sep = "\n"
  )
)

ggscatter(
  coupling_sal_right,
  x = "DR_1",
  y = "Number_Events",
  title = "Saline - number of coupled events (EEG right)",
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
    label.x = -0.1,
    label.sep = "\n"
  )
)

# 4 - Regression analysis -------------------------------------------------
# Two Way ANCOVA, DV: DR 1st min; IV: Drug, Number of coupled events;
# Covariate: Non-coupled spindles (number)

# Linearity assumption
ggscatter(
  coupling_SO_left,
  x = "Non_coupled_Spi",
  y = "DR_1",
  facet.by  = c("Drug"),
  short.panel.labs = FALSE
)

ggscatter(
  coupling_SO_right,
  x = "Non_coupled_Spi",
  y = "DR_1",
  facet.by  = c("Drug"),
  short.panel.labs = FALSE
)

# Homogeneity of regression slopes
coupling_SO_right %>%
  anova_test(
    DR_1 ~ Non_coupled_Spi + Drug + Number_Events +
      Drug * Number_Events + Non_coupled_Spi * Drug +
      Non_coupled_Spi * Number_Events + Non_coupled_Spi * Number_Events *
      Drug
  )

# Normality of residuals

# Homogeneity of variances

# Stepwise regressions
lm.basic <- lm(DR_1 ~ Number_Events, data = coupling_sal_left)
summary(lm.basic)

lm.adjusted <- lm(DR_1 ~ Non_coupled_Spi + Number_Events, data = coupling_sal_left)
summary(lm.adjusted)

# Standardised:
lm.adjusted <- lm(scale(DR_1) ~ scale(Non_coupled_Spi) + scale(Number_Events), data = coupling_sal_left)
summary(lm.adjusted)

# Regression analysis
res.aov_left <- coupling_SO_left %>%
  anova_test(DR_1 ~ Non_coupled_Spi + Drug * Number_Events)
get_anova_table(res.aov_left)
summary(res.aov_left)

res.aov_left_new <- coupling_SO_left %>%
  anova_test(DR_1 ~ Number_Events + Drug * Non_coupled_Spi)
get_anova_table(res.aov_left_new)

res.aov_left_drug <- coupling_SO_left %>%
  anova_test(DR_1 ~ Number_Events + Drug)
get_anova_table(res.aov_left_drug)

write.csv(
  get_anova_table(res.aov_left),
  file.path(
    dataPath,
    "Data/For Anuck/EEGLeft_Regression_with_covariate_Non_coupled_Spi.csv"
  ),
  row.names = FALSE
)

write.csv(
  get_anova_table(res.aov_left_new),
  file.path(
    dataPath,
    "Data/For Anuck/EEGLeft_Regression_with_covariate_Coupled_events.csv"
  ),
  row.names = FALSE
)

write.csv(
  get_anova_table(res.aov_left_drug),
  file.path(
    dataPath,
    "Data/For Anuck/EEGLeft_Effect_Drug_with_covariate_Coupled_events.csv"
  ),
  row.names = FALSE
)

res.aov_right <- coupling_SO_right %>%
  anova_test(DR_1 ~ Non_coupled_Spi + Drug * Number_Events)
get_anova_table(res.aov_right)

res.aov_right_new <- coupling_SO_right %>%
  anova_test(DR_1 ~ Number_Events + Drug * Non_coupled_Spi)
get_anova_table(res.aov_right_new)

res.aov_right_drug <- coupling_SO_right %>%
  anova_test(DR_1 ~ Number_Events + Drug)
get_anova_table(res.aov_right_drug)

write.csv(
  get_anova_table(res.aov_right),
  file.path(
    dataPath,
    "Data/For Anuck/EEGright_Regression_with_covariate_Non_coupled_Spi.csv"
  ),
  row.names = FALSE
)

write.csv(
  get_anova_table(res.aov_right_new),
  file.path(
    dataPath,
    "Data/For Anuck/EEGright_Regression_with_covariate_Coupled_events.csv"
  ),
  row.names = FALSE
)

write.csv(
  get_anova_table(res.aov_right_drug),
  file.path(
    dataPath,
    "Data/For Anuck/EEGright_Effect_Drug_with_covariate_Coupled_events.csv"
  ),
  row.names = FALSE
)

# Post-hoc tests
# Effect of Coupled Events for each drug condition

Left_post_hoc <- coupling_SO_left %>%
  group_by(Drug) %>%
  anova_test(DR_1 ~ Non_coupled_Spi + Number_Events)

Left_post_hoc_new <- coupling_SO_left %>%
  group_by(Drug) %>%
  anova_test(DR_1 ~ Number_Events + Non_coupled_Spi)

write.csv(
  Left_post_hoc,
  file.path(
    dataPath,
    "Data/For Anuck/EEGLeft_Post_hoc_with_covariate_Non_coupled_Spi.csv"
  ),
  row.names = FALSE
)

write.csv(
  Left_post_hoc_new,
  file.path(
    dataPath,
    "Data/For Anuck/EEGLeft_Post_hoc_with_covariate_Coupled_Events.csv"
  ),
  row.names = FALSE
)

right_post_hoc <- coupling_SO_right %>%
  group_by(Drug) %>%
  anova_test(DR_1 ~ Non_coupled_Spi + Number_Events)

right_post_hoc_new <- coupling_SO_right %>%
  group_by(Drug) %>%
  anova_test(DR_1 ~ Number_Events + Non_coupled_Spi)

write.csv(
  right_post_hoc,
  file.path(
    dataPath,
    "Data/For Anuck/EEGright_Post_hoc_with_covariate_Non_coupled_Spi.csv"
  ),
  row.names = FALSE
)

write.csv(
  right_post_hoc_new,
  file.path(
    dataPath,
    "Data/For Anuck/EEGright_Post_hoc_with_covariate_Coupled_Events.csv"
  ),
  row.names = FALSE
)

# Mean differences
t.test(coupling_mus_left$Number_Events,
       coupling_sal_left$Number_Events,
       paired = TRUE)

t.test(coupling_mus_right$Number_Events,
       coupling_sal_right$Number_Events,
       paired = TRUE)

# Differences between rearing and NOR performance correlations

# Coupled SO-spindles 
r.SOspiRear = -0.061
r.SOspiDR = 0.76
r.RearDR = 0.055
cocor.dep.groups.overlap(r.SOspiRear, r.SOspiDR, r.RearDR, n = 10, alternative = "two.sided")

# Coupled Spindles 
r.spiRear = 0.15
r.spiDR = 0.78
r.RearDR = 0.055
cocor.dep.groups.overlap(r.spiRear, r.spiDR, r.RearDR, n = 10, alternative = "two.sided")


# Comparing correlation
df_wide <- coupling_SO_left[,c("Animal", "Drug", "DR_1", "dRear_1", "Number_Events")]
df_wide <- reshape(df_wide, idvar = "Animal", timevar = "Drug", direction = "wide")


MUS_DR_Events_res <- cor.test(df_wide$DR_1.MUS, df_wide$Number_Events.MUS, 
                           method = "spearman")
SAL_DR_Events_res <- cor.test(df_wide$DR_1.SAL, df_wide$Number_Events.SAL, 
                           method = "spearman")
DR_cor_res <- cor.test(df_wide$DR_1.MUS, df_wide$DR_1.SAL, 
                       method = "spearman")
MUS_DR_SAL_Events_res <- cor.test(df_wide$DR_1.MUS, df_wide$Number_Events.SAL, 
                               method = "spearman")
MUS_Events_SAL_DR_res <- cor.test(df_wide$Number_Events.MUS, df_wide$DR_1.SAL, 
                               method = "spearman")
MUS_Events_SAL_Events_res <- cor.test(df_wide$Number_Events.MUS, df_wide$Number_Events.SAL, 
                                method = "spearman")

r.jk = unname(MUS_DR_Events_res$estimate)
r.hm = unname(SAL_DR_Events_res$estimate)
r.jh = unname(DR_cor_res$estimate)
r.jm = unname(MUS_DR_SAL_Events_res$estimate)
r.kh = unname(MUS_Events_SAL_DR_res$estimate)
r.km = unname(MUS_Events_SAL_Events_res$estimate)


cocor.dep.groups.nonoverlap(r.jk, r.hm, r.jh, r.jm, r.kh, r.km, n = 10,
                            alternative = "two.sided", test = "all", alpha = 0.05,
                            conf.level = 0.95, null.value = 0, data.name = NULL,
                            var.labels = c("DR_Mus", "Events_Mus", "DR_Sal", "Events_Sal"), return.htest = FALSE)

# Coupled events and Change in Rearing

MUS_dRear_Events_res <- cor.test(df_wide$dRear_1.MUS, df_wide$Number_Events.MUS, 
                              method = "spearman")
SAL_dRear_Events_res <- cor.test(df_wide$dRear_1.SAL, df_wide$Number_Events.SAL, 
                              method = "spearman")
dRear_cor_res <- cor.test(df_wide$dRear_1.MUS, df_wide$dRear_1.SAL, 
                          method = "spearman")
MUS_dRear_SAL_Events_res <- cor.test(df_wide$dRear_1.MUS, df_wide$Number_Events.SAL, 
                                  method = "spearman")
MUS_Events_SAL_dRear_res <- cor.test(df_wide$Number_Events.MUS, df_wide$dRear_1.SAL, 
                                  method = "spearman")
MUS_Events_SAL_Events_res <- cor.test(df_wide$Number_Events.MUS, df_wide$Number_Events.SAL, 
                                method = "spearman")

r.jk = unname(MUS_dRear_Events_res$estimate)
r.hm = unname(SAL_dRear_Events_res$estimate)
r.jh = unname(dRear_cor_res$estimate)
r.jm = unname(MUS_dRear_SAL_Events_res$estimate)
r.kh = unname(MUS_Events_SAL_dRear_res$estimate)
r.km = unname(MUS_Events_SAL_Events_res$estimate)


cocor.dep.groups.nonoverlap(r.jk, r.hm, r.jh, r.jm, r.kh, r.km, n = 10,
                            alternative = "two.si ed", test = "all", alpha = 0.05,
                            conf.level = 0.95, null.value = 0, data.name = NULL,
                            var.labels = c("dRear_Mus", "Events_Mus", "dRear_Sal", "Events_Sal"), return.htest = FALSE)




