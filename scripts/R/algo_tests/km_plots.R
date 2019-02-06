# install.packages(pkgs = "survival")
# BiocManager::install("survminer")

library(survival)
library(survminer)
library(tidyverse)

# Ovarian cancer data
# futime - died or were lost to follow up
# fustat - whether censored or not
data(ovarian)
glimpse(ovarian)
help(ovarian)

# rx - therapy type
# resid.ds - tumor regression
# ecog.ps - patient performance according to ecog criteria

# Dichotomise age, change data variables
ovarian$rx <- factor(ovarian$rx,
                     levels = c("1", "2"),
                     labels = c("A", "B"))

ovarian$resid.ds <- factor(ovarian$resid.ds,
                           levels = c("1", "2"),
                           labels = c("no", "yes"))

ovarian$ecog.ps <- factor(ovarian$ecog.ps,
                          levels = c("1", "2"),
                          labels = c("good", "bad"))

# Data seems bimodal? 
# [One would imagine if clear distinction, but histogram looks to be ver much in 55-60]

hist(ovarian$age)

# mutate adds new column, preserves existing variables. transmute creates new column and drops old.
ovarian <- ovarian %>% mutate(age_group = ifelse(age >= 50, "old", "young"))

ovarian$age_group <- factor(ovarian$age_group)

# Fit survival data using Kaplan-Meier method
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object

# Stratify by treatment
fit1 <- survfit(surv_object ~ rx, data = ovarian)

# Shows treatment groups and stats at various timepoints.

ggsurvplot(fit1, data = ovarian, pval = TRUE)
summary(fit1)

# Stratify bytumor regression
fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
ggsurvplot(fit2, data = ovarian, pval = TRUE)
# Survival significant relative to tumor regression - may want follow-up study

# Cox PH - Hazard ratios can utilise multitude of variables for more informant model

# Fit coxph model
fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps,
                   data = ovarian)
ggforest(fit.coxph, data = ovarian)
