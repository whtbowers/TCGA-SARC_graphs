library(survival)
library(survminer)

# for coxph, formula is linear model with survival object as response variable.
# Survival object created using Surv(time, event)

data("lung")
head(lung)

# Status: 1=censored, 2 = death
# ph.karno: Karnofsky performance score - 0-100, bad-good

# Univariate cox regression

res.cox <- coxph(Surv(time, status) ~ sex + ph.ecog, data = lung)
res.cox

# Call:
#   coxph(formula = Surv(time, status) ~ sex, data = lung)
# 
# coef exp(coef) se(coef)      z       p
# sex -0.5310    0.5880   0.1672 -3.176 0.00149
# 
# Likelihood ratio test=10.63  on 1 df, p=0.001111
# n= 228, number of events= 165 

summary(res.cox)

# Multivariate Cox regression analysis
# Visualising estimated survival time distribution

# Application of univariate coxph to multiple covariates simultaneously

covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

# Simply takes results from multiple tests and stacks into single data frame
