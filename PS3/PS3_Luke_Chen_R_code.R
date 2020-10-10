#Setting libraries
library(here)
library(dplyr)
library(stargazer)
library(Hmisc)
library(haven)
library(sandwich)
library(lmtest)
library(xtable)
library(arsenal)
library(ggplot2)
library(AER)
library(rdd)
library(DescTools)

#Set working directory
setwd(here("/PS3"))

#Import Data
ps3.dat_og <- data.frame(read_dta("ps3_2020.dta"))

#Create a copy for future use.
ps3.dat <- data.frame(ps3.dat_og)

#### Problem 1 ####
## Part a) Fixing random values
#Recode missing values as "NA"
ps3.dat[,c(34:42,45)] <- apply(ps3.dat_og[,c(34:42,45)], MARGIN = 2,
                               function(x) {ifelse(x==8 | x ==9, NA, x)})
ps3.dat[,c(43,46,48)] <- apply(ps3.dat_og[,c(43,46,48)], MARGIN = 2,
                               function(x) {ifelse(x==99, NA, x)})                         
ps3.dat$cigar6 <- ifelse(ps3.dat_og$cigar6 == 6, NA, ps3.dat_og$cigar6)
ps3.dat$drink5 <- ifelse(ps3.dat_og$drink5 == 5, NA, ps3.dat_og$drink5)


## Part b) Address potential sampling bias by exclusion of NA observations
#Create a dataset with an indicator varible for having an NA (=1) or not (=0)
ps3.dat_na <- ps3.dat %>% mutate(has_na = ifelse(rowSums(is.na(ps3.dat))>0, 1,0))

#See what kind of reasons for the excluded observations are prevalent (count NA's)
summary(subset(ps3.dat_na, has_na == 1)[34:48])
#This suggests that there is selection into missing values

#See a summary table with ANOVA results for difference based on NA or no NA 
#(omitted from write-up but interesting)
na_test_all <- tableby(has_na ~ ., data = ps3.dat_na)
summary(na_test_all, text=TRUE)

#See and print a summary table with ANOVA results for difference in outcome variables
na_test_out <- tableby(has_na ~ dbrwt + omaps +fmaps, data = ps3.dat_na)
summary(na_test_out, text = TRUE)
#write2html(na_test_out, file = "PS3Q1b.html")

#Remove all rows that have any NA values from the original dataset
ps3.dat <- data.frame(ps3.dat[complete.cases(ps3.dat), ]) 
labels(ps3.dat) <- labels(ps3.dat_og)

## Part c) Create a basic summary table of the final analysis dataset
stargazer(ps3.dat, type = "text", style = "qje",
          font.size = "small",
          omit.summary.stat = c("p25", "p75", "max","min","n"),
          notes = "N = 92,828 for all variables",
          covariate.labels = paste(unlist(labels(ps3.dat)),sep = ""),
          out = "PS3Q1c.tex")

#### Problem 2 ####
##Part a : Compute mean difference in outcome variables by smoking status
dbrwt_means <- ps3.dat %>% group_by(tobacco) %>% summarise(mean(dbrwt), sd(dbrwt))
omaps_means <- ps3.dat %>% group_by(tobacco) %>% summarise(mean(omaps), sd(omaps))
fmaps_means <- ps3.dat %>% group_by(tobacco) %>% summarise(mean(fmaps), sd(fmaps))
all_means <- left_join(dbrwt_means,omaps_means, by = "tobacco")
all_means <- left_join(all_means, fmaps_means, by = "tobacco")
all_means[3,c(2,4,6)] <- all_means[1,c(2,4,6)] - all_means[2, c(2,4,6)]
all_means[3,c(3,5,7)] <- c(t.test(dbrwt ~ tobacco, data = ps3.dat)$stderr,
                           t.test(omaps ~ tobacco, data = ps3.dat)$stderr,
                           t.test(fmaps ~ tobacco, data = ps3.dat)$stderr)
all_means[,1] <- NULL
row.names(all_means) <- c("Non Smoker","Smoker","Mean Difference")
all_means_tbl <- xtable(all_means, caption = c("Mean Difference defined as: $E[Outcome|Non-Smoker] - E[Outcome|Smoker]$"), digits = 5)
#print.xtable(all_means_tbl, type = "latex", file = "PS3Q2a.tex")

##Part b : Create a summary table comparing means by smoking status
vbls_by_tobacco <- tableby(tobacco ~ ., data = ps3.dat)
summary(vbls_by_tobacco, text=T)
#write2html(vbls_by_tobacco, file = "PS3Q2b.html")


## Part c: In the write-up

## Part d: Linear regression of birth weight on smoking status and observables
# First, clean up the covariates of interest a little more
ps3.dat$ormoth <- as.factor(ps3.dat$ormoth)
ps3.dat$orfath <- as.factor(ps3.dat$orfath)
ps3.dat$mrace3 <- as.factor(ps3.dat$mrace3)
ps3.dat$stresfip <- as.factor(ps3.dat$stresfip)
ps3.dat$tobacco <- relevel(as.factor(ps3.dat$tobacco),"2")

# Run the regression
naive_smoke_lm <- lm(dbrwt ~ tobacco, data = ps3.dat)
smoke_lm <- lm(data = ps3.dat, dbrwt ~ tobacco + stresfip +dmage + ormoth +
                 mrace3 + dmeduc + dmar + nlbnl + dlivord + 
                 monpre + nprevist + disllb + dfage + orfath +
                 dfeduc + dplural + anemia + cardiac + lung +
                 diabetes + herpes + chyper + pre4000 +
                 preterm + drink)
# There is no clear cluster structure so let us just use Heterskedastic-Robust SE
naive_smoke_lm.se <- sqrt(diag(vcovHC(naive_smoke_lm, type = "HC1")))
smoke_lm.se <- sqrt(diag(vcovHC(smoke_lm, type = "HC1")))

# Create results table
smoke_lm_results <- stargazer(naive_smoke_lm,smoke_lm, 
                              se = list(naive_smoke_lm.se,smoke_lm.se), 
                              type = "text",
                              covariate.labels = "Smoked While Pregnant",
                              dep.var.labels = "Birth Weight",
                              keep = c("tobacco", "Constant"),
                              add.lines = list(c("Observables", "No", "Yes"),
                                               c("", "", "")),
                              style = "qje",
                              omit.stat = c("ser","f"),
                              out = "PS3Q2d.tex")

#### Problem 3 ####
## Part a ##
#Create propensity score through logit model
pscore <- glm(tobacco ~ stresfip + dmage + ormoth +
                mrace3 + dmeduc + dmar + nlbnl + dlivord + 
                monpre + nprevist + disllb + dfage + orfath +
                dfeduc + dplural + anemia + cardiac + lung +
                diabetes + herpes + chyper + pre4000 +
                preterm + drink,
              data = ps3.dat,
              family = "binomial")

pscore2 <- glm(tobacco ~ stresfip + dmage + ormoth +
                 mrace3 + dmeduc + dmar + 
                 monpre + nprevist + disllb + dfage + orfath +
                 dfeduc + dplural +
                 diabetes + pre4000 +
                 preterm + drink,
               data = ps3.dat,
               family = "binomial")

pscore.se <- sqrt(diag(vcovHC(pscore, type = "HC1")))
pscore2.se <- sqrt(diag(vcovHC(pscore2, type = "HC1")))

pscore_pr2 <- signif(PseudoR2(pscore),digits = 3)
pscore2_pr2 <- signif(PseudoR2(pscore), digits = 3)

mean_pscore <- signif(mean(fitted.values(pscore)), digits = 3)
mean_pscore2 <- signif(mean(fitted.values(pscore2)), digits =3)

mean_sq_diff <- signif(mean(fitted.values(pscore) - fitted.values(pscore2))^2, digits =3)

stargazer(pscore, pscore2, type = "text", style = "qje",
          se = list(pscore.se, pscore2.se),
          add.lines = list(c("Mean Propensity", mean_pscore, mean_pscore2),
                           c("Mean Sq Difference", "--", mean_sq_diff),
                           c("","",""),
                           c("Model Pseudo R2", pscore_pr2, pscore2_pr2)),
          omit.stat = c("aic", "ll"),
          column.labels = c("Unrestricted", "Restricted"),
          dep.var.labels = c("Propensity Score"),
          keep = 0, 
          out = "PS3Q3a.tex")

# Results table with covariate coefficients shown.
stargazer(pscore, pscore2, type = "text", style = "qje",
          se = list(pscore.se, pscore2.se),
          add.lines = list(c("Mean Propensity", mean_pscore, mean_pscore2),
                           c("Mean Sq Difference", "--", mean_sq_diff),
                           c("","",""),
                           c("Model Pseudo R2", pscore_pr2, pscore2_pr2)),
          omit.stat = c("aic", "ll"),
          column.labels = c("Unrestricted", "Restricted"),
          dep.var.labels = c("Propensity Score"))

## Part b ##
ps3.dat <- ps3.dat %>% mutate(pscore = fitted.values(pscore))
pscore_lm <- lm(dbrwt ~ tobacco + pscore, data = ps3.dat)
pscore_lm.se <- sqrt(diag(vcovHC(pscore_lm, type = "HC1")))

stargazer(pscore_lm, se = list(pscore_lm.se), 
          style = "qje", type = "text",
          omit.stat = c("f","ser"),
          dep.var.labels = c("Birth Weight"),
          covariate.labels = c("Smoked While Pregnant", "Propensity"),
          out = "PS3Q3b.tex")

## Part c ##
# Create weights as new column variable
ps3.dat <-ps3.dat %>% mutate(lm_weight = ifelse(tobacco == 1,
                                         sqrt(1/pscore),
                                         sqrt(1/(1-pscore))))
#Create ATT weights as another column variable
ps3.dat <- ps3.dat %>% mutate(treated_weight = ifelse(tobacco==1,
                                                      1,
                                                      sqrt(pscore/(1-pscore))))

# Run the WLS regression
smoke_wls <- lm(dbrwt ~ tobacco,
                data = ps3.dat,
                weights = ps3.dat$lm_weight)


smoke_wls2 <- lm(dbrwt ~ tobacco + stresfip + dmage + ormoth +
                   mrace3 + dmeduc + dmar + 
                   monpre + nprevist + disllb + dfage + orfath +
                   dfeduc + dplural +
                   diabetes + pre4000 +
                   preterm + drink,
                 data = ps3.dat, weights = ps3.dat$lm_weight)

smoke_treat_wls <- lm(dbrwt~tobacco, 
                      data = ps3.dat,
                      weights = ps3.dat$treated_weight)

smoke_treat_wls2 <- lm(dbrwt ~ tobacco + stresfip + dmage + ormoth +
                         mrace3 + dmeduc + dmar + 
                         monpre + nprevist + disllb + dfage + orfath +
                         dfeduc + dplural +
                         diabetes + pre4000 +
                         preterm + drink, 
                       data = ps3.dat, 
                       weights = ps3.dat$treated_weight)

stargazer(smoke_wls, smoke_wls2, smoke_treat_wls, smoke_treat_wls2,
          type = "text",
          style = "qje",
          keep = c("Constant", "tobacco"),
          covariate.labels = c("Smoked While Pregnant"),
          column.labels = c("ATE", "ATET"),
          column.separate = c(2,2),
          dep.var.labels = c("Birth Weight"),
          omit.stat = c("ser","f"),
          add.lines = list(c("Observables","No", "Yes", "No", "Yes"),
                           c("","","")),
          out = "PS3Q3c.tex")

#### Problem 4 ####
## Create bins
# Equally spaced bins will be of width:
(max(ps3.dat$pscore) - min(ps3.dat$pscore))/100

# Round this value to 0.01 for computational ease.
# Create column variable for bins
bintervals <- seq(0,1,0.01)
ps3.dat <- ps3.dat %>% mutate(bin = cut(pscore, breaks = bintervals, labels = FALSE))

# Calculate mean difference in outcome for each bin/ block
block.dat <- ps3.dat %>% group_by(bin, tobacco) %>% summarise(mean(dbrwt), sd(dbrwt), n())
block.dat <- left_join(subset(block.dat,tobacco == 1),subset(block.dat,tobacco == 2), by = c("bin"))
colnames(block.dat) <- c("bin", "tobacco1", "meanbwt1", "sd_bwt1", "n1", "tobacco2", "meanbwt2", "sd_bwt2", "n2")
block.dat <- block.dat %>% mutate(bwt_dif = meanbwt1 - meanbwt2) %>% 
  mutate(bwt_dif_sd = sqrt(sd_bdmswt1^2/n1 + sd_bwt2^2/n2)) %>%
  mutate(sample_wt = ((n1 + n2)/ nrow(ps3.dat))) %>%
  mutate(treat_wt = (n1 / sum(block.dat$n1, na.rm = T)))


# Calulcate ATE by taking mean block differences weighted by proportion of total sample
ATE <- weighted.mean(x = block.dat$bwt_dif, w =  block.dat$sample_wt, na.rm = T)
ATE_sd <- sqrt(sum((block.dat$sample_wt)^2 * (block.dat$bwt_dif_sd)^2, na.rm = T))
ATE_sd

#Calculate ATET by taking mean block differences weighted by proportion of treated sample
ATT <- weighted.mean(x = block.dat$bwt_dif, w =  block.dat$treat_wt, na.rm = T)
ATT_sd <- sqrt(sum((block.dat$treat_wt)^2 * (block.dat$bwt_dif_sd)^2, na.rm = T))

block_results <- cbind(c(ATE, ATE_sd, dt(x=(ATE/ATE_sd), df = 92826)), 
                       c(ATT, ATT_sd, dt(x=(ATT/ATT_sd), df = 98826)))
row.names(block_results) <- c("Smoked While Pregnant","Coeff. Std. Error","P-Value")
block_results <- as.table(block_results)
colnames(block_results) <- c("ATE", "ATET")
#xtable(block_results, 
#       caption = "Results of Non-Parametric Propensity Score Block Estimation of Effect of Smoking on Birth Weight",
#       align = "ccc", digits = 4)
block_results

#### Problem 5 ####
#Create a new column variable for "low birth weight"
ps3.dat <- ps3.dat %>% mutate(low_bwt = ifelse(dbrwt < 2500, 1, 0))

# Re-do the block-style estimation
block.dat2 <- ps3.dat %>% group_by(bin, tobacco) %>% summarise(mean(low_bwt), sd(low_bwt), n())
block.dat2 <- left_join(subset(block.dat2,tobacco == 1),subset(block.dat2,tobacco == 2), by = c("bin"))
colnames(block.dat2) <- c("bin", "tobacco1", "low_bwt1", "sd_bwt1", "n1", "tobacco2", "low_bwt2", "sd_bwt2", "n2")
block.dat2 <- block.dat2 %>% mutate(bwt_dif = low_bwt1 - low_bwt2) %>% 
  mutate(bwt_dif_sd = sqrt(sd_bwt1^2/n1 + sd_bwt2^2/n2)) %>%
  mutate(sample_wt = ((n1 + n2)/ (sum(block.dat$n1) + sum(block.dat$n2,na.rm = T)))) %>%
  mutate(treat_wt = (n1)/sum(block.dat2$n1, na.rm = T))

# Calulcate ATE by taking mean block differences weighted by proportion of total sample
ATE2 <- weighted.mean(x = block.dat2$bwt_dif, w =  block.dat2$sample_wt, na.rm = T)
ATE_sd2 <- sqrt(sum(block.dat2$sample_wt^2*block.dat2$bwt_dif_sd^2, na.rm = T))

#Calculate ATET by taking mean block differences weighted by proportion of treated sample
ATT2 <- weighted.mean(x = block.dat2$bwt_dif, w =  block.dat2$treat_wt, na.rm = T)
ATT_sd2 <- sqrt(sum((block.dat2$treat_wt)^2 * (block.dat2$bwt_dif_sd)^2, na.rm = T))

block2_results <- cbind(c(ATE2, ATE_sd2, dt(x=(ATE2/ATE_sd2), df = nrow(ps3.dat))),
                        c(ATT2, ATT_sd2, dt(x=(ATT2/ATT_sd2), df = nrow(ps3.dat))))
                        
row.names(block2_results) <- c("Smoked While Pregnant","Coeff. Std. Error","P-Value")

block2_results <- as.table(block2_results)
colnames(block2_results) <- c("ATE", "ATET")
#xtable(block2_results, 
#       caption = "Results of Non-Parametric Propensity Score Block Estimation of Effect of Smoking on Likelihood of Low Birth Weight",
#       align = "ccc", digits = 4)
block2_results
