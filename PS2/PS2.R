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
library(estimatr)

#Set working directory
setwd(here("/PS2"))

#Import datasets
twomile <- data.frame(read_dta("2miledata_2020.dta"))
sitecov <- data.frame(read_dta("sitecovariates_2020.dta"))
allcov <- data.frame(read_dta("allcovariates_2020.dta"))
allsites <- data.frame(read_dta("allsites_2020.dta"))

###Question 1 
##part a
reg1a1 <- lm(lnmdvalhs0 ~ npl2000 + lnmeanhs8, data= allsites)
reg1a1se <- sqrt(diag(vcovHC(reg1a1, type = "HC1")))

allsites_hous_char <- allsites %>% select(fips, lnmdvalhs0, npl2000, lnmeanhs8, firestoveheat80, noaircond80, nofullkitchen80, 
                                 zerofullbath80,bedrms1_80occ,bedrms2_80occ,
                                 bedrms3_80occ, bedrms4_80occ, bedrms5_80occ,
                                 blt2_5yrs80occ,blt6_10yrs80occ,blt10_20yrs80occ,
                                 blt20_30yrs80occ, blt30_40yrs80occ, blt40_yrs80occ,
                                 detach80occ, mobile80occ)

reg1a2 <- lm(lnmdvalhs0 ~ npl2000+ lnmeanhs8+ firestoveheat80+ noaircond80+ nofullkitchen80+ 
             zerofullbath80+bedrms1_80occ+bedrms2_80occ+
             bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
             blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
             blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
             detach80occ+ mobile80occ, data = allsites)
reg1a2se <- sqrt(diag(vcovHC(reg1a2, type="HC1")))

allsites_econ_demo <- allsites %>% select(lnmdvalhs0, npl2000, lnmeanhs8, firestoveheat80, noaircond80, nofullkitchen80, 
                                            zerofullbath80,bedrms1_80occ,bedrms2_80occ,
                                            bedrms3_80occ, bedrms4_80occ, bedrms5_80occ,
                                            blt2_5yrs80occ,blt6_10yrs80occ,blt10_20yrs80occ,
                                            blt20_30yrs80occ, blt30_40yrs80occ, blt40_yrs80occ,
                                            detach80occ, mobile80occ , fips, pop_den8,
                                          shrblk8,shrhsp8,child8,old8,shrfor8,
                                          ffh8, smhse8,hsdrop8,no_hs_diploma8,
                                          ba_or_better8,unemprt8,povrat8,welfare8,
                                          avhhin8, tothsun8, ownocc8, occupied80)

reg1a3 <- lm(lnmdvalhs0 ~ lnmdvalhs0+ npl2000+ lnmeanhs8+ firestoveheat80+ noaircond80+ nofullkitchen80+ 
             zerofullbath80+bedrms1_80occ+bedrms2_80occ+
             bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
             blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
             blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
             detach80occ+ mobile80occ + fips+ pop_den8+
             shrblk8+shrhsp8+child8+old8+shrfor8+
             ffh8+ smhse8+hsdrop8+no_hs_diploma8+
             ba_or_better8+unemprt8+povrat8+welfare8+
             avhhin8+ tothsun8+ ownocc8+ occupied80, data = allsites)
reg1a3se <- sqrt(diag(vcovHC(reg1a3, type = "HC1")))

allsites_state <- allsites %>% select(lnmdvalhs0, npl2000, lnmeanhs8,statefips) %>% left_join(allsites_econ_demo)
allsites_state$statefips <- as.factor(allsites_state$statefips)

reg1a4 <- lm(lnmdvalhs0 ~ npl2000+ lnmeanhs8+ firestoveheat80+ noaircond80+ nofullkitchen80+ 
               zerofullbath80+bedrms1_80occ+bedrms2_80occ+
               bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
               blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
               blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
               detach80occ+ mobile80occ + fips+ pop_den8+
               shrblk8+shrhsp8+child8+old8+shrfor8+
               ffh8+ smhse8+hsdrop8+no_hs_diploma8+
               ba_or_better8+unemprt8+povrat8+welfare8+
               avhhin8+ tothsun8+ ownocc8+ occupied80 +
               as.factor(statefips), data = allsites)
reg1a4se <- sqrt(diag(vcovHC(reg1a4, type = "HC1")))

stargazer(reg1a1, reg1a2, reg1a3, reg1a4,
          se = list(reg1a1se,reg1a2se,reg1a3se, reg1a4se),
          type = 'text', style = "aer", 
          keep = c("npl2000","lnmeanhs8","Constant"),
          omit.stat = c("F"),
          dep.var.labels = c("ln median housing price 2000"),
          covariate.labels = c("NPL by 2000", "ln mean housing price 1980"),
          add.lines = list(c("Housing Characteristics", "N", "Y", "Y", "Y"),
                           c("Economic and Demographic", "N", "N","Y","Y"),
                           c("State Fixed Effects", "N", "N", "N", "Y"),
                           c("")),
          out = "PS2Q1a.tex")

##part b
tbl1b1 <- tableby(npl2000 ~ lnmeanhs8 + firestoveheat80 + noaircond80 + 
                    nofullkitchen80 + zerofullbath80 + 
                    blt0_1yrs80occ + blt2_5yrs80occ + blt6_10yrs80occ + blt10_20yrs80occ +
                    blt20_30yrs80occ + blt30_40yrs80occ  +
                    npl2000 + pop_den8 + shrblk8 + shrhsp8 +
                    child8 + old8 + shrfor8 + ffh8 + 
                    smhse8 + hsdrop8 + no_hs_diploma8 +
                    ba_or_better8 + unemprt8 + povrat8 + welfare8 + avhhin8, 
                  data = allsites)

#write2html(tbl1b1, file = "PS2Q1b1")


#Create a new column with dummy =1 for hrs_82 >= 28.5, = 0 for hrs_82 < 28.5
sitecov <- sitecov %>% mutate(npl83 = ifelse(hrs_82 >= 28.5, 1,0))

tbl1b2 <- tableby(npl83 ~ firestoveheat80 + noaircond80 + 
                    nofullkitchen80 + zerofullbath80 + 
                    blt0_1yrs80occ + blt2_5yrs80occ + blt6_10yrs80occ + blt10_20yrs80occ +
                    blt20_30yrs80occ + blt30_40yrs80occ  +
                    npl2000 + pop_den8 + shrblk8 + shrhsp8 +
                    child8 + shrfor8 + ffh8 + 
                    smhse8 + hsdrop8 + no_hs_diploma8 +
                    ba_or_better8 + unemprt8 + povrat8 + welfare8 + avhhin8, 
                  data = sitecov)
#write2html(tbl1b2, file = "PS2Q1b2")

#Create a new column with dummy =0 if 16.5 <= hrs < 28.5, =1 if 28.5 <= hrs < 40.5. (NA everywhere else)
sitecov <- sitecov %>% mutate(hrs_range = ifelse(hrs_82 >= 16.5 & hrs_82 < 28.5, 0,
                                                 ifelse(hrs_82 < 40.5 & hrs_82 >= 28.5, 1, NA)))

tbl1b3 <- tableby(hrs_range ~ firestoveheat80 + noaircond80 + 
                    nofullkitchen80 + zerofullbath80 + 
                    blt0_1yrs80occ + blt2_5yrs80occ + blt6_10yrs80occ + blt10_20yrs80occ +
                    blt20_30yrs80occ + blt30_40yrs80occ  +
                    npl2000 + pop_den8 + shrblk8 + shrhsp8 +
                    child8 + shrfor8 + ffh8 + 
                    smhse8 + hsdrop8 + no_hs_diploma8 +
                    ba_or_better8 + unemprt8 + povrat8 + welfare8 + avhhin8, 
                  data = sitecov)
#write2html(tbl1b3, file = "PS2Q1b3")

### Question 2
## part a in write-up
## part b
#Create density plot
ggplot(data=twomile, aes(x=hrs_82)) + stat_density(alpha = 0.5) +
  geom_vline(xintercept = 28.5, linetype='dashed', color = "blue") +
  xlab("1982 HRS Score") + ylab("Density") + ggtitle("Density of 1982 HRS Scores") + 
  scale_x_continuous(breaks = c(28.5, seq(0,85,20))) + theme_minimal()

#Run local linear regressions on either side of 28.5 hrs_82



                     