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

#Set working directory
setwd(here())

#Import datasets
twomile <- data.frame(read_dta("2miledata_2020.dta"))
sitecov <- data.frame(read_dta("sitecovariates_2020.dta"))
allcov <- data.frame(read_dta("allcovariates_2020.dta"))
allsites <- data.frame(read_dta("allsites_2020.dta"))

###Question 1 
##part a
reg1a1 <- lm(lnmdvalhs0 ~ npl2000 + lnmeanhs8, data= allsites)
reg1a1se <- sqrt(diag(vcovHC(reg1a1, type = "HC1")))

allsites_hous_char <- allsites %>% select(lnmdvalhs0, npl2000, lnmeanhs8, firestoveheat80, noaircond80, nofullkitchen80, 
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
                                            detach80occ, mobile80occ, pop_den8,
                                          shrblk8,shrhsp8,child8,old8,shrfor8,
                                          ffh8, smhse8,hsdrop8,no_hs_diploma8,
                                          ba_or_better8,unemprt8,povrat8,welfare8,
                                          avhhin8, tothsun8, ownocc8, occupied80)
#Starting here, R require 12+ GB of RAM to process the models
#I 'cheated' by getting the results from Stata and plugging
#Them into the table in my write-up.

reg1a3 <- lm(lnmdvalhs0 ~., data = allsites_econ_demo)
reg1a3se <- sqrt(diag(vcovHC(reg1a3, type = "HC1")))

allsites_state <- allsites %>% select(lnmdvalhs0, npl2000, lnmeanhs8, firestoveheat80, noaircond80, nofullkitchen80, 
                                      zerofullbath80,bedrms1_80occ,bedrms2_80occ,
                                      bedrms3_80occ, bedrms4_80occ, bedrms5_80occ,
                                      blt2_5yrs80occ,blt6_10yrs80occ,blt10_20yrs80occ,
                                      blt20_30yrs80occ, blt30_40yrs80occ, blt40_yrs80occ,
                                      detach80occ, mobile80occ, pop_den8,
                                      shrblk8,shrhsp8,child8,old8,shrfor8,
                                      ffh8, smhse8,hsdrop8,no_hs_diploma8,
                                      ba_or_better8,unemprt8,povrat8,welfare8,
                                      avhhin8, tothsun8, ownocc8, occupied80, statefips)
allsites_state$statefips <- as.factor(allsites_state$statefips)

reg1a4 <- lm(lnmdvalhs0 ~ .,data = allsites_state)
reg1a4se <- sqrt(diag(vcovHC(reg1a4, type = "HC1")))

stargazer(reg1a1, reg1a2, reg1a3, reg1a4,
          se = list(reg1a1se,reg1a2se,reg1a3se, reg1a4se),
          type = 'text', style = "aer", 
          keep = c("npl2000","lnmeanhs8","Constant"),
          omit.stat = c("F","ser","adj.rsq"),
          dep.var.labels = c("ln median housing price 2000"),
          covariate.labels = c("NPL by 2000", "ln mean housing price 1980"),
          add.lines = list(c("Housing Characteristics", "N", "Y", "Y", "Y"),
                           c("Economic and Demographic", "N", "N","Y","Y"),
                           c("State Fixed Effects", "N", "N", "N", "Y"),
                           c("")),
          out = "PS2Q1a.tex")

##part b
tbl1b1 <- tableby(npl2000 ~  firestoveheat80+ noaircond80+ nofullkitchen80+ 
                  zerofullbath80+bedrms1_80occ+bedrms2_80occ+
                  bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
                  blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
                  blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
                  detach80occ+ mobile80occ + pop_den8+
                    shrblk8+shrhsp8+child8+old8+shrfor8+
                    ffh8+ smhse8+hsdrop8+no_hs_diploma8+
                    ba_or_better8+unemprt8+povrat8+welfare8+
                    avhhin8+ tothsun8+ ownocc8+ occupied80,
                  data = allsites)

#write2html(tbl1b1, file = "PS2Q1b1")
summary(tbl1b1, text = T)

#Create a new column with dummy =1 for hrs_82 >= 28.5, = 0 for hrs_82 < 28.5
sitecov <- sitecov %>% mutate(npl83 = ifelse(hrs_82 >= 28.5, 1,0))

tbl1b2 <- tableby(npl83 ~ firestoveheat80+ noaircond80+ nofullkitchen80+ 
                    zerofullbath80+bedrms1_80occ+bedrms2_80occ+
                    bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
                    blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
                    blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
                    detach80occ+ mobile80occ + pop_den8+
                    shrblk8+shrhsp8+child8+shrfor8+
                    ffh8+ smhse8+hsdrop8+no_hs_diploma8+
                    ba_or_better8+unemprt8+povrat8+welfare8+
                    avhhin8+ tothsun8+ ownocc8+ occupied80, 
                  data = sitecov)
#write2html(tbl1b2, file = "PS2Q1b2")
summary(tbl1b2, text = T)

#Create a new column with dummy =0 if 16.5 <= hrs < 28.5, =1 if 28.5 <= hrs < 40.5. (NA everywhere else)
sitecov <- sitecov %>% mutate(hrs_range = ifelse(hrs_82 >= 16.5 & hrs_82 < 28.5, 0,
                                                 ifelse(hrs_82 < 40.5 & hrs_82 >= 28.5, 1, NA)))

tbl1b3 <- tableby(hrs_range ~ firestoveheat80+ noaircond80+ nofullkitchen80+ 
                    zerofullbath80+bedrms1_80occ+bedrms2_80occ+
                    bedrms3_80occ+ bedrms4_80occ+ bedrms5_80occ+
                    blt2_5yrs80occ+blt6_10yrs80occ+blt10_20yrs80occ+
                    blt20_30yrs80occ+ blt30_40yrs80occ+ blt40_yrs80occ+
                    detach80occ+ mobile80occ + pop_den8+
                    shrblk8+shrhsp8+child8+shrfor8+
                    ffh8+ smhse8+hsdrop8+no_hs_diploma8+
                    ba_or_better8+unemprt8+povrat8+welfare8+
                    avhhin8+ tothsun8+ ownocc8+ occupied80, 
                  data = sitecov)
#write2html(tbl1b3, file = "PS2Q1b3")
summary(tbl1b3, text = T)

### Question 2
## part a in write-up
## part b
#Create density plot
ggplot(data=twomile, aes(x=hrs_82)) +
  geom_histogram(aes(y=stat(density)),binwidth = 3, alpha = 0.75, fill = "grey") + 
  geom_line(aes(x = hrs_82, y = ..density..), stat = 'density', color = "dark green")  + 
  geom_vline(xintercept = 28.5, linetype='dashed', color = "blue") +
  xlab("1982 HRS Score") + ylab("Density") + ggtitle("Density and Histogram of 1982 HRS Scores") + 
  scale_x_continuous(breaks = c(28.5, seq(0,85,20))) + theme_minimal()

#Create bins by hrs_82 with width = 3
bintervals <- seq(0,81,3)
twomile <- twomile %>% mutate(bins = cut(hrs_82, breaks = bintervals, labels = FALSE))
#Create dataset of midpoints and a threshold dummy
twomile_mids <- twomile %>% group_by(bins) %>% summarise(mean(npl2000)) %>% 
  ungroup
twomile_mids$bins <- twomile_mids$bins*3 - 1.5
twomile_mids <- twomile_mids %>% mutate(above = ifelse(bins >=28.5,1,0))
colnames(twomile_mids) <- c("hrs_82", "mean_npl_2000", "above")

#Run local linear regressions on either side of 28.5 hrs_82
reg2b <- lm(mean_npl_2000~ above*hrs_82 , data = twomile_mids)

stargazer(reg2b, se = list(sqrt(diag(vcovHC(reg2b,type="HC1")))),
          type = 'text', out = "PS2Q2b3.tex")

ggplot(data = twomile_mids, aes(x=hrs_82, y = mean_npl_2000, group = above)) + 
  geom_point(color = "dark green") + geom_smooth(method = "lm") + 
  ggtitle("Probability of NPL 2000 =1 with respect to HRS 82")
  

### Question 3
## part a
#Createinstrument dummy variable based on hrs_82
twomile <- twomile %>% mutate(above = ifelse(hrs_82 >= 28.5, 1, 0))

#Run first-stage regression with no controls
reg3a1 <- lm(npl2000 ~ above*hrs_82, data = twomile)
reg3a1.se <- sqrt(diag(vcovHC(reg3a1, type = "HC1")))

#Run first-stage regression with housing controls
reg3a2 <- lm(npl2000 ~ above*hrs_82 + firestoveheat80_nbr+ + lnmeanhs8_nbr +
               noaircond80_nbr+ nofullkitchen80_nbr+ 
             zerofullbath80_nbr+bedrms1_80occ_nbr+
               bedrms2_80occ_nbr+
             bedrms3_80occ_nbr+ bedrms4_80occ_nbr
             + bedrms5_80occ_nbr+
             blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
               blt10_20yrs80occ_nbr+
             blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
             detach80occ_nbr+ mobile80occ_nbr, data = twomile)
reg3a2.se <- sqrt(diag(vcovHC(reg3a2, type="HC1")))

reg3a3 <-  lm(npl2000 ~ above*hrs_82 + firestoveheat80_nbr+ + lnmeanhs8_nbr +
                noaircond80_nbr+ nofullkitchen80_nbr+ 
                zerofullbath80_nbr+bedrms1_80occ_nbr+
                bedrms2_80occ_nbr+
                bedrms3_80occ_nbr+ bedrms4_80occ_nbr
              + bedrms5_80occ_nbr+
                blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
                blt10_20yrs80occ_nbr+
                blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
                detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
                shrblk8_nbr+shrhsp8_nbr+child8_nbr+shrfor8_nbr+
                ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
                ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+
                avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr
                , data = twomile)

reg3a3.se <- sqrt(diag(vcovHC(reg3a3,type="HC1")))

reg3a4 <- lm(npl2000 ~ above*hrs_82 + firestoveheat80_nbr+ + lnmeanhs8_nbr +
               noaircond80_nbr+ nofullkitchen80_nbr+ 
               zerofullbath80_nbr+bedrms1_80occ_nbr+
               bedrms2_80occ_nbr+
               bedrms3_80occ_nbr+ bedrms4_80occ_nbr
             + bedrms5_80occ_nbr+
               blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
               blt10_20yrs80occ_nbr+
               blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
               detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
               shrblk8_nbr+shrhsp8_nbr+child8_nbr+ shrfor8_nbr+
               ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
               ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+
               avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips)
             , data = twomile)
reg3a4.se <- sqrt(diag(vcovHC(reg3a4,type="HC1")))

#Create a new dataset of subsetted values of hrs_82
twomile1 <- twomile %>% filter(hrs_82 >= 16.5 & hrs_82 <40.5) 

#Run fully-controlled first-stage on new subset
reg3a5 <- lm(npl2000 ~ above*hrs_82 + firestoveheat80_nbr+ lnmeanhs8_nbr +
               noaircond80_nbr+ nofullkitchen80_nbr+ 
               zerofullbath80_nbr+bedrms1_80occ_nbr+
               bedrms2_80occ_nbr+
               bedrms3_80occ_nbr+ bedrms4_80occ_nbr
             + bedrms5_80occ_nbr+
               blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
               blt10_20yrs80occ_nbr+
               blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
               detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
               shrblk8_nbr+shrhsp8_nbr+child8_nbr+shrfor8_nbr+
               ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
               ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+
               avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips)
             , data = twomile1)
reg3a5.se <- sqrt(diag(vcovHC(reg3a5, type="HC1")))

stargazer(reg3a1, reg3a2, reg3a3, reg3a4, reg3a5,
          se = list(reg3a1.se,reg3a2.se,reg3a3.se, reg3a4.se, reg3a5.se),
          type = 'text', style = "aer", 
          keep = c("above","above:hrs_82","hrs_82", "Constant"),
          omit.stat = c("F","ser","adj.rsq"),
          column.labels = c("Unrestricted Sample","Restricted Sample"),
          column.separate = c(4,1),
          dep.var.labels = c("Probability of Being on NPL by 2000"),
          add.lines = list(c("Housing Characteristics", "N", "Y", "Y", "Y","Y"),
                           c("Economic and Demographic", "N", "N","Y","Y","Y"),
                           c("State Fixed Effects", "N", "N", "N", "Y","Y"),
                           c("")),
          out = "PS2Q3a2.tex")
## part b
#Create plot of HRS score (x) and NPL2000 (y)
ggplot(data = twomile, aes(x= hrs_82, y=npl2000)) + geom_point(alpha = 0.5, color = "dark green") +
  xlab("HRS Score in 1982") + ylab("NPL Status in 2000") +
  ggtitle("Plot of NPL Status in 2000 Against HRS Score in 1982") +
  geom_vline(xintercept = 28.5, linetype='dashed', color = "dark blue") +
  scale_x_continuous(breaks = c(28.5, seq(0,85,20))) + theme_minimal()
  
ggplot(data = twomile, aes(x = hrs_82, y = fitted.values(reg3a4))) + geom_point(color = "dark green")+
  xlab("HRS Score in 1982") + ylab("Fitted Value of NPL Status in 2000") +
  ggtitle("Plot of Predicted NPL Status in 2000 Against HRS Score in 1982") +
  geom_vline(xintercept = 28.5, linetype='dashed', color = "dark blue") +
  scale_x_continuous(breaks = c(28.5, seq(0,85,20))) + theme_minimal()

ggplot(data = twomile, aes(x=hrs_82, y = lnmeanhs8_nbr)) + geom_point(color = "dark green") +
  xlab("HRS Score in 1982") + ylab("ln Mean House Price 1980") +
  ggtitle("Plot of 1980 House Prices Against HRS Score in 1982") +
  geom_vline(xintercept = 28.5, linetype='dashed', color = "dark blue") +
  scale_x_continuous(breaks = c(28.5, seq(0,85,20))) + theme_minimal()

### Question 4
# Run 2SLS using restricted sample and all controls
reg4a <- ivreg(lnmdvalhs0_nbr ~ npl2000 + hrs_82*above - above +
                 firestoveheat80_nbr+ 
                 noaircond80_nbr+ nofullkitchen80_nbr+ 
                 zerofullbath80_nbr+bedrms1_80occ_nbr+
                 bedrms2_80occ_nbr+
                 bedrms3_80occ_nbr+ bedrms4_80occ_nbr +
                 bedrms5_80occ_nbr+
                 blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
                 blt10_20yrs80occ_nbr+
                 blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
                 detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
                 shrblk8_nbr+shrhsp8_nbr+child8_nbr+shrfor8_nbr+
                 ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
                 ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+ lnmeanhs8_nbr +
                 avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips)|
                 above*hrs_82 + firestoveheat80_nbr+ 
                 noaircond80_nbr+ nofullkitchen80_nbr+ 
                 zerofullbath80_nbr+bedrms1_80occ_nbr+
                 bedrms2_80occ_nbr+
                 bedrms3_80occ_nbr+ bedrms4_80occ_nbr
               + bedrms5_80occ_nbr+
                 blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
                 blt10_20yrs80occ_nbr+
                 blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
                 detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
                 shrblk8_nbr+shrhsp8_nbr+child8_nbr+lnmeanhs8_nbr+shrfor8_nbr+
                 ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
                 ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+
                 avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips),
               data = twomile1)
reg4a.se <- sqrt(diag(vcovHC(reg4a,type="HC1")))

#Now same for unrestricted sample
reg4b <-ivreg(lnmdvalhs0_nbr ~ npl2000 + hrs_82*above - above +
                firestoveheat80_nbr+ 
                noaircond80_nbr+ nofullkitchen80_nbr+ 
                zerofullbath80_nbr+bedrms1_80occ_nbr+
                bedrms2_80occ_nbr+
                bedrms3_80occ_nbr+ bedrms4_80occ_nbr +
                bedrms5_80occ_nbr+
                blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
                blt10_20yrs80occ_nbr+
                blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
                detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
                shrblk8_nbr+shrhsp8_nbr+child8_nbr+shrfor8_nbr+
                ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
                ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+ lnmeanhs8_nbr +
                avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips)|
                above*hrs_82 + firestoveheat80_nbr+ 
                noaircond80_nbr+ nofullkitchen80_nbr+ 
                zerofullbath80_nbr+bedrms1_80occ_nbr+
                bedrms2_80occ_nbr+
                bedrms3_80occ_nbr+ bedrms4_80occ_nbr
              + bedrms5_80occ_nbr+
                blt2_5yrs80occ_nbr+blt6_10yrs80occ_nbr+
                blt10_20yrs80occ_nbr+
                blt20_30yrs80occ_nbr+ blt30_40yrs80occ_nbr+ blt40_yrs80occ_nbr+
                detach80occ_nbr+ mobile80occ_nbr +  pop_den8_nbr+
                shrblk8_nbr+shrhsp8_nbr+child8_nbr+lnmeanhs8_nbr+shrfor8_nbr+
                ffh8_nbr+ smhse8_nbr+hsdrop8_nbr+no_hs_diploma8_nbr+
                ba_or_better8_nbr+unemprt8_nbr+povrat8_nbr+welfare8_nbr+
                avhhin8_nbr+ tothsun8_nbr+ ownocc8_nbr+ occupied80_nbr + as.factor(statefips),
              data = twomile)
reg4b.se <- sqrt(diag(vcovHC(reg4b,type="HC1")))

stargazer(reg4b, reg4a,
          se = list(reg4b.se, reg4a.se), type = 'text' ,style = "aer",
          keep = c("npl2000","hrs_82:above","hrs_82","Constant"),
          omit.stat = c("F","ser","adj.rsq"),
          column.labels = c("Unrestricted Sample","Restricted Sample"),
          column.separate = c(1,1),
          dep.var.labels = c("Natural Log of Median Housing Value 2000"),
          covariate.labels = c("NPL 2000 (Fitted)"),
          add.lines = list(c("Housing Characteristics", "Y","Y"),
                           c("Economic and Demographic", "Y","Y"),
                           c("State Fixed Effects", "Y","Y"),
                           c("")),
          out = "PS2Q4a.tex")

#Random plot showing the non-trend lol
#ggplot(data = twomile1, aes(x=fitted.values(reg3a5), y=lnmdvalhs0_nbr)) + geom_point() + geom_smooth(method = 'lm')
