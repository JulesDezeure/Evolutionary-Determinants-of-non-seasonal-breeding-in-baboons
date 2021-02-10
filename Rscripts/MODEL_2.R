library(lme4)
library(blme)
library(MuMIn)
library(DHARMa)
library(car)

#for the selection of the best phase of the sine term (birth timing), of the best non-seasonal environmental variation, 
#of the best group birth synchrony metrics on offspring mortality probability, see the scripts ran for Model 1 (same procedure)



TAB=read.csv2('file_path of the excel tab Infant Mortality') 

#response variable: we only kept the infants whose survival outcome is known. 
TAB=TAB[which(TAB$Death_18months=="1"|TAB$Death_18months=="0"), ] #195 donnees
TAB$Death_18months=as.factor(TAB$Death_18months)

#we scale the quantitative fixed effects:
TAB$Relative_Rank=scale(TAB$Relative_Rank, center=TRUE, scale=TRUE)
TAB$Dif_NDVI=scale(TAB$Dif_NDVI)
TAB$Dif_Rainfall=scale(TAB$Dif_Rainfall)
TAB$Nb_Birth_Before1=scale(TAB$Nb_Birth_Before1)
TAB$Nb_Birth_Before3=scale(TAB$Nb_Birth_Before3)
TAB$Nb_Birth_Before6=scale(TAB$Nb_Birth_Before6)
TAB$Nb_Birth_After1=scale(TAB$Nb_Birth_After1)
TAB$Nb_Birth_After3=scale(TAB$Nb_Birth_After3)
TAB$Nb_Birth_After6=scale(TAB$Nb_Birth_After6)
TAB$Nb_Birth_Both1=scale(TAB$Nb_Birth_Both1)
TAB$Nb_Birth_Both3=scale(TAB$Nb_Birth_Both3)
TAB$Nb_Birth_Both6=scale(TAB$Nb_Birth_Both6)
TAB$NB_Adult_Females_V2=scale(TAB$NB_Adult_Females_V2)

#randomized birth dates drawn here (1000 times)
for (i in 1:999) {
  print(i)
  vec=c()
  for (j in 1:nrow(TAB)) {
    ff=runif(1,TAB$Birth_radian[j]-TAB$Uncertainty_Rad[j],TAB$Birth_radian[j]+TAB$Uncertainty_Rad[j])
    
    if (ff>2*pi) {ff=ff-2*pi}
    if (ff<0) {ff=ff+2*pi}
    vec[j]=ff
  }
  TAB=cbind(TAB,vec)
}



#here we run the full model, without any female rank interactions as these interactions were not significant... 
#taking into account birth date uncertainty. 
Estimate_sin=c()
Estimate_DifNDVI=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Sex=c()
Estimate_Group_Size=c()
Estimate_Compet_Soc=c()

CILOW_sin=c()
CIUP_sin=c()
CILOW_DifNDVI=c()
CIUP_DifNDVI=c()
CILOW_TroopL=c()
CIUP_TroopL=c()
CILOW_TroopM=c()
CIUP_TroopM=c()
CILOW_Rank=c()
CIUP_Rank=c()
CILOW_Parity=c()
CIUP_Parity=c()
CILOW_Sex=c()
CIUP_Sex=c()
CILOW_Group_Size=c()
CIUP_Group_Size=c()
CILOW_Compet_Social=c()
CIUP_Compet_Social=c()

Chisq_sin=c()
Chisq_DifNDVI=c()
Chisq_Troop=c()
Chisq_Rank=c()
Chisq_Parity=c()
Chisq_Sex=c()
Chisq_Group_Size=c()
Chisq_Compet_Social=c()


p_sin=c()
p_DifNDVI=c()
p_Troop=c()
p_Rank=c()
p_Parity=c()
p_Sex=c()
p_Group_Size=c()
p_Compet_Social=c()

j=ncol(TAB)-999
for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD=glmer(Death_18months ~ sin(TAB[,j]+pi/2) + Dif_NDVI + Troop + Relative_Rank + Parity + 
              Sex + NB_Adult_Females_V2 + Nb_Birth_After6 + (1|Mother), 
            glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
  SS=summary(MOD)$coefficients
  CC=confint(MOD,method='Wald')
  AA=Anova(MOD,type=3)
  
  Estimate_sin=c(Estimate_sin,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[4,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[5,1])
  Estimate_Rank=c(Estimate_Rank,SS[6,1])
  Estimate_Parity=c(Estimate_Parity,SS[7,1])
  Estimate_Sex=c(Estimate_Sex,SS[8,1])
  Estimate_Group_Size=c(Estimate_Group_Size,SS[9,1])
  Estimate_Compet_Soc=c(Estimate_Compet_Soc,SS[10,1])
  
  CILOW_sin=c(CILOW_sin,CC[3,1])
  CIUP_sin=c(CIUP_sin,CC[3,2])
  CILOW_DifNDVI=c(CILOW_DifNDVI,CC[4,1])
  CIUP_DifNDVI=c(CIUP_DifNDVI,CC[4,2])
  CILOW_TroopL=c(CILOW_TroopL,CC[5,1])
  CIUP_TroopL=c(CIUP_TroopL,CC[5,2])
  CILOW_TroopM=c(CILOW_TroopM,CC[6,1])
  CIUP_TroopM=c(CIUP_TroopM,CC[6,2])
  CILOW_Rank=c(CILOW_Rank,CC[7,1])
  CIUP_Rank=c(CIUP_Rank,CC[7,2])
  CILOW_Parity=c(CILOW_Parity,CC[8,1])
  CIUP_Parity=c(CIUP_Parity,CC[8,2])
  CILOW_Sex=c(CILOW_Sex,CC[9,1])
  CIUP_Sex=c(CIUP_Sex,CC[9,2])
  CILOW_Group_Size=c(CILOW_Group_Size,CC[10,1])
  CIUP_Group_Size=c(CIUP_Group_Size,CC[10,2])
  CILOW_Compet_Social=c(CILOW_Compet_Social,CC[11,1])
  CIUP_Compet_Social=c(CIUP_Compet_Social,CC[11,2])
  
  
  
  Chisq_sin=c(Chisq_sin,AA[2,1])
  Chisq_DifNDVI=c(Chisq_DifNDVI,AA[3,1])
  Chisq_Troop=c(Chisq_Troop,AA[4,1])
  Chisq_Rank=c(Chisq_Rank,AA[5,1])
  Chisq_Parity=c(Chisq_Parity,AA[6,1])
  Chisq_Sex=c(Chisq_Sex,AA[7,1])
  Chisq_Group_Size=c(Chisq_Group_Size,AA[8,1])
  Chisq_Compet_Social=c(Chisq_Compet_Social,AA[9,1])
  
  
  
  p_sin=c(p_sin,AA[2,3])
  p_DifNDVI=c(p_DifNDVI,AA[3,3])
  p_Troop=c(p_Troop,AA[4,3])
  p_Rank=c(p_Rank,AA[5,3])
  p_Parity=c(p_Parity,AA[6,3])
  p_Sex=c(p_Sex,AA[7,3])
  p_Group_Size=c(p_Group_Size,AA[8,3])
  p_Compet_Social=c(p_Compet_Social,AA[9,3])
  
  
}

# > mean(Estimate_DifNDVI)  #
# [1] -0.6627889
# > mean(Estimate_Compet_Soc)  #
# [1] 0.1589601
# > mean(Estimate_sin) #
# [1] -1.180106
# > mean(Estimate_Group_Size)  #
# [1] -0.08548049
# > mean(Estimate_TroopL)  #
# [1] -1.17803
# > mean(Estimate_TroopM) #
# [1] 0.6148467
# > mean(Estimate_Rank)  #
# [1] -0.5135817
# > mean(Estimate_Parity)  #
# [1] -0.4241944
# > mean(Estimate_Sex)  #
# [1] 0.2000266
# > 
#   > mean(CILOW_DifNDVI)  #
# [1] -1.241685
# > mean(CIUP_DifNDVI)  #
# [1] -0.08389264
# > mean(CILOW_Compet_Social)  #
# [1] -0.3081505
# > mean(CIUP_Compet_Social)  #
# [1] 0.6260708
# > mean(CILOW_sin)  #
# [1] -1.883936
# > mean(CIUP_sin)  #
# [1] -0.4762763
# > mean(CILOW_Group_Size)  #
# [1] -0.6528377
# > mean(CIUP_Group_Size)  #
# [1] 0.4818767
# > mean(CILOW_TroopL)  #
# [1] -2.241199
# > mean(CIUP_TroopL)  #
# [1] -0.1148604
# > mean(CILOW_TroopM)  #
# [1] -2.846721
# > mean(CIUP_TroopM)  #
# [1] 4.076415
# > mean(CILOW_Rank)  #
# [1] -1.015702
# > mean(CIUP_Rank)  #
# [1] -0.01146129
# > mean(CILOW_Parity)  #
# [1] -1.877992
# > mean(CIUP_Parity)  #
# [1] 1.029603
# > mean(CILOW_Sex) #
# [1] -0.7153621
# > mean(CIUP_Sex)  #
# [1] 1.115415
# > 
#   > mean(Chisq_DifNDVI)  #
# [1] 5.037682
# > mean(Chisq_Compet_Social)  #
# [1] 0.44798
# > mean(Chisq_sin)  #
# [1] 10.83922
# > mean(Chisq_Group_Size)  #
# [1] 0.08984627
# > mean(Chisq_Troop)  #
# [1] 5.143207
# > mean(Chisq_Rank)  #
# [1] 4.023335
# > mean(Chisq_Parity)  #
# [1] 0.3332665
# > mean(Chisq_Sex)  #
# [1] 0.1890407
# > 
#   > mean(p_DifNDVI)  #
# [1] 0.0250068
# > mean(p_Compet_Social)  #
# [1] 0.5055968
# > mean(p_sin)  #
# [1] 0.001694085
# > mean(p_Group_Size)  #
# [1] 0.7683711
# > mean(p_Troop)  #
# [1] 0.07705976
# > mean(p_Rank)  #
# [1] 0.04541731
# > mean(p_Parity)  #
# [1] 0.5684502
# > mean(p_Sex)  #
# [1] 0.6692348
# > 
#   > 

# Here we computed the 95% confidence interval of the pvalues of our fixed effect of interest (non-seasonal NDVI variation)
#   > wilcox.test(p_DifNDVI,conf.int = TRUE, conf.level = 0.95)
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  p_DifNDVI
# V = 500500, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   0.02470517 0.02511315
# sample estimates:
#   (pseudo)median 
# 0.0248796 
