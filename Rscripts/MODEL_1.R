library(lme4)
library(blme)
library(MuMIn)
library(DHARMa)
library(car)

### Scripts ran for Model 1, on interbirth intervals (IBI)

TAB=read.csv2('file_path of the excel tab IBI') #149 observations

#first we select the good IBI only, i.e. the ones for which the infant opening the IBI survived until weaning
#the infant closing the IBI is not born dead
#and no abortions/infants died in between the two infants of the IBI
TAB=TAB[which(TAB$DeathInBw=='0'),]
TAB=TAB[which(TAB$DeathFirstYear1=="n"),]
TAB=TAB[which(TAB$Borndead=='n'),] #on se retrouve avec 120 IBI



#####we scale all quantitative response variables: 

#our proxies of group reproductive synchrony:
TAB$Nb_Birth_After1=scale(TAB$Nb_Birth_After1)
TAB$Nb_Birth_After3=scale(TAB$Nb_Birth_After3)
TAB$Nb_Birth_After6=scale(TAB$Nb_Birth_After6)
TAB$Nb_Birth_Before1=scale(TAB$Nb_Birth_Before1)
TAB$Nb_Birth_Before3=scale(TAB$Nb_Birth_Before3)
TAB$Nb_Birth_Before6=scale(TAB$Nb_Birth_Before6)
TAB$Nb_Birth_Both1=scale(TAB$Nb_Birth_Both1)
TAB$Nb_Birth_Both3=scale(TAB$Nb_Birth_Both3)
TAB$Nb_Birth_Both6=scale(TAB$Nb_Birth_Both6)

TAB$Age_Year_Mother_CARRE=TAB$Age_Year_Mother*TAB$Age_Year_Mother #quadratic effect of age put in the model. 
TAB$Age_Year_Mother=scale(TAB$Age_Year_Mother)
TAB$Age_Year_Mother_CARRE=scale(TAB$Age_Year_Mother_CARRE)

TAB$NB_Adult_Females_V2=scale(TAB$NB_Adult_Females_V2)
TAB$Relative_Rank=scale(TAB$Relative_Rank)
TAB$Dif_NDVI=scale(TAB$Dif_NDVI)
TAB$Dif_Rainfall=scale(TAB$Dif_Rainfall)



#here we draw randomized birth dates according to the number of days of uncertainty for each birth date. 
#and store this new set of birth dates, in radians, in a vector merged in our table. 
j=1
for (i in 1:999) { 
  print(i)
  vec=c()
  for (j in 1:nrow(TAB)) { 
    ff=runif(1,TAB$Birth_Rad[j]-TAB$Uncertainty_Birth_DayRad[j],TAB$Birth_Rad[j]+TAB$Uncertainty_Birth_DayRad[j])
    
    if (ff>2*pi) {ff=ff-2*pi}
    if (ff<0) {ff=ff+2*pi}
    vec[j]=ff
  }
  TAB=cbind(TAB,vec)
}


###### I°) Determine the best sine phase effect on IBI: 

AIC0=c()
AIC1=c()
AIC2=c()
AIC3=c()
AIC4=c()
AIC5=c()


for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD0=lmer(IBIDays ~ sin(TAB[,j]+0*pi/6) + (1|Mother), data=TAB)
  MOD1=lmer(IBIDays ~ sin(TAB[,j]+1*pi/6) + (1|Mother), data=TAB)
  MOD2=lmer(IBIDays ~ sin(TAB[,j]+2*pi/6) + (1|Mother), data=TAB)
  MOD3=lmer(IBIDays ~ sin(TAB[,j]+3*pi/6) + (1|Mother), data=TAB)
  MOD4=lmer(IBIDays ~ sin(TAB[,j]+4*pi/6) + (1|Mother), data=TAB)
  MOD5=lmer(IBIDays ~ sin(TAB[,j]+5*pi/6) + (1|Mother), data=TAB)
  
  AIC0=c(AIC0,AIC(MOD0))
  AIC1=c(AIC1,AIC(MOD1))
  AIC2=c(AIC2,AIC(MOD2))
  AIC3=c(AIC3,AIC(MOD3))
  AIC4=c(AIC4,AIC(MOD4))
  AIC5=c(AIC5,AIC(MOD5))

}
mean(AIC0)
mean(AIC1)
mean(AIC2)
mean(AIC3)
mean(AIC4)
mean(AIC5)
#we select here the best phase as the one that minimize AIC, here AIC1 is minimum. 


###### II°) Determine the best non-seasonal variation to explain IBI variation. 

AIC1=c()
AIC2=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD1=lmer(IBIDays ~ sin(TAB[,j]+1*pi/6) + Dif_NDVI + (1|Mother), data=TAB)
  MOD2=lmer(IBIDays ~ sin(TAB[,j]+1*pi/6) + Dif_Rainfall + (1|Mother), data=TAB)
  
  AIC1=c(AIC1,AIC(MOD1))
  AIC2=c(AIC2,AIC(MOD2))
  
}

mean(AIC1)
mean(AIC2)
#we select the best one as the one minimizing AIC, here it is Dif_NDVI



###### III°) Determine the best group birth synchrony metric to explain IBI variation. 


MOD0=lmer(IBIDays ~ Nb_Birth_After1 + (1|Mother), data=TAB)
MOD1=lmer(IBIDays ~ Nb_Birth_After3 + (1|Mother), data=TAB)
MOD2=lmer(IBIDays ~ Nb_Birth_After6 + (1|Mother), data=TAB)
MOD3=lmer(IBIDays ~ Nb_Birth_Before1 + (1|Mother), data=TAB)
MOD4=lmer(IBIDays ~ Nb_Birth_Before3 + (1|Mother), data=TAB)
MOD5=lmer(IBIDays ~ Nb_Birth_Before6 + (1|Mother), data=TAB)
MOD6=lmer(IBIDays ~ Nb_Birth_Both1 + (1|Mother), data=TAB)
MOD7=lmer(IBIDays ~ Nb_Birth_Both3 + (1|Mother), data=TAB)
MOD8=lmer(IBIDays ~ Nb_Birth_Both6 + (1|Mother), data=TAB)
  
AIC(MOD0,MOD1,MOD2,MOD3,MOD4,MOD5,MOD6,MOD7,MOD8)
#we select the best group reprdductive synchrony as the one minimizing the AIC, here NB_Birth_Before3



#here we ran our full models with 3 interactions, and we detect if the 3 interactions of rank with resp.
#seasonal environmental variation, non-seasonal environmental variation, and group reproductive synchrony
#are significant or not

Estimate_Interaction_NDVI=c()
Estimate_Interaction_sin=c()
CILOW_Interaction_sin=c()
CIUP_Interaction_sin=c()
CILOW_Interaction_NDVI=c()
CIUP_Interaction_NDVI=c()
Chisq_Interaction_sin=c()
Chisq_Interaction_NDVI=c()
p_Interaction_sin=c()
p_Interaction_NDVI=c()


j=ncol(TAB)-999

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + Troop + Relative_Rank + Parity + Sex1 + 
             NB_Adult_Females_V2 +  Age_Year_Mother + Age_Year_Mother_CARRE +
             Nb_Birth_Before3 + (sin(TAB[,j]+pi/6) : Relative_Rank) + (Dif_NDVI:Relative_Rank) +
             (Nb_Birth_Before3:Relative_Rank) + (1|Mother), data=TAB)
  SS=summary(MOD)$coefficients
  CC=confint(MOD,method='Wald')
  AA=Anova(MOD,type=3)
  
  Estimate_Interaction_sin=c(Estimate_Interaction_sin,SS[13,1])
  Estimate_Interaction_NDVI=c(Estimate_Interaction_NDVI,SS[14,1])
  
  CILOW_Interaction_sin=c(CILOW_Interaction_sin,CC[15,1])
  CIUP_Interaction_sin=c(CIUP_Interaction_sin,CC[15,2])
  CILOW_Interaction_NDVI=c(CILOW_Interaction_NDVI,CC[16,1])
  CIUP_Interaction_NDVI=c(CIUP_Interaction_NDVI,CC[16,2])
  
  Chisq_Interaction_sin=c(Chisq_Interaction_sin,AA[12,1])
  Chisq_Interaction_NDVI=c(Chisq_Interaction_NDVI,AA[13,1])
  p_Interaction_sin=c(p_Interaction_sin,AA[12,3])
  p_Interaction_NDVI=c(p_Interaction_NDVI,AA[13,3])
  
  
}

mean(Estimate_Interaction_sin) #
mean(Estimate_Interaction_NDVI) #

mean(CILOW_Interaction_sin)  #
mean(CIUP_Interaction_sin)  #
mean(CILOW_Interaction_NDVI)  #
mean(CIUP_Interaction_NDVI) #

mean(Chisq_Interaction_sin) #
mean(Chisq_Interaction_NDVI) #

mean(p_Interaction_sin) #
mean(p_Interaction_NDVI) #



#Here we ran our final global model, taking into account birth date uncertainty
#by actually running 1000 models, with randomized birth dates drawn in the possible birth dates
#with the significant inetractions only, i.e. the one between female rank and group reproductive synchrony. 


Estimate_sin=c()
Estimate_DifNDVI=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Group_Size=c()
Estimate_Age1=c()
Estimate_Age2=c()
Estimate_Sex=c()
Estimate_Compet_Soc=c()
Estimate_Interaction=c()

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
CILOW_Age1=c()
CIUP_Age1=c()
CILOW_Age2=c()
CIUP_Age2=c()
CILOW_Compet_Social=c()
CIUP_Compet_Social=c()
CILOW_Interaction=c()
CIUP_Interaction=c()

Chisq_sin=c()
Chisq_DifNDVI=c()
Chisq_Troop=c()
Chisq_Rank=c()
Chisq_Parity=c()
Chisq_Sex=c()
Chisq_Group_Size=c()
Chisq_Age1=c()
Chisq_Age2=c()
Chisq_Compet_Social=c()
Chisq_Interaction=c()

p_sin=c()
p_DifNDVI=c()
p_Troop=c()
p_Rank=c()
p_Parity=c()
p_Group_Size=c()
p_Age1=c()
p_Age2=c()
p_Sex=c()
p_Compet_Social=c()
p_Interaction=c()

j=ncol(TAB)-999

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + Troop + Relative_Rank + Parity + Sex1 + 
             NB_Adult_Females_V2 +  Age_Year_Mother + Age_Year_Mother_CARRE +
             Nb_Birth_Before3 + 
             (Nb_Birth_Before3:Relative_Rank) + (1|Mother), data=TAB) #full model here, TAB[,j] is the drawn randomized birth dates of the infant opening the IBI
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
  Estimate_Age1=c(Estimate_Age1,SS[10,1])
  Estimate_Age2=c(Estimate_Age2,SS[11,1])
  Estimate_Compet_Soc=c(Estimate_Compet_Soc,SS[12,1])
  Estimate_Interaction=c(Estimate_Interaction,SS[13,1])
  
  CILOW_sin=c(CILOW_sin,CC[4,1])
  CIUP_sin=c(CIUP_sin,CC[4,2])
  CILOW_DifNDVI=c(CILOW_DifNDVI,CC[5,1])
  CIUP_DifNDVI=c(CIUP_DifNDVI,CC[5,2])
  CILOW_TroopL=c(CILOW_TroopL,CC[6,1])
  CIUP_TroopL=c(CIUP_TroopL,CC[6,2])
  CILOW_TroopM=c(CILOW_TroopM,CC[7,1])
  CIUP_TroopM=c(CIUP_TroopM,CC[7,2])
  CILOW_Rank=c(CILOW_Rank,CC[8,1])
  CIUP_Rank=c(CIUP_Rank,CC[8,2])
  CILOW_Parity=c(CILOW_Parity,CC[9,1])
  CIUP_Parity=c(CIUP_Parity,CC[9,2])
  CILOW_Sex=c(CILOW_Sex,CC[10,1])
  CIUP_Sex=c(CIUP_Sex,CC[10,2])
  CILOW_Group_Size=c(CILOW_Group_Size,CC[11,1])
  CIUP_Group_Size=c(CIUP_Group_Size,CC[11,2])
  CILOW_Age1=c(CILOW_Age1,CC[12,1])
  CIUP_Age1=c(CIUP_Age1,CC[12,2])
  CILOW_Age2=c(CILOW_Age2,CC[13,1])
  CIUP_Age2=c(CIUP_Age2,CC[13,2])
  CILOW_Compet_Social=c(CILOW_Compet_Social,CC[14,1])
  CIUP_Compet_Social=c(CIUP_Compet_Social,CC[14,2])
  CILOW_Interaction=c(CILOW_Interaction,CC[15,1])
  CIUP_Interaction=c(CIUP_Interaction,CC[15,2])
  
  Chisq_sin=c(Chisq_sin,AA[2,1])
  Chisq_DifNDVI=c(Chisq_DifNDVI,AA[3,1])
  Chisq_Troop=c(Chisq_Troop,AA[4,1])
  Chisq_Rank=c(Chisq_Rank,AA[5,1])
  Chisq_Sex=c(Chisq_Sex,AA[6,1])
  Chisq_Parity=c(Chisq_Parity,AA[7,1])
  Chisq_Group_Size=c(Chisq_Group_Size,AA[8,1])
  Chisq_Age1=c(Chisq_Age1,AA[9,1])
  Chisq_Age2=c(Chisq_Age2,AA[10,1])
  Chisq_Compet_Social=c(Chisq_Compet_Social,AA[11,1])
  Chisq_Interaction=c(Chisq_Interaction,AA[12,1])
  
  p_sin=c(p_sin,AA[2,3])
  p_DifNDVI=c(p_DifNDVI,AA[3,3])
  p_Troop=c(p_Troop,AA[4,3])
  p_Rank=c(p_Rank,AA[5,3])
  p_Parity=c(p_Parity,AA[6,3])
  p_Sex=c(p_Sex,AA[7,3])
  p_Group_Size=c(p_Group_Size,AA[8,3])
  p_Age1=c(p_Age1,AA[9,3])
  p_Age2=c(p_Age2,AA[10,3])
  p_Compet_Social=c(p_Compet_Social,AA[11,3])
  p_Interaction=c(p_Interaction,AA[12,3])
  
}


# > summary(Estimate_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -16.846  -7.605  -6.828  -6.841  -6.084  -3.094 
# > summary(Estimate_Compet_Soc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.52   16.33   17.00   16.98   17.65   20.38 
# > summary(Estimate_Interaction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -31.69  -28.89  -28.23  -28.26  -27.54  -23.43 
# > 
#   > summary(Estimate_sin)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -14.17   49.19   52.35   52.16   55.09   64.34 
# > summary(Estimate_Group_Size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.475  -3.347  -2.219  -2.208  -1.087   3.036 
# > summary(Estimate_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -70.19  -62.67  -60.90  -61.02  -59.24  -49.83 
# > summary(Estimate_TroopM)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -39.243 -23.493 -19.239 -19.411 -15.427   9.951 
# > summary(Estimate_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -36.88  -33.25  -32.65  -32.64  -32.02  -29.63 
# > summary(Estimate_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -24.961 -15.649 -13.119 -13.283 -10.899   2.314 
# > summary(Estimate_Sex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 30.32   35.92   37.65   37.57   39.22   44.84 
# > summary(Estimate_Age1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -224.5  -215.8  -213.2  -213.2  -210.8  -188.3 
# > summary(Estimate_Age2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 169.0   187.3   189.9   189.8   192.2   201.8 
# > 
#   > 
#   > summary(CILOW_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -40.03  -30.51  -29.76  -29.75  -28.95  -25.83 
# > summary(CIUP_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.34   15.29   16.10   16.07   16.80   19.65 
# > summary(CILOW_Compet_Social)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.1180 -4.3437 -3.6354 -3.6654 -3.0118 -0.3874 
# > summary(CIUP_Compet_Social)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.16   36.95   37.64   37.62   38.31   41.14 
# > summary(CILOW_Interaction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -52.77  -50.14  -49.44  -49.46  -48.76  -45.40 
# > summary(CIUP_Interaction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -10.613  -7.696  -7.020  -7.057  -6.328  -1.458 
# > 
#   > summary(CILOW_sin)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -142.09   18.23   21.40   21.17   24.28   33.02 
# > summary(CIUP_sin)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 71.75   80.17   83.12   83.15   86.03  113.75 
# > summary(CILOW_Group_Size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -31.10  -27.97  -26.77  -26.78  -25.65  -21.37 
# > summary(CIUP_Group_Size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.11   21.27   22.34   22.37   23.50   27.44 
# > summary(CILOW_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -130.8  -124.2  -122.7  -122.7  -121.3  -116.0 
# > summary(CIUP_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -9.6672 -1.7134  0.8960  0.6947  3.2397 16.8230 
# > summary(CILOW_TroopM)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -189.1  -172.7  -168.4  -168.5  -164.2  -145.5 
# > summary(CIUP_TroopM)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 110.6   125.9   129.7   129.7   133.6   165.4 
# > summary(CILOW_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -66.72  -61.21  -60.49  -60.51  -59.83  -56.64 
# > summary(CIUP_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -7.795  -5.442  -4.755  -4.762  -4.092  -2.075 
# > summary(CILOW_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -96.36  -86.95  -84.61  -84.79  -82.55  -76.41 
# > summary(CIUP_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.16   55.59   58.47   58.22   60.86   83.69 
# > summary(CILOW_Sex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -11.033  -5.055  -3.298  -3.400  -1.842   4.368 
# > summary(CIUP_Sex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 71.56   76.90   78.60   78.54   80.25   85.31 
# > summary(CILOW_Age1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -532.2  -361.0  -358.8  -359.0  -356.6  -347.9 
# > summary(CIUP_Age1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -81.24  -70.70  -67.53  -67.47  -64.55  155.64 
# > summary(CILOW_Age2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -129.98   48.00   51.06   50.83   53.94   62.78 
# > summary(CIUP_Age2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 318.1   326.4   328.8   328.7   330.9   467.9 
# > 
#   > 
#   > summary(Chisq_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07114 0.27277 0.34299 0.35183 0.42327 2.02789 
# > summary(Chisq_Compet_Social)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.648   2.396   2.608   2.608   2.807   3.699 
# > summary(Chisq_Interaction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.368   6.477   6.810   6.835   7.149   8.684 
# > 
#   > summary(Chisq_sin)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.04711  9.71110 10.96279 11.05224 12.38531 16.26930 
# > summary(Chisq_Group_Size)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.000e-08 8.849e-03 3.232e-02 4.747e-02 7.054e-02 2.656e-01 
# > summary(Chisq_Troop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.308   3.479   3.756   3.795   4.079   5.204 
# > summary(Chisq_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.416   5.053   5.262   5.273   5.493   6.289 
# > summary(Chisq_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.065   2.956   3.246   3.244   3.504   4.716 
# > summary(Chisq_Sex)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003106 0.088332 0.129023 0.142339 0.186941 0.473232 
# > summary(Chisq_Age1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.151   7.978   8.234   8.244   8.506   9.435 
# > summary(Chisq_Age2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.227   6.946   7.181   7.184   7.417   8.128 
# > 
#   > 
#   > summary(p_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1544  0.5153  0.5581  0.5603  0.6015  0.7897 
# > summary(p_Compet_Social)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05443 0.09385 0.10633 0.10856 0.12161 0.19916 
# > summary(p_Interaction)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003210 0.007501 0.009063 0.009332 0.010929 0.036618 
# > 
#   > summary(p_sin)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000549 0.0004327 0.0009296 0.0022085 0.0018316 0.8281689 
# > summary(p_Group_Size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6063  0.7906  0.8573  0.8537  0.9251  0.9999 
# > summary(p_Troop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07413 0.13010 0.15292 0.15295 0.17559 0.31533 
# > summary(p_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01215 0.01909 0.02179 0.02205 0.02458 0.03560 
# > summary(p_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4915  0.6655  0.7194  0.7167  0.7663  0.9556 
# > summary(p_Sex)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02989 0.06122 0.07160 0.07408 0.08554 0.15071 
# > summary(p_Age1)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002129 0.003540 0.004112 0.004437 0.004736 0.283258 
# > summary(p_Age2)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.004359 0.006460 0.007366 0.007721 0.008403 0.267936 
# > 


#here we finally computed the 95% confidence intervals of the pvalue of interest (i.e. the interaction between female rank and group reproductive synchrony)
wilcox.test(p_Interaction,conf.int = TRUE, conf.level = 0.95)

# 95 percent confidence interval:
#   0.009055183 0.009365206
# sample estimates:
#   (pseudo)median 
# 0.00918441






## lastly, we run model diagnostics:


MOD=lmer(IBIDays ~ sin(Birth_Rad+pi/6) + Dif_NDVI + Troop + Relative_Rank + Parity + Sex1 + 
                 NB_Adult_Females_V2 +  Age_Year_Mother + Age_Year_Mother_CARRE +
                 Nb_Birth_Before3 + 
                 (Nb_Birth_Before3:Relative_Rank) + (1|Mother), data=TAB) 

#if singular (not the case here, we ran the same model with the function 'blmer')

r.squaredGLMM(MOD)
vif(MOD) #only the 2 fixed effects of female age have vif >2.5, which was expected. ==> no collinearity
qqPlot(residuals(MOD))





