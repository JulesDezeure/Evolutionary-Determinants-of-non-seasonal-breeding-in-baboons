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

TAB$Mother=factor(TAB$Mother)


#####we scale all quantitative response variables: 

#our proxies of group reproductive synchrony:
TAB$NB_Births_15days_Around=scale(TAB$NB_Births_15days_Around)
TAB$Nb_Birth_After1=scale(TAB$Nb_Birth_After1)
TAB$Nb_Birth_Before1=scale(TAB$Nb_Birth_Before1)

TAB$Nb_Birth_Both1=scale(TAB$Nb_Birth_Both1)
TAB$NB_Births_2months_before=scale(TAB$NB_Births_2months_before)
TAB$NB_Births_2months_after=scale(TAB$NB_Births_2months_after)

TAB$NB_Births_2months_Around=scale(TAB$NB_Births_2months_Around)
TAB$NB_Births_4months_before=scale(TAB$NB_Births_4months_before)
TAB$NB_Births_4months_after=scale(TAB$NB_Births_4months_after)

TAB$Nb_Birth_Both3=scale(TAB$Nb_Birth_Both3)
TAB$Nb_Birth_Before6=scale(TAB$Nb_Birth_Before6)
TAB$Nb_Birth_After6=scale(TAB$Nb_Birth_After6)


#female age
TAB$Age_Year_Mother_CARRE=TAB$Age_Year_Mother*TAB$Age_Year_Mother #quadratic effect of female age 
TAB$Age_Year_Mother=scale(TAB$Age_Year_Mother)
TAB$Age_Year_Mother_CARRE=scale(TAB$Age_Year_Mother_CARRE)

#group soze (number of adult females in the group)
TAB$NB_Adult_Females_V2carre=TAB$NB_Adult_Females_V2*TAB$NB_Adult_Females_V2 # quadratic effect of group size
TAB$NB_Adult_Females_V2carre=scale(TAB$NB_Adult_Females_V2carre)
TAB$NB_Adult_Females_V2=scale(TAB$NB_Adult_Females_V2)

#female rank
TAB$Relative_Rank=scale(TAB$Relative_Rank)

#non seasonal environmental variation
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


### III°) Determine if we would put a simple or a quadratic effect of group size in the full model. 

AICGroup1=c()
AICGroup2=c()
for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MODGroup1=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + Troop + Relative_Rank + Parity + Sex1 + 
                   NB_Adult_Females_V2 + Age_Year_Mother + Age_Year_Mother_CARRE +
                   (1|Mother), data=TAB)
  MODGroup2=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + Troop + Relative_Rank + Parity + Sex1 + 
                   NB_Adult_Females_V2 + NB_Adult_Females_V2carre +  Age_Year_Mother + Age_Year_Mother_CARRE +
                   (1|Mother), data=TAB)
  AICGroup1=c(AICGroup1,AIC(MODGroup1))
  AICGroup2=c(AICGroup2,AIC(MODGroup2))
  
}

mean(AICGroup1) #
mean(AICGroup2) #

#AIC group2 is lower so we put the quadratic effect. 


###### IV°) Determine the best group birth synchrony metric to explain IBI variation. 



AIC0=c()
AIC1=c()
AIC2=c()
AIC3=c()
AIC4=c()
AIC5=c()
AIC6=c()
AIC7=c()
AIC8=c()
AIC9=c()
AIC10=c()
AIC11=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD0=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_15days_Around + 
              (1|Mother), data=TAB)
  MOD1=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              Nb_Birth_Before1 + 
              (1|Mother), data=TAB)
  MOD2=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              Nb_Birth_After1 + 
              (1|Mother), data=TAB)
  MOD3=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              Nb_Birth_Both1 + 
              (1|Mother), data=TAB)
  MOD4=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_2months_before +
              (1|Mother), data=TAB)
  MOD5=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_2months_after + 
              (1|Mother), data=TAB)
  MOD6=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_2months_Around + 
              (1|Mother), data=TAB)
  MOD7=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_4months_before + 
              (1|Mother), data=TAB)
  MOD8=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              NB_Births_4months_after + 
              (1|Mother), data=TAB)
  MOD9=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
              Troop + Relative_Rank + Parity + Sex1 + 
              Age_Year_Mother + Age_Year_Mother_CARRE +
              Nb_Birth_Both3 + 
              (1|Mother), data=TAB)
  MOD10=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
               Troop + Relative_Rank + Parity + Sex1 + 
               Age_Year_Mother + Age_Year_Mother_CARRE +
               Nb_Birth_Before6 +
               (1|Mother), data=TAB)
  MOD11=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
               Troop + Relative_Rank + Parity + Sex1 + 
               Age_Year_Mother + Age_Year_Mother_CARRE +
               Nb_Birth_After6 +
               (1|Mother), data=TAB)
  
  AIC0=c(AIC0,AIC(MOD0))
  AIC1=c(AIC1,AIC(MOD1))
  AIC2=c(AIC2,AIC(MOD2))
  AIC3=c(AIC3,AIC(MOD3))
  AIC4=c(AIC4,AIC(MOD4))
  AIC5=c(AIC5,AIC(MOD5))
  AIC6=c(AIC6,AIC(MOD6))
  AIC7=c(AIC7,AIC(MOD7))
  AIC8=c(AIC8,AIC(MOD8))
  AIC9=c(AIC9,AIC(MOD9))
  AIC10=c(AIC10,AIC(MOD10))
  AIC11=c(AIC11,AIC(MOD11))
  
  
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2) 
mean(AIC3) 
mean(AIC4) 
mean(AIC5) 
mean(AIC6) 
mean(AIC7) 
mean(AIC8) 
mean(AIC9) 
mean(AIC10) 
mean(AIC11) 

#we select the best group reprdductive synchrony as the one minimizing the AIC, here AIC6 (2 months around)


## V°) FULL MODELS

# for these full models (and all the ones therafter), we applied Rubin's rule to extract
# pooled estimates, pooled SD and then to compute their confidence intervals and Pvalues associated. 

#here we ran our full models with 3 interactions, and we detect if the 3 interactions of rank with resp.
#seasonal environmental variation, non-seasonal environmental variation, and group reproductive synchrony
#are significant or not


Estimate_Interaction_Rank_Synchro=c()
Estimate_Interaction_Rank_Sinus=c()
Estimate_Interaction_Rank_DifNDVI=c()
VAR_Interaction_Rank_Synchro=c()
VAR_Interaction_Rank_Sinus=c()
VAR_Interaction_Rank_DifNDVI=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  MOD=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
             Troop + Relative_Rank + Parity + Sex1 + 
             Age_Year_Mother + Age_Year_Mother_CARRE +
             Nb_Birth_Before3 + (Nb_Birth_Before3:Relative_Rank) +
             ((sin(TAB[,j]+pi/6)):Relative_Rank) + (Dif_NDVI:Relative_Rank) +
             (1|Mother), data=TAB)
  
  SS=summary(MOD)$coefficients
  VV=vcov(MOD)
  
  Estimate_Interaction_Rank_Synchro=c(Estimate_Interaction_Rank_Synchro,SS[14,1])
  Estimate_Interaction_Rank_Sinus=c(Estimate_Interaction_Rank_Sinus,SS[15,1])
  Estimate_Interaction_Rank_DifNDVI=c(Estimate_Interaction_Rank_DifNDVI,SS[16,1])
  
  VAR_Interaction_Rank_Synchro=c(VAR_Interaction_Rank_Synchro,VV[14, 14])
  VAR_Interaction_Rank_Sinus=c(VAR_Interaction_Rank_Sinus,VV[15, 15])
  VAR_Interaction_Rank_DifNDVI=c(VAR_Interaction_Rank_DifNDVI,VV[16, 16])
  
}

#from the 1000 estimated and covariance matrix values extracted above, 
# we applied the Rubin's rule manually as follow: 
m=1000 
n=1000
k=length(coef(MOD))
#Interaction rank:synchronie
MAT=data.frame(Estimate_Interaction_Rank_Synchro,VAR_Interaction_Rank_Synchro)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))

if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VECVEC

#Interaction rank:sinus
MAT=data.frame(Estimate_Interaction_Rank_Sinus,VAR_Interaction_Rank_Sinus)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))

if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VECVEC

#interaction rank:NDVI_NS
MAT=data.frame(Estimate_Interaction_Rank_DifNDVI,VAR_Interaction_Rank_DifNDVI)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VECVEC


#we saw that the only significant interaction was the one between rank and reproductive synchrony, 
# so we only kept this one in the full model. We considered birth date uncertainty by drawing randomely
# 1000 birth dates in the possible birth date range for each IBI value, 
# as follow (the results are presented Table 2 of the main text): 



Estimate_NDVI_Seasonality=c()
Estimate_DifNDVI=c()
Estimate_Nb_Females=c()
Estimate_Nb_Females_Carre=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Sex=c()
Estimate_Age=c()
Estimate_Age_Carre=c()
Estimate_Compet=c()
Estimate_Interaction_Rank_Compet=c()

VAR_NDVI_Seasonality=c()
VAR_DifNDVI=c()
VAR_Nb_Females=c()
VAR_Nb_Females_Carre=c()
VAR_TroopL=c()
VAR_TroopM=c()
VAR_Rank=c()
VAR_Parity=c()
VAR_Sex=c()
VAR_Age=c()
VAR_Age_Carre=c()
VAR_Compet=c()
VAR_Interaction_Rank_Compet=c()

p_NDVI_Seasonality=c()
p_DifNDVI=c()
p_Nb_Females=c()
p_Nb_Females_Carre=c()
p_TroopL=c()
p_TroopM=c()
p_Rank=c()
p_Parity=c()
p_Sex=c()
p_Age=c()
p_Age_Carre=c()
p_Compet=c()
p_Interaction_Rank_Compet=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  
  MOD=lmer(IBIDays ~ sin(TAB[,j]+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
             Troop + Relative_Rank + Parity + Sex1 + 
             Age_Year_Mother + Age_Year_Mother_CARRE +
             NB_Births_2months_Around + (NB_Births_2months_Around:Relative_Rank) +
             (1|Mother), data=TAB)
  SS=summary(MOD)$coefficients
  VV=vcov(MOD)
  SS
  
  Estimate_NDVI_Seasonality=c(Estimate_NDVI_Seasonality,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_Nb_Females=c(Estimate_Nb_Females,SS[4,1])
  Estimate_Nb_Females_Carre=c(Estimate_Nb_Females_Carre,SS[5,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[6,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[7,1])
  Estimate_Rank=c(Estimate_Rank,SS[8,1])
  Estimate_Parity=c(Estimate_Parity,SS[9,1])
  Estimate_Sex=c(Estimate_Sex,SS[10,1])
  Estimate_Age=c(Estimate_Age,SS[11,1])
  Estimate_Age_Carre=c(Estimate_Age_Carre,SS[12,1])
  Estimate_Compet=c(Estimate_Compet,SS[13,1])
  Estimate_Interaction_Rank_Compet=c(Estimate_Interaction_Rank_Compet,SS[14,1])
  
  VAR_NDVI_Seasonality=c(VAR_NDVI_Seasonality,VV[2, 2])
  VAR_DifNDVI=c(VAR_DifNDVI,VV[3, 3])
  VAR_Nb_Females=c(VAR_Nb_Females,VV[4, 4])
  VAR_Nb_Females_Carre=c(VAR_Nb_Females_Carre,VV[5, 5])
  VAR_TroopL=c(VAR_TroopL,VV[6, 6])
  VAR_TroopM=c(VAR_TroopM,VV[7, 7])
  VAR_Rank=c(VAR_Rank,VV[8, 8])
  VAR_Parity=c(VAR_Parity,VV[9, 9])
  VAR_Sex=c(VAR_Sex,VV[10, 10])
  VAR_Age=c(VAR_Age,VV[11, 11])
  VAR_Age_Carre=c(VAR_Age_Carre,VV[12, 12])
  VAR_Compet=c(VAR_Compet,VV[13, 13])
  VAR_Interaction_Rank_Compet=c(VAR_Interaction_Rank_Compet,VV[14, 14])
  
  p_NDVI_Seasonality=c(p_NDVI_Seasonality,SS[2,4])
  p_DifNDVI=c(p_DifNDVI,SS[3,4])
  p_Nb_Females=c(p_Nb_Females,SS[4,4])
  p_Nb_Females_Carre=c(p_Nb_Females_Carre,SS[5,4])
  p_TroopL=c(p_TroopL,SS[6,4])
  p_TroopM=c(p_TroopM,SS[7,4])
  p_Rank=c(p_Rank,SS[8,4])
  p_Parity=c(p_Parity,SS[9,4])
  p_Sex=c(p_Sex,SS[10,4])
  p_Age=c(p_Age,SS[11,4])
  p_Age_Carre=c(p_Age_Carre,SS[12,4])
  p_Compet=c(p_Compet,SS[13,4])
  p_Interaction_Rank_Compet=c(p_Interaction_Rank_Compet,SS[14,4])
  
  
}



m=1000 #je pense que m ici est le nombre de imputed dataset, cad de simulations, donc ici 1000
n=1000
k=length(coef(MOD))

# NDVI_NS (non-seasonal environmental variation)
MAT=data.frame(Estimate_DifNDVI,VAR_DifNDVI,p_DifNDVI)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))

if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}

VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Non-seasonal environmental variation', VECVEC)
tableau=data.frame(VEC)
tableau=as.data.frame(t(tableau))
colnames(tableau) <- c("Fixed effects","Estimates",'SD Estimates','P-values')

# Reproductive synchrony
MAT=data.frame(Estimate_Compet,VAR_Compet,p_Compet)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Reproductive synchrony', VECVEC)
tableau=rbind(tableau,VEC)
VECVEC

# SYNCHRONY:Rank interaction
MAT=data.frame(Estimate_Interaction_Rank_Compet,VAR_Interaction_Rank_Compet,p_Interaction_Rank_Compet)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Reproductive synchrony : Rank', VECVEC)
tableau=rbind(tableau,VEC)


# NDVI SEASONALITY
MAT=data.frame(Estimate_NDVI_Seasonality,VAR_NDVI_Seasonality,p_NDVI_Seasonality)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Seasonal environmental variation', VECVEC)
tableau=rbind(tableau,VEC)



# Number of adult females
MAT=data.frame(Estimate_Nb_Females,VAR_Nb_Females,p_Nb_Females)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Number of adult females', VECVEC)
tableau=rbind(tableau,VEC)



# Number of adult females (quadratic)
MAT=data.frame(Estimate_Nb_Females_Carre,VAR_Nb_Females_Carre,p_Nb_Females_Carre)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Number of adult females carre', VECVEC)
tableau=rbind(tableau,VEC)



# GROUP L
MAT=data.frame(Estimate_TroopL,VAR_TroopL,p_TroopL)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Group (L)', VECVEC)
tableau=rbind(tableau,VEC)


# GROUP M
MAT=data.frame(Estimate_TroopM,VAR_TroopM,p_TroopM)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Group (M)', VECVEC)
tableau=rbind(tableau,VEC)


# Female rank
MAT=data.frame(Estimate_Rank,VAR_Rank,p_Rank)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Rank', VECVEC)
tableau=rbind(tableau,VEC)



# Parity primiparous
MAT=data.frame(Estimate_Parity,VAR_Parity,p_Parity)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Parity (primiparous)', VECVEC)
tableau=rbind(tableau,VEC)


# Sex male
MAT=data.frame(Estimate_Sex,VAR_Sex,p_Sex)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
#3/ pooled pvalues: carrément difficile
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Sex (Male)', VECVEC)
tableau=rbind(tableau,VEC)


# Female age (simple)
MAT=data.frame(Estimate_Age,VAR_Age,p_Age)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Age', VECVEC)
tableau=rbind(tableau,VEC)


# Female age (quadratic)
MAT=data.frame(Estimate_Age_Carre,VAR_Age_Carre,p_Age_Carre)  
pooledMean<- mean(MAT[,1])
(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
lambda <- (withinVar + (withinVar/m))/totVar
nu_old <- (m-1)/lambda^2  
nu_com <- n-k
nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
(nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs))
if (pooledMean > 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = FALSE) * 2
}
if (pooledMean < 0) {
  pooledP=pt(q = pooledMean / pooledSE, df = nu_BR, lower.tail = TRUE) * 2
}
VECVEC=round(c(pooledMean, pooledSE, pooledP),5)
VEC=c('Age carre', VECVEC)
tableau=rbind(tableau,VEC)

tableau



## lastly, we run model diagnostics:

MOD=lmer(IBIDays ~ sin(Birth_Rad+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
           Troop + Relative_Rank + Parity + Sex1 + 
           Age_Year_Mother + Age_Year_Mother_CARRE +
           NB_Births_2months_Around + (NB_Births_2months_Around:Relative_Rank) +
           (1|Mother), data=TAB)

#if singular (not the case here, we ran the same model with the function 'blmer')

r.squaredGLMM(MOD)
vif(MOD) #only the 2 fixed effects of female age have vif >2.5, which was expected. ==> no collinearity
qqPlot(residuals(MOD))





