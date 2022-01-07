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

###we scale the quantitative fixed effects:
# female rank
TAB$Relative_Rank=scale(TAB$Relative_Rank, center=TRUE, scale=TRUE)
# non seasonal NDVI 
TAB$Dif_NDVI=scale(TAB$Dif_NDVI)
# non seasonal rainfall 
TAB$Dif_Rainfall=scale(TAB$Dif_Rainfall)
# our metrics of group reproductive synchrony: number of infants born before, around and after the focal birth                   
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
# number of adult females in the group (proxy of group size)
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

#Following the exact same procedure as the one explained for MODEL_1 (IBI), we determined : 
# the best sine phase, which is pi/2 here
# the best non-seasonal environmental variation, i.e. NDVI_NS
# the best group size effect, i.e. simple effect (not quadratic)
# the best group reproductive synchrony metric, i.e. number of infants born 1 month after the focal birth


#here we run the full model, without any female rank interactions as these interactions were not significant... 
#taking into account birth date uncertainty. and applying Rubin's rule to compute pooled estimate, SD and pvalues

Estimate_NDVI_Seasonality=c()
Estimate_DifNDVI=c()
Estimate_Nb_Females=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Sex=c()
Esimate_Compet=c()

VAR_NDVI_Seasonality=c()
VAR_DifNDVI=c()
VAR_Nb_Females=c()
VAR_TroopL=c()
VAR_TroopM=c()
VAR_Rank=c()
VAR_Parity=c()
VAR_Sex=c()
VAR_Compet=c()

p_NDVI_Seasonality=c()
p_DifNDVI=c()
p_Nb_Females=c()
p_TroopL=c()
p_TroopM=c()
p_Rank=c()
p_Parity=c()
p_Sex=c()
p_Compet=c()

for (j in (ncol(TAB)-999):ncol(TAB)) {
  
  print(j)
  
  
  MOD=glmer(Death_18months ~ sin(TAB[,j]+pi/2) + Dif_NDVI + NB_Adult_Females_V2 +
              Troop + Relative_Rank + Parity + 
              Sex +  Nb_Birth_After1 + (1|Mother), 
            glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
  
  
  SS=summary(MOD)$coefficients
  VV=vcov(MOD)
  
  Estimate_NDVI_Seasonality=c(Estimate_NDVI_Seasonality,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_Nb_Females=c(Estimate_Nb_Females,SS[4,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[5,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[6,1])
  Estimate_Rank=c(Estimate_Rank,SS[7,1])
  Estimate_Parity=c(Estimate_Parity,SS[8,1])
  Estimate_Sex=c(Estimate_Sex,SS[9,1])
  Esimate_Compet=c(Esimate_Compet,SS[10,1])
  
  VAR_NDVI_Seasonality=c(VAR_NDVI_Seasonality,VV[2, 2])
  VAR_DifNDVI=c(VAR_DifNDVI,VV[3, 3])
  VAR_Nb_Females=c(VAR_Nb_Females,VV[4, 4])
  VAR_TroopL=c(VAR_TroopL,VV[5, 5])
  VAR_TroopM=c(VAR_TroopM,VV[6, 6])
  VAR_Rank=c(VAR_Rank,VV[7, 7])
  VAR_Parity=c(VAR_Parity,VV[8, 8])
  VAR_Sex=c(VAR_Sex,VV[9, 9])
  VAR_Compet=c(VAR_Compet,VV[10, 10])
  
  p_NDVI_Seasonality=c(p_NDVI_Seasonality,SS[2,4])
  p_DifNDVI=c(p_DifNDVI,SS[3,4])
  p_Nb_Females=c(p_Nb_Females,SS[4,4])
  p_TroopL=c(p_TroopL,SS[5,4])
  p_TroopM=c(p_TroopM,SS[6,4])
  p_Rank=c(p_Rank,SS[7,4])
  p_Parity=c(p_Parity,SS[8,4])
  p_Sex=c(p_Sex,SS[9,4])
  p_Compet=c(p_Compet,SS[10,4])
  
}

m=1000 # m ici est le nombre de imputed dataset, cad de simulations, donc ici 1000
n=1000
k=length(coef(MOD))

# DIF NDVI: non seasonal environmental variation
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

# Group reproductive synchrony
MAT=data.frame(Esimate_Compet,VAR_Compet,p_Compet)  
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


# Sine term of birth date (seasonality)
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


# Rank
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



# Parity primi
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





