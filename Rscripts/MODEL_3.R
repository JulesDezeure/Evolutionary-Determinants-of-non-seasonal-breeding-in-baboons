library(lubridate)
library(lme4)
library(blme)
library(MuMIn)
library(DHARMa)
library(car)


#The first step is to take into account the uncertainty around cycle resumption dates
#to do so, we ran 1000 randomized dates, drawn in the possible dates

TABFPPO=read.csv2('file_path of the excel tab Cycle resumption')


I=-250 #initialization (can be another value that -250)
for (k in 0:999) { # 1000 randomisations
  k=k+1
  print(k)
  
  vec=c()
  
  for (i in 1:(nrow(TAB)) ) { 
    
    ###### 1st case: we have an uncertain cycle resumption date (with uncertainty > 0 day)
    if (is.na(TAB$FPPO[i])==FALSE && TAB$FPPO[i]=="1" && TAB$FPPO_Uncertainty[i]>0) { 
      I=i    #on note la ligne de cette incertitude I
      NEW_DATE=TAB$FPPO_Date[i]+runif(1,-TAB$FPPO_Uncertainty[i]/2,TAB$FPPO_Uncertainty[i]/2) #here we draw a random date
      
      j=I-5
      for (j in (I-5):(I+5) ) { 
        #Here we move within a smaller time window, +- 5 months around the initial cycle resumption month
        
        #case 1 A: there was a O
        if (is.na(TAB$FPPO[j])==FALSE && TAB$FPPO[j]=="0") {
          if ( TAB$Month[j] == month(NEW_DATE) ) {vec[j]="1"} 
          if(!TAB$Month[j] == month(NEW_DATE)) {vec[j]='0'} 
        }
        #case 1B: there was a 1
        if (is.na(TAB$FPPO[j])==FALSE && TAB$FPPO[j]=="1") {
          if ( TAB$Month[j] == month(NEW_DATE) ) {vec[j]="1"} 
          if(!TAB$Month[j] == month(NEW_DATE)) {vec[j]='0'} 
        }
        #case C: it was a NA
        if (is.na(TAB$FPPO[j])==TRUE) {
          if ( TAB$Month[j] == month(NEW_DATE) ) {vec[j]="1"} 
          if(!TAB$Month[j] == month(NEW_DATE)) {vec[j]=NA}  
        }
        j=j+1
        
      }
      i=i+5
    }
    
    ###### 2nde case: if we are after an uncertain cycle resumption, we don't want to delete what we just did, so we go 5 rows later
    if ((i-I)<5) {i=i+5}
    
    ###### 3rd case: No cycle resumption, or an accurate cycle resumption date (with 0 day of uncertainty)
    
    else { vec[i]=TAB$FPPO[i] 
    i=i+1}
  }
  
  TAB=cbind(TAB,vec)
}

#that's how we obtained the 1000 clumns named 'vec' in the Table 'cycle resumption'





### FIRST, we select the various environmental fixed effects to retain for the full model: 


TABFPPO=read.csv2('file_path of the excel tab Cycle resumption')

#we scale the numerous environmental predictors:
TABFPPO$Rainfall_Seasonality0=scale(TABFPPO$Rainfall_Seasonality0)
TABFPPO$Rainfall_Seasonality1=scale(TABFPPO$Rainfall_Seasonality1)
TABFPPO$Rainfall_Seasonality2=scale(TABFPPO$Rainfall_Seasonality2)
TABFPPO$Rainfall_Seasonality3=scale(TABFPPO$Rainfall_Seasonality3)
TABFPPO$Rainfall_Seasonality4=scale(TABFPPO$Rainfall_Seasonality4)
TABFPPO$Rainfall_Seasonality5=scale(TABFPPO$Rainfall_Seasonality5)
TABFPPO$Rainfall_Seasonality6=scale(TABFPPO$Rainfall_Seasonality6)
TABFPPO$Rainfall_Seasonality7=scale(TABFPPO$Rainfall_Seasonality7)
TABFPPO$Rainfall_Seasonality8=scale(TABFPPO$Rainfall_Seasonality8)
TABFPPO$Rainfall_Seasonality9=scale(TABFPPO$Rainfall_Seasonality9)
TABFPPO$Rainfall_Seasonality10=scale(TABFPPO$Rainfall_Seasonality10)
TABFPPO$Rainfall_Seasonality11=scale(TABFPPO$Rainfall_Seasonality11)
TABFPPO$Rainfall_Seasonality12=scale(TABFPPO$Rainfall_Seasonality12)

TABFPPO$Dif_Rainfall0=scale(TABFPPO$Dif_Rainfall0)
TABFPPO$Dif_Rainfall1=scale(TABFPPO$Dif_Rainfall1)
TABFPPO$Dif_Rainfall2=scale(TABFPPO$Dif_Rainfall2)
TABFPPO$Dif_Rainfall3=scale(TABFPPO$Dif_Rainfall3)
TABFPPO$Dif_Rainfall4=scale(TABFPPO$Dif_Rainfall4)
TABFPPO$Dif_Rainfall5=scale(TABFPPO$Dif_Rainfall5)
TABFPPO$Dif_Rainfall6=scale(TABFPPO$Dif_Rainfall6)
TABFPPO$Dif_Rainfall7=scale(TABFPPO$Dif_Rainfall7)
TABFPPO$Dif_Rainfall8=scale(TABFPPO$Dif_Rainfall8)
TABFPPO$Dif_Rainfall9=scale(TABFPPO$Dif_Rainfall9)
TABFPPO$Dif_Rainfall10=scale(TABFPPO$Dif_Rainfall10)
TABFPPO$Dif_Rainfall11=scale(TABFPPO$Dif_Rainfall11)
TABFPPO$Dif_Rainfall12=scale(TABFPPO$Dif_Rainfall12)

TABFPPO$NDVI_Seasonality0=scale(TABFPPO$NDVI_Seasonality0)
TABFPPO$NDVI_Seasonality1=scale(TABFPPO$NDVI_Seasonality1)
TABFPPO$NDVI_Seasonality2=scale(TABFPPO$NDVI_Seasonality2)
TABFPPO$NDVI_Seasonality3=scale(TABFPPO$NDVI_Seasonality3)
TABFPPO$NDVI_Seasonality4=scale(TABFPPO$NDVI_Seasonality4)
TABFPPO$NDVI_Seasonality5=scale(TABFPPO$NDVI_Seasonality5)
TABFPPO$NDVI_Seasonality6=scale(TABFPPO$NDVI_Seasonality6)
TABFPPO$NDVI_Seasonality7=scale(TABFPPO$NDVI_Seasonality7)
TABFPPO$NDVI_Seasonality8=scale(TABFPPO$NDVI_Seasonality8)
TABFPPO$NDVI_Seasonality9=scale(TABFPPO$NDVI_Seasonality9)
TABFPPO$NDVI_Seasonality10=scale(TABFPPO$NDVI_Seasonality10)
TABFPPO$NDVI_Seasonality11=scale(TABFPPO$NDVI_Seasonality11)
TABFPPO$NDVI_Seasonality12=scale(TABFPPO$NDVI_Seasonality12)

TABFPPO$Dif_NDVI0=scale(TABFPPO$Dif_NDVI0)
TABFPPO$Dif_NDVI1=scale(TABFPPO$Dif_NDVI1)
TABFPPO$Dif_NDVI2=scale(TABFPPO$Dif_NDVI2)
TABFPPO$Dif_NDVI3=scale(TABFPPO$Dif_NDVI3)
TABFPPO$Dif_NDVI4=scale(TABFPPO$Dif_NDVI4)
TABFPPO$Dif_NDVI5=scale(TABFPPO$Dif_NDVI5)
TABFPPO$Dif_NDVI6=scale(TABFPPO$Dif_NDVI6)
TABFPPO$Dif_NDVI7=scale(TABFPPO$Dif_NDVI7)
TABFPPO$Dif_NDVI8=scale(TABFPPO$Dif_NDVI8)
TABFPPO$Dif_NDVI9=scale(TABFPPO$Dif_NDVI9)
TABFPPO$Dif_NDVI10=scale(TABFPPO$Dif_NDVI10)
TABFPPO$Dif_NDVI11=scale(TABFPPO$Dif_NDVI11)
TABFPPO$Dif_NDVI12=scale(TABFPPO$Dif_NDVI12)

TAB$FPPO=as.factor(TAB$FPPO)
TABFPPO$Relative.rank=scale(TABFPPO$Relative.rank, center=TRUE, scale=TRUE)
TABFPPO$NB_Adult_Females_V2=as.numeric(scale(TABFPPO$NB_Adult_Females_V2))
TABFPPO$Number_Conc_Update0=scale(TABFPPO$Number_Conc_Update0)
TABFPPO$Number_Conc_Update1=scale(TABFPPO$Number_Conc_Update1)
TABFPPO$Number_Conc_Update2=scale(TABFPPO$Number_Conc_Update2)
TABFPPO$Number_Conc_Update3=scale(TABFPPO$Number_Conc_Update3)
TABFPPO$Number_Conc_Update4=scale(TABFPPO$Number_Conc_Update4)
TABFPPO$Number_Conc_Update5=scale(TABFPPO$Number_Conc_Update5)
TABFPPO$Number_Conc_Update6=scale(TABFPPO$Number_Conc_Update6)
TABFPPO$NB_Cyc_Precis0=as.numeric(scale(TABFPPO$NB_Cyc_Precis0))


#here we select the best sine phase for the sine term fixed effect
AIC0=c()
AIC1=c()
AIC2=c()
AIC3=c()
AIC4=c()
AIC5=c()

for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)
  
  MOD0=bglmer(y ~ sin(Month_Rad+0*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ sin(Month_Rad+1*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ sin(Month_Rad+2*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ sin(Month_Rad+3*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ sin(Month_Rad+4*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ sin(Month_Rad+5*pi/6) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)#the best model: minimizing AIC
mean(AIC3) 
mean(AIC4) 
mean(AIC5)



#here we select the best time window for rainfall seasonality (i.e. Rain_S)
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
AIC12=c()

for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)

  MOD0=bglmer(y ~ Rainfall_Seasonality0 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ Rainfall_Seasonality1 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ Rainfall_Seasonality2 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ Rainfall_Seasonality3 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ Rainfall_Seasonality4 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ Rainfall_Seasonality5 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
  MOD6=bglmer(y ~ Rainfall_Seasonality6 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC6=c(AIC6,AIC(MOD6))
  
  MOD7=bglmer(y ~ Rainfall_Seasonality7 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC7=c(AIC7,AIC(MOD7))
  
  MOD8=bglmer(y ~ Rainfall_Seasonality8 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC8=c(AIC8,AIC(MOD8))
  
  MOD9=bglmer(y ~ Rainfall_Seasonality9 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC9=c(AIC9,AIC(MOD9))
  
  MOD10=bglmer(y ~ Rainfall_Seasonality10 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC10=c(AIC10,AIC(MOD10))
  
  MOD11=bglmer(y ~ Rainfall_Seasonality11 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC11=c(AIC11,AIC(MOD11))
  
  MOD12=bglmer(y ~ Rainfall_Seasonality12 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC12=c(AIC12,AIC(MOD12))
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
mean(AIC10) #the best model: minimizing AIC
mean(AIC11)
mean(AIC12)


#here we select the best time window for NDVI seasonality (i.e. NDVI_S)
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
AIC12=c()

for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)
  
  MOD0=bglmer(y ~ NDVI_Seasonality0 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ NDVI_Seasonality1 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ NDVI_Seasonality2 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ NDVI_Seasonality3 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ NDVI_Seasonality4 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ NDVI_Seasonality5 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
  MOD6=bglmer(y ~ NDVI_Seasonality6 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC6=c(AIC6,AIC(MOD6))
  
  MOD7=bglmer(y ~ NDVI_Seasonality7 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC7=c(AIC7,AIC(MOD7))
  
  MOD8=bglmer(y ~ NDVI_Seasonality8 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC8=c(AIC8,AIC(MOD8))
  
  MOD9=bglmer(y ~ NDVI_Seasonality9 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC9=c(AIC9,AIC(MOD9))
  
  MOD10=bglmer(y ~ NDVI_Seasonality10 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC10=c(AIC10,AIC(MOD10))
  
  MOD11=bglmer(y ~ NDVI_Seasonality11 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC11=c(AIC11,AIC(MOD11))
  
  MOD12=bglmer(y ~ NDVI_Seasonality12 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC12=c(AIC12,AIC(MOD12))
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
mean(AIC3) 
mean(AIC4) #the best model: minimizing AIC
mean(AIC5)
mean(AIC6) 
mean(AIC7) 
mean(AIC8)
mean(AIC9) 
mean(AIC10) 
mean(AIC11)
mean(AIC12)




#here we select the best time window for non-seasonal variation in rainfall (i.e. Rain_NS)
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
AIC12=c()

for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)
  
  MOD0=bglmer(y ~ Dif_Rainfall0 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ Dif_Rainfall1 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ Dif_Rainfall2 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ Dif_Rainfall3 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ Dif_Rainfall4 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ Dif_Rainfall5 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
  MOD6=bglmer(y ~ Dif_Rainfall6 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC6=c(AIC6,AIC(MOD6))
  
  MOD7=bglmer(y ~ Dif_Rainfall7 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC7=c(AIC7,AIC(MOD7))
  
  MOD8=bglmer(y ~ Dif_Rainfall8 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC8=c(AIC8,AIC(MOD8))
  
  MOD9=bglmer(y ~ Dif_Rainfall9 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC9=c(AIC9,AIC(MOD9))
  
  MOD10=bglmer(y ~ Dif_Rainfall10 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC10=c(AIC10,AIC(MOD10))
  
  MOD11=bglmer(y ~ Dif_Rainfall11 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC11=c(AIC11,AIC(MOD11))
  
  MOD12=bglmer(y ~ Dif_Rainfall12 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC12=c(AIC12,AIC(MOD12))
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
mean(AIC3) 
mean(AIC4) #the best model: minimizing AIC
mean(AIC5)
mean(AIC6) 
mean(AIC7) 
mean(AIC8)
mean(AIC9) 
mean(AIC10) 
mean(AIC11)
mean(AIC12)


#here we select the best time window for non-seasonal variation in NDVI (i.e. NDVI_NS)
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
AIC12=c()

for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)
  
  MOD0=bglmer(y ~ Dif_NDVI0 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ Dif_NDVI1 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ Dif_NDVI2 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ Dif_NDVI3 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ Dif_NDVI4 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ Dif_NDVI5 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
  MOD6=bglmer(y ~ Dif_NDVI6 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC6=c(AIC6,AIC(MOD6))
  
  MOD7=bglmer(y ~ Dif_NDVI7 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC7=c(AIC7,AIC(MOD7))
  
  MOD8=bglmer(y ~ Dif_NDVI8 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC8=c(AIC8,AIC(MOD8))
  
  MOD9=bglmer(y ~ Dif_NDVI9 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC9=c(AIC9,AIC(MOD9))
  
  MOD10=bglmer(y ~ Dif_NDVI10 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC10=c(AIC10,AIC(MOD10))
  
  MOD11=bglmer(y ~ Dif_NDVI11 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC11=c(AIC11,AIC(MOD11))
  
  MOD12=bglmer(y ~ Dif_NDVI12 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC12=c(AIC12,AIC(MOD12))
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
mean(AIC3) #the best model: minimizing AIC
mean(AIC4) 
mean(AIC5)
mean(AIC6) 
mean(AIC7) 
mean(AIC8)
mean(AIC9) 
mean(AIC10) 
mean(AIC11)
mean(AIC12)


### here we chose between (i) the sine term (proxy of temperature, photoperiod, or any other seasonal variation)
# (ii) the NDVI terms (seasonal + non-seasonal), our proxy of food availability, 
# (iii) the rainfall terms (seasonal + non-seasonal)


AIC0=c()
AIC1=c()
AIC2=c()
j=ncol(TABFPPO)-1000
for (j in (ncol(TABFPPO)-1000):ncol(TABFPPO) ) {
  y=TABFPPO[,j]
  print(j)
  
  MOD0=bglmer(y ~ sin(Month_Rad+pi/3) + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ Rainfall_Seasonality10 + Dif_Rainfall4 + (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2)) 
  
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
#here we select Model 1 as it has the lowest AIC => NDVI values. 


#then we have to select the best group reproductive synchrony effect
#first metric: the mean number of conception over the past X months (X going from 0 to 6)

AIC0=c()
AIC1=c()
AIC2=c()
AIC3=c()
AIC4=c()
AIC5=c()
AIC6=c()


for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD0=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update0 + (Number_Conc_Update0:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update1 + (Number_Conc_Update1:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update2 + (Number_Conc_Update2:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
  MOD3=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update3 + (Number_Conc_Update3:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC3=c(AIC3,AIC(MOD3))
  
  MOD4=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update4 + (Number_Conc_Update4:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC4=c(AIC4,AIC(MOD4))
  
  MOD5=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update5 + (Number_Conc_Update5:Relative.rank) +
                (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC5=c(AIC5,AIC(MOD5))
  
  MOD6=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
               Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) +
               (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) +
               (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC6=c(AIC6,AIC(MOD6))
  
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
mean(AIC3) 
mean(AIC4) 
mean(AIC5)
mean(AIC6) #the best model: minimizing AIC


#our second metric on group reproductive synchrony: number of cycling females.

AIC0=c()
AIC1=c()
AIC2=c()

for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD0=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                NB_Cyc_Precis0 + (NB_Cyc_Precis0:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC0=c(AIC0,AIC(MOD0))
  
  MOD1=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                NB_Cyc_Precis1 + (NB_Cyc_Precis1:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC1=c(AIC1,AIC(MOD1))
  
  MOD2=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                NB_Cyc_Precis2 + (NB_Cyc_Precis2:Relative.rank) +
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  AIC2=c(AIC2,AIC(MOD2))
  
}

mean(AIC0) 
mean(AIC1) 
mean(AIC2)
#actually for these models (4.2), we did not use AIC as we did not have the same number of observations per model
#we checked the significance of the effect Number of cycling female in the 3 models, and showed the one with the srongest effects (O here)





### we ran full models with all interactions, but removed therafter the interactions between female rank and both seasonal
#and non-seaosnal NDVI variation. 

#Full MODEL 3.1 here : 
Estimate_NDVI_Seasonality=c()
Estimate_DifNDVI=c()
Estimate_Nb_Females=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Parity_primi=c()
Esimate_Compet=c()
Estimate_Interaction_Rank_Compet=c()

CILOW_NDVI_Seasonality=c()
CIUP_NDVI_Seasonality=c()
CILOW_DifNDVI=c()
CIUP_DifNDVI=c()
CILOW_Nb_Females=c()
CIUP_Nb_Females=c()
CILOW_TroopL=c()
CIUP_TroopL=c()
CILOW_TroopM=c()
CIUP_TroopM=c()
CILOW_Rank=c()
CIUP_Rank=c()
CILOW_Parity=c()
CIUP_Parity=c()
CILOW_Parity_primi=c()
CIUP_Parity_primi=c()
CILOW_Compet=c()
CIUP_Compet=c()
CILOW_Interaction_Rank_Compet=c()
CIUP_Interaction_Rank_Compet=c()

Chisq_NDVI_Seasonality=c()
Chisq_DifNDVI=c()
Chisq_Nb_Females=c()
Chisq_Troop=c()
Chisq_Rank=c()
Chisq_Parity=c()
Chisq_Compet=c()
Chisq_Interaction_Rank_Compet=c()

p_NDVI_Seasonality=c()
p_DifNDVI=c()
p_Nb_Females=c()
p_Troop=c()
p_Rank=c()
p_Parity=c()
p_Compet=c()
p_Interaction_Rank_Compet=c()


for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
               Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) + 
               (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  SS=summary(MOD)$coefficients
  CC=confint(MOD,method='Wald')
  AA=Anova(MOD)
  
  Estimate_NDVI_Seasonality=c(Estimate_NDVI_Seasonality,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_Nb_Females=c(Estimate_Nb_Females,SS[4,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[5,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[6,1])
  Estimate_Rank=c(Estimate_Rank,SS[7,1])
  Estimate_Parity_primi=c(Estimate_Parity_primi,SS[8,1])
  Esimate_Compet=c(Esimate_Compet,SS[9,1])
  Estimate_Interaction_Rank_Compet=c(Estimate_Interaction_Rank_Compet,SS[10,1])
  
  CILOW_NDVI_Seasonality=c(CILOW_NDVI_Seasonality,CC[3,1])
  CIUP_NDVI_Seasonality=c(CIUP_NDVI_Seasonality,CC[3,2])
  CILOW_DifNDVI=c(CILOW_DifNDVI,CC[4,1])
  CIUP_DifNDVI=c(CIUP_DifNDVI,CC[4,2])
  CILOW_Nb_Females=c(CILOW_Nb_Females,CC[5,1])
  CIUP_Nb_Females=c(CIUP_Nb_Females,CC[5,2])
  CILOW_TroopL=c(CILOW_TroopL,CC[6,1])
  CIUP_TroopL=c(CIUP_TroopL,CC[6,2])
  CILOW_TroopM=c(CILOW_TroopM,CC[7,1])
  CIUP_TroopM=c(CIUP_TroopM,CC[7,2])
  CILOW_Rank=c(CILOW_Rank,CC[8,1])
  CIUP_Rank=c(CIUP_Rank,CC[8,2])
  CILOW_Parity_primi=c(CILOW_Parity_primi,CC[9,1])
  CIUP_Parity_primi=c(CIUP_Parity_primi,CC[9,2])
  CILOW_Compet=c(CILOW_Compet,CC[10,1])
  CIUP_Compet=c(CIUP_Compet,CC[10,2])
  CILOW_Interaction_Rank_Compet=c(CILOW_Interaction_Rank_Compet,CC[11,1])
  CIUP_Interaction_Rank_Compet=c(CIUP_Interaction_Rank_Compet,CC[11,2])
  
  
  Chisq_NDVI_Seasonality=c(Chisq_NDVI_Seasonality,AA[1,1])
  Chisq_DifNDVI=c(Chisq_DifNDVI,AA[2,1])
  Chisq_Nb_Females=c(Chisq_Nb_Females,AA[3,1])
  Chisq_Troop=c(Chisq_Troop,AA[4,1])
  Chisq_Rank=c(Chisq_Rank,AA[5,1])
  Chisq_Parity=c(Chisq_Parity,AA[6,1])
  Chisq_Compet=c(Chisq_Compet,AA[7,1])
  Chisq_Interaction_Rank_Compet=c(Chisq_Interaction_Rank_Compet,AA[8,1])
  
  p_NDVI_Seasonality=c(p_NDVI_Seasonality,AA[1,3])
  p_DifNDVI=c(p_DifNDVI,AA[2,3])
  p_Nb_Females=c(p_Nb_Females,AA[3,3])
  p_Troop=c(p_Troop,AA[4,3])
  p_Rank=c(p_Rank,AA[5,3])
  p_Parity=c(p_Parity,AA[6,3])
  p_Compet=c(p_Compet,AA[7,3])
  p_Interaction_Rank_Compet=c(p_Interaction_Rank_Compet,AA[8,3])
  
}

# > summary(Estimate_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04124 0.11508 0.13538 0.13465 0.15551 0.23670 
# > summary(Esimate_Compet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.32257 -0.23151 -0.19985 -0.19951 -0.16864 -0.01749 
# > summary(Estimate_Interaction_Rank_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07687 0.19319 0.21774 0.21843 0.24528 0.32919 
# > 
#   > summary(Estimate_NDVI_Seasonality)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.29575 -0.17133 -0.13501 -0.13450 -0.10000  0.02616 
# > summary(Estimate_Nb_Females)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.210249 -0.065157 -0.027367 -0.028192  0.007478  0.152331 
# > summary(Estimate_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.08248 0.12392 0.13504 0.13561 0.14718 0.19694 
# > summary(Estimate_TroopM)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.35740 -0.72874 -0.59620 -0.61467 -0.49131 -0.05425 
# > summary(Estimate_Rank)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.01383  0.02771  0.03836  0.03840  0.04914  0.09314 
# > summary(Estimate_Parity_primi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3218 -0.2392 -0.2244 -0.2257 -0.2097 -0.1653 
# > 
#   > 
#   > 
#   > summary(CILOW_DifNDVI)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.130676 -0.051059 -0.029633 -0.029704 -0.007415  0.078873 
# > summary(CIUP_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2132  0.2805  0.2992  0.2990  0.3186  0.3945 
# > summary(CILOW_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.5379 -0.4387 -0.4065 -0.4058 -0.3741 -0.2128 
# > summary(CIUP_Compet)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.107217 -0.023491  0.006285  0.006781  0.036127  0.177786 
# > summary(CILOW_Interaction_Rank_Compet)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.105226  0.002932  0.028544  0.028491  0.055050  0.135197 
# > summary(CIUP_Interaction_Rank_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2590  0.3827  0.4079  0.4084  0.4358  0.5232 
# > 
#   > summary(CILOW_NDVI_Seasonality)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.4809 -0.3518 -0.3146 -0.3142 -0.2786 -0.1508 
# > summary(CIUP_NDVI_Seasonality)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.110628  0.009094  0.044605  0.045233  0.078968  0.203109 
# > summary(CILOW_Nb_Females)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.4982 -0.3518 -0.3134 -0.3143 -0.2787 -0.1357 
# > summary(CIUP_Nb_Females)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07774 0.22135 0.25836 0.25792 0.29338 0.44037 
# > summary(CILOW_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.2738 -0.2329 -0.2224 -0.2222 -0.2111 -0.1652 
# > summary(CIUP_TroopL)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4388  0.4809  0.4928  0.4934  0.5056  0.5591 
# > summary(CILOW_TroopM)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.790  -2.057  -1.921  -1.954  -1.815  -1.391 
# > summary(CIUP_TroopM)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.07505 0.61262 0.73040 0.72506 0.83694 1.28263 
# > summary(CILOW_Rank)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.18438 -0.14477 -0.13479 -0.13457 -0.12398 -0.08233 
# > summary(CIUP_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1567  0.2003  0.2115  0.2114  0.2222  0.2686 
# > summary(CILOW_Parity_primi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.7814 -0.6911 -0.6755 -0.6773 -0.6609 -0.6145 
# > summary(CIUP_Parity_primi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1378  0.2127  0.2269  0.2259  0.2415  0.2842 
# > 
#   > 
#   > summary(Chisq_DifNDVI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2211  1.8416  2.5877  2.7412  3.5056  8.6402 
# > summary(Chisq_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01248 1.31258 2.01120 2.15881 2.90975 6.63459 
# > summary(Chisq_Interaction_Rank_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6845  3.9605  5.0820  5.2324  6.3686 11.1264 
# > summary(Chisq_NDVI_Seasonality)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00006 1.20045 2.17174 2.46046 3.46497 9.80470 
# > summary(Chisq_Nb_Females)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.0148016 0.0760535 0.1762156 0.2349216 2.0474308 
# > summary(Chisq_Troop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4596  1.2476  1.5733  1.6553  2.0104  4.3663 
# > summary(Chisq_Rank)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000002 0.0218647 0.0639554 0.0929936 0.1362620 0.5241990 
# > summary(Chisq_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5196  0.8302  0.9500  0.9687  1.0776  1.8831 
# > 
#   > 
#   > 
#   > summary(p_DifNDVI)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003288 0.061163 0.107697 0.132145 0.174759 0.638208 
# > wilcox.test(p_DifNDVI, conf.int = TRUE, conf.level = 0.95)
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  p_DifNDVI
# V = 500500, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   0.1131514 0.1248814
# sample estimates:
#   (pseudo)median 
# 0.1189271 
# 
# > 
#   > summary(p_Compet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01000 0.08805 0.15614 0.19155 0.25193 0.91104 
# > wilcox.test(p_Compet, conf.int = TRUE, conf.level = 0.95)
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  p_Compet
# V = 500500, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   0.1648843 0.1810289
# sample estimates:
#   (pseudo)median 
# 0.1727599 
# 
# > 
#   > summary(p_Interaction_Rank_Compet)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008511 0.0116157 0.0241761 0.0362945 0.0465794 0.4080389 
# > wilcox.test(p_Interaction_Rank_Compet, conf.int = TRUE, conf.level = 0.95)
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  p_Interaction_Rank_Compet
# V = 500500, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   0.02718558 0.03080286
# sample estimates:
#   (pseudo)median 
# 0.02896238 
# 
# > 
#   > summary(p_NDVI_Seasonality)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001741 0.062682 0.140567 0.201855 0.273231 0.993827 
# > wilcox.test(p_NDVI_Seasonality, conf.int = TRUE, conf.level = 0.95)
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  p_NDVI_Seasonality
# V = 500500, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0
# 95 percent confidence interval:
#   0.1596545 0.1823030
# sample estimates:
#   (pseudo)median 
# 0.1707584 
# 
# > 
#   > summary(p_Nb_Females)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1525  0.6279  0.7827  0.7491  0.9032  0.9997 
# > summary(p_Troop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1127  0.3660  0.4554  0.4533  0.5359  0.7947 
# > summary(p_Rank)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4691  0.7120  0.8004  0.7947  0.8824  0.9997 
# > summary(p_Parity)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1700  0.2992  0.3297  0.3298  0.3622  0.4710




##now the full model 3.2: 

Estimate_NDVI_Seasonality=c()
Estimate_DifNDVI=c()
Estimate_Nb_Females=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Estimate_Parity_primi=c()
Esimate_Compet=c()

CILOW_NDVI_Seasonality=c()
CIUP_NDVI_Seasonality=c()
CILOW_DifNDVI=c()
CIUP_DifNDVI=c()
CILOW_Nb_Females=c()
CIUP_Nb_Females=c()
CILOW_TroopL=c()
CIUP_TroopL=c()
CILOW_TroopM=c()
CIUP_TroopM=c()
CILOW_Rank=c()
CIUP_Rank=c()
CILOW_Parity=c()
CIUP_Parity=c()
CILOW_Parity_primi=c()
CIUP_Parity_primi=c()
CILOW_Compet=c()
CIUP_Compet=c()


Chisq_NDVI_Seasonality=c()
Chisq_DifNDVI=c()
Chisq_Nb_Females=c()
Chisq_Troop=c()
Chisq_Rank=c()
Chisq_Parity=c()
Chisq_Compet=c()

p_NDVI_Seasonality=c()
p_DifNDVI=c()
p_Nb_Females=c()
p_Troop=c()
p_Rank=c()
p_Parity=c()
p_Compet=c()


for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
               NB_Cyc_Precis0 + 
               (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  SS=summary(MOD)$coefficients
  CC=confint(MOD,method='Wald')
  AA=Anova(MOD)
  
  Estimate_NDVI_Seasonality=c(Estimate_NDVI_Seasonality,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_Nb_Females=c(Estimate_Nb_Females,SS[4,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[5,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[6,1])
  Estimate_Rank=c(Estimate_Rank,SS[7,1])
  Estimate_Parity_primi=c(Estimate_Parity_primi,SS[8,1])
  Esimate_Compet=c(Esimate_Compet,SS[9,1])
  
  
  CILOW_NDVI_Seasonality=c(CILOW_NDVI_Seasonality,CC[3,1])
  CIUP_NDVI_Seasonality=c(CIUP_NDVI_Seasonality,CC[3,2])
  CILOW_DifNDVI=c(CILOW_DifNDVI,CC[4,1])
  CIUP_DifNDVI=c(CIUP_DifNDVI,CC[4,2])
  CILOW_Nb_Females=c(CILOW_Nb_Females,CC[5,1])
  CIUP_Nb_Females=c(CIUP_Nb_Females,CC[5,2])
  CILOW_TroopL=c(CILOW_TroopL,CC[6,1])
  CIUP_TroopL=c(CIUP_TroopL,CC[6,2])
  CILOW_TroopM=c(CILOW_TroopM,CC[7,1])
  CIUP_TroopM=c(CIUP_TroopM,CC[7,2])
  CILOW_Rank=c(CILOW_Rank,CC[8,1])
  CIUP_Rank=c(CIUP_Rank,CC[8,2])
  CILOW_Parity_primi=c(CILOW_Parity_primi,CC[9,1])
  CIUP_Parity_primi=c(CIUP_Parity_primi,CC[9,2])
  CILOW_Compet=c(CILOW_Compet,CC[10,1])
  CIUP_Compet=c(CIUP_Compet,CC[10,2])
  
  Chisq_NDVI_Seasonality=c(Chisq_NDVI_Seasonality,AA[1,1])
  Chisq_DifNDVI=c(Chisq_DifNDVI,AA[2,1])
  Chisq_Nb_Females=c(Chisq_Nb_Females,AA[3,1])
  Chisq_Troop=c(Chisq_Troop,AA[4,1])
  Chisq_Rank=c(Chisq_Rank,AA[5,1])
  Chisq_Parity=c(Chisq_Parity,AA[6,1])
  Chisq_Compet=c(Chisq_Compet,AA[7,1])
  
  
  p_NDVI_Seasonality=c(p_NDVI_Seasonality,AA[1,3])
  p_DifNDVI=c(p_DifNDVI,AA[2,3])
  p_Nb_Females=c(p_Nb_Females,AA[3,3])
  p_Troop=c(p_Troop,AA[4,3])
  p_Rank=c(p_Rank,AA[5,3])
  p_Parity=c(p_Parity,AA[6,3])
  p_Compet=c(p_Compet,AA[7,3])
  
  
}

# > mean(Estimate_NDVI_Seasonality) #
# [1] -0.3020513
# > mean(Estimate_DifNDVI) # 
# [1] 0.2510957
# > mean(Estimate_Nb_Females) #
# [1] -0.1031005
# > mean(Estimate_TroopL) #
# [1] 0.106887
# > mean(Estimate_TroopM) #
# [1] -6.711624
# > mean(Estimate_Rank) #
# [1] -0.08883971
# > mean(Estimate_Parity_primi) #
# [1] -0.1226803
# > mean(Esimate_Compet) #
# [1] 0.0470487
# > 
#   > mean(CILOW_NDVI_Seasonality) #
# [1] -0.5578969
# > mean(CIUP_NDVI_Seasonality) #
# [1] -0.04620573
# > mean(CILOW_DifNDVI) #
# [1] 0.01174567
# > mean(CIUP_DifNDVI) #
# [1] 0.4904457
# > mean(CILOW_Nb_Females) #
# [1] -0.4545546
# > mean(CIUP_Nb_Females) #
# [1] 0.2483537
# > mean(CILOW_TroopL) #
# [1] -0.4437014
# > mean(CIUP_TroopL) #
# [1] 0.6574754
# > mean(CILOW_TroopM) #
# [1] -949.9452
# > mean(CIUP_TroopM) #
# [1] 936.522
# > mean(CILOW_Rank) #
# [1] -0.3459627
# > mean(CIUP_Rank) #
# [1] 0.1682833
# > mean(CILOW_Parity_primi) #
# [1] -0.8223633
# > mean(CIUP_Parity_primi) #
# [1] 0.5770028
# > mean(CILOW_Compet) #
# [1] -0.2295269
# > mean(CIUP_Compet) #
# [1] 0.3236243
# > 
#   > 
#   > mean(Chisq_NDVI_Seasonality) #
# [1] 5.493864
# > mean(Chisq_DifNDVI) #
# [1] 4.502669
# > mean(Chisq_Nb_Females) #
# [1] 0.4134709
# > mean(Chisq_Troop) #
# [1] 0.589187
# > mean(Chisq_Rank) #
# [1] 0.5787618
# > mean(Chisq_Parity) #
# [1] 0.2547517
# > mean(Chisq_Compet) #
# [1] 0.2861497
# > 
#   > 
#   > mean(p_NDVI_Seasonality) #
# [1] 0.03097167
# > mean(p_DifNDVI) #
# [1] 0.0641455
# > mean(p_Nb_Females) #
# [1] 0.5756969
# > mean(p_Troop) #
# [1] 0.7628713
# > mean(p_Rank) #0
# [1] 0.5173142
# > mean(p_Parity) #
# [1] 0.6925928
# > mean(p_Compet) #
# [1] 0.6851783




## lastly, we run model diagnostics for both models 3.1 and 3.2, as follow:

#MODEL 3.1
MOD=bglmer(FPPO ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) + 
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')

r.squaredGLMM(MOD)
vif(MOD) #only the 2 fixed effects of female age have vif >2.5, which was expected. ==> no collinearity
simulationOutput <- simulateResiduals(fittedModel = MOD, plot = F)
plot(simulationOutput)
hist(simulationOutput)


#MODEL 3.2
MOD=bglmer(FPPO ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             NB_Cyc_Precis0 + 
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')

r.squaredGLMM(MOD)
vif(MOD) #only the 2 fixed effects of female age have vif >2.5, which was expected. ==> no collinearity
simulationOutput <- simulateResiduals(fittedModel = MOD, plot = F)
plot(simulationOutput)
hist(simulationOutput)






#we ran similar models, following the same procedure,
# but with different outcomes, for conception probability models (MODELS 4.1 and 4.2)




