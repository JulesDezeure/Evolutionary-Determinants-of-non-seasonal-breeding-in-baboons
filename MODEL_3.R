library(lubridate)
library(lme4)
library(blme)
library(MuMIn)
library(DHARMa)
library(car)


#The first step is to take into account the uncertainty around cycle resumption dates
#to do so, we ran 1000 randomized dates, drawn in the possible dates

TABFPPO=read.csv2('file_path of the excel tab Cycle resumption')


######################################################################################################

# This section of code is just to show how we obtained our randomized date of reproductive event
# drawing 1000 dates in the pool of possible date given their uncertainty. 

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

######################################################################################################



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


#we chose which of the group size effect (simple or quadratic) is the best now:

TABFPPO$NB_Adult_Females_V2carre=TABFPPO$NB_Adult_Females_V2*TABFPPO$NB_Adult_Females_V2 # quadratic effect of group size
TABFPPO$NB_Adult_Females_V2carre=scale(TABFPPO$NB_Adult_Females_V2carre)
TABFPPO$NB_Adult_Females_V2=scale(TABFPPO$NB_Adult_Females_V2)
#Also we scale all the quantitative fixed effects:
#female rank
TABFPPO$Relative.rank=scale(TABFPPO$Relative.rank)


AIC1=c()
AIC2=c()
for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD1=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + 
                Troop + Relative.rank + Parity + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD2=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + 
                NB_Adult_Females_V2carre + Troop + Relative.rank + Parity + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  
  
  AIC1=c(AIC1,AIC(MOD1))
  AIC2=c(AIC2,AIC(MOD2))
  
}


mean(AIC1)
mean(AIC2)
#we chose a simple effect of group size as mean(AIC1) < mean(AIC2)


#########  then we have to select the best group reproductive synchrony effect


#metrics of group reproductive synchrony to scale
TABFPPO$Number_Conc_Update1=scale(TABFPPO$Number_Conc_Update1)
TABFPPO$Number_Conc_Update2=scale(TABFPPO$Number_Conc_Update2)
TABFPPO$Number_Conc_Update4=scale(TABFPPO$Number_Conc_Update4)
TABFPPO$Number_Conc_Update6=scale(TABFPPO$Number_Conc_Update6)
TABFPPO$Number_Conc_Update0=scale(TABFPPO$Number_Conc_Update0)
TABFPPO$SynchConc1=scale(TABFPPO$SynchConc1)
TABFPPO$SynchConc2=scale(TABFPPO$SynchConc2)
TABFPPO$SynchConc3=scale(TABFPPO$SynchConc3)
TABFPPO$AfterConc1=scale(TABFPPO$AfterConc1)
TABFPPO$AfterConc2=scale(TABFPPO$AfterConc2)
TABFPPO$AfterConc4=scale(TABFPPO$AfterConc4)
TABFPPO$AfterConc6=scale(TABFPPO$AfterConc6)


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


for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  
  print(j)
  y=TABFPPO[,j]
  
  MOD0=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update0 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD1=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update1 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD2=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                AfterConc1 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD3=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                SynchConc1 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD4=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update2 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD5=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                AfterConc2 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD6=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                SynchConc2 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD7=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update4 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD8=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                AfterConc4 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD9=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                SynchConc3 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD10=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                Number_Conc_Update6 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  MOD11=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
                AfterConc6 + 
                (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  

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
# minimum AIC for Model 10, i.e. we will include in the full model the reprodutive synchrony metric : 
# mean  number of conceptions over the past 6 months (i.e. Number_Conc_Update6)



### we ran full models with all interactions, to check which of the rank interactions are significant 
# and should be kept in our final model. 

Estimate_Interaction_Rank_Synchro=c()
Estimate_Interaction_Rank_NDVIS=c()
Estimate_Interaction_Rank_DifNDVI=c()
VAR_Interaction_Rank_Synchro=c()
VAR_Interaction_Rank_NDVIS=c()
VAR_Interaction_Rank_DifNDVI=c()

for (j in (ncol(TABFPPO)-999):ncol(TABFPPO)) {
  print(j)
  y=TABFPPO[,j]
  
  MOD=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
               Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) + 
               (NDVI_Seasonality4:Relative.rank) + (Dif_NDVI3:Relative.rank) + 
               (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
  SS=summary(MOD)$coefficients
  VV=vcov(MOD)
  Estimate_Interaction_Rank_Synchro=c(Estimate_Interaction_Rank_Synchro,SS[10,1])
  Estimate_Interaction_Rank_NDVIS=c(Estimate_Interaction_Rank_NDVIS,SS[11,1])
  Estimate_Interaction_Rank_DifNDVI=c(Estimate_Interaction_Rank_DifNDVI,SS[12,1])
  VAR_Interaction_Rank_Synchro=c(VAR_Interaction_Rank_Synchro,VV[10, 10])
  VAR_Interaction_Rank_NDVIS=c(VAR_Interaction_Rank_NDVIS,VV[11, 11])
  VAR_Interaction_Rank_DifNDVI=c(VAR_Interaction_Rank_DifNDVI,VV[12, 12])
  
}
m=1000 
n=1000
k=length(coef(MOD))

# Interaction rank : synchrony (significant, kept in our final model)
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

# Interaction rank : seasonal NDVI (not significant, not kept in our final model)
MAT=data.frame(Estimate_Interaction_Rank_NDVIS,VAR_Interaction_Rank_NDVIS)  
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

# Interaction Rank : non-seasonal NDVI (not significant, not kept in final model)
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




# Here is the full model, with results presented in Table 4 of the main text
# we controlled cycle resumption date uncertainty with our 1000 randomizations
# and applied Rubin's rule to extract pooled estimates, SD and Pvalues. 

Estimate_NDVI_Seasonality=c()
Estimate_DifNDVI=c()
Estimate_Nb_Females=c()
Estimate_TroopL=c()
Estimate_TroopM=c()
Estimate_Rank=c()
Estimate_Parity=c()
Esimate_Compet=c()
Estimate_Interaction_Rank_Compet=c()

VAR_NDVI_Seasonality=c()
VAR_DifNDVI=c()
VAR_Nb_Females=c()
VAR_TroopL=c()
VAR_TroopM=c()
VAR_Rank=c()
VAR_Parity=c()
VAR_Compet=c()
VAR_Interaction_Rank_Compet=c()

p_NDVI_Seasonality=c()
p_DifNDVI=c()
p_Nb_Females=c()
p_TroopL=c()
p_TroopM=c()
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
  VV=vcov(MOD)
  
  Estimate_NDVI_Seasonality=c(Estimate_NDVI_Seasonality,SS[2,1])
  Estimate_DifNDVI=c(Estimate_DifNDVI,SS[3,1])
  Estimate_Nb_Females=c(Estimate_Nb_Females,SS[4,1])
  Estimate_TroopL=c(Estimate_TroopL,SS[5,1])
  Estimate_TroopM=c(Estimate_TroopM,SS[6,1])
  Estimate_Rank=c(Estimate_Rank,SS[7,1])
  Estimate_Parity=c(Estimate_Parity,SS[8,1])
  Esimate_Compet=c(Esimate_Compet,SS[9,1])
  Estimate_Interaction_Rank_Compet=c(Estimate_Interaction_Rank_Compet,SS[10,1])
  
  VAR_NDVI_Seasonality=c(VAR_NDVI_Seasonality,VV[2, 2])
  VAR_DifNDVI=c(VAR_DifNDVI,VV[3, 3])
  VAR_Nb_Females=c(VAR_Nb_Females,VV[4, 4])
  VAR_TroopL=c(VAR_TroopL,VV[5, 5])
  VAR_TroopM=c(VAR_TroopM,VV[6, 6])
  VAR_Rank=c(VAR_Rank,VV[7, 7])
  VAR_Parity=c(VAR_Parity,VV[8, 8])
  VAR_Compet=c(VAR_Compet,VV[9, 9])
  VAR_Interaction_Rank_Compet=c(VAR_Interaction_Rank_Compet,VV[10, 10])
  
  p_NDVI_Seasonality=c(p_NDVI_Seasonality,SS[2,4])
  p_DifNDVI=c(p_DifNDVI,SS[3,4])
  p_Nb_Females=c(p_Nb_Females,SS[4,4])
  p_TroopL=c(p_TroopL,SS[5,4])
  p_TroopM=c(p_TroopM,SS[6,4])
  p_Rank=c(p_Rank,SS[7,4])
  p_Parity=c(p_Parity,SS[8,4])
  p_Compet=c(p_Compet,SS[9,4])
  p_Interaction_Rank_Compet=c(p_Interaction_Rank_Compet,SS[10,4])
  
}

m=1000 
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


# Interaction synchrony : rank
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



# Number of adult females in the group
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



# GROUP (L)
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


# GROUP (M)
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
-(betweenVar <- mean(MAT[,2])) # mean of variances
(withinVar <- sd(MAT[,1])^2) # variance of variances
(dfCorrection <- (nrow(MAT)+1)/(nrow(MAT))) # dfCorrection
(totVar <- betweenVar+ withinVar*dfCorrection) # total variance
(pooledSE <- sqrt(totVar)) 
-lambda <- (withinVar + (withinVar/m))/totVar
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

# Female parity (primiparous)
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
tableau



## lastly, we run model diagnostics for both models 3.1 and 3.2, as follow:

MOD=bglmer(FPPO ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) + 
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')

r.squaredGLMM(MOD)
vif(MOD) #only the 2 fixed effects of female age have vif >2.5, which was expected. ==> no collinearity
simulationOutput <- simulateResiduals(fittedModel = MOD, plot = F)
plot(simulationOutput)
hist(simulationOutput)



#we ran similar models, following the same procedure,
# but with different outcomes, for conception probability models (MODEL 4)




