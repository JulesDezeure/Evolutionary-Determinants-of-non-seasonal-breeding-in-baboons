#This is code to simulate the effect of group synchrony on birth seasonality (Appendix E)
#Code written by Lugdiwine Burtschell
#Code developed by Lugdiwine Burtschell, Jules Dezeure and Elise Huchard

library(CircStats)
library(ggplot2)



input_file1<-read.csv2("Extract_ModelIBI_Rank_Synchro_Interaction.csv")
input_file2<-read.csv2("Births_Distribution.csv")

######### Parameters for the simulation ######### 

nbF<-18 #number of female in the group (multiple of 3)
rank<-c(rep("high",nbF/3),rep("mid",nbF/3),rep("low",nbF/3)) #ranks associated to each female
tmax<-15*360 #length of simulation, in days (1 month = 30 days / 1 year = 12 months = 360 days)
IBI<-c(654.5,691,696.3)#mean of IBI for high, mid and low-rank females

#mean birthDate in our population
meanBD<-322 #November 18

#variation of IBI with group synchrony
# probaSynchrony<-read.xlsx(file=input_file1,1)
probaSynchrony<-input_file1
interceptRank<-estimateRank<-vector()
for (i in 1:3){
  probaSynchronyRank<-subset(probaSynchrony,probaSynchrony$group==i)
  model<-lm(probaSynchronyRank$y~probaSynchronyRank$x)
interceptRank[i]<-model$coefficients[1]
  estimateRank[i]<-model$coefficients[2]
}
f_high<-function(x){
  return(interceptRank[1]+estimateRank[1]*x)}
f_mid<-function(x){
  return(interceptRank[2]+estimateRank[2]*x)}
f_low<-function(x){
  return(interceptRank[3]+estimateRank[3]*x)}


#Observed birth(for initiation)
#birthsDistribution<-read.xlsx(file=input_file2,1)
birthsDistribution<-input_file2


birthsObserved<-vector()
for (i in 1:12){
  birthsObserved<-c(birthsObserved,rep(i,birthsDistribution$Count[i]))
}



######### Simulation function ######### 

simulation<-function(random,s){#random= True for random model and False for synchrony model 
  #s =degree of birth seasonality (0= no seasonality, 1 = complete seasonality)
  
  #Initiation: each female gives birth within the first 2 years 
  lastBirth<-sample(birthsObserved,nbF)+ #birth month randomly picked from all observed values
    sample(c(0,1),nbF,replace=T)*12  #and randomly assigned to year 1 or 2
  lastBirth<-lastBirth*30-15#last birth date for each female, in days, at day 15 of each month
  
  #addition of seasonality
  BD<-lastBirth%%360 #day of year
  BD[BD<=(meanBD-180)]<-BD[BD<=(meanBD-180)]+360 #birth date centered on meanBD
  deviation<-meanBD-BD #positive if Birth should occur in the 180 days before meanBD and negative if after
  lastBirth<-lastBirth+round(s*deviation) #birth becomes closer to mean BirthDate (births can become negative)
  
  births<-lastBirth #record of all births in the simulation
  lastIBI<-rep(tmax,nbF) #last IBI computed for each female
  
  for (t in min(lastBirth):tmax){
    for (female in 1:nbF){
      if((t-lastBirth[female])>=round(lastIBI[female])){#if the female reaches the end of IBI, it gives birth again
        lastBirth[female]<-t
        births<-c(births,t)
      }
      
      if (t==lastBirth[female]){ #if an infant is born, we compute the associated IBI
        if(random){ #random model
          if (rank[female]=="high"){
            lastIBI[female]<-IBI[1]
          }
          if (rank[female]=="mid"){
            lastIBI[female]<-IBI[2]
          }
          if (rank[female]=="low"){
            lastIBI[female]<-IBI[3]
          }
        }
        else{#synchrnoy model
          Birth3mois<-length(births[births>=(t-90)&births<=t]) #Group synchrony: number of birth in the past 3 months
          
          if (rank[female]=="high"){
            lastIBI[female]<-f_high(Birth3mois)
          }
          if (rank[female]=="mid"){
            lastIBI[female]<-f_mid(Birth3mois)
          }
          if (rank[female]=="low"){
            lastIBI[female]<-f_low(Birth3mois)
          }
        }
        #addition of seasonality
        nextBD<-(t+lastIBI[female])%%360 #day of year
        if(nextBD<=(meanBD-180)){nextBD=nextBD+360} #next birth date centered on meanBD
        deviation<-meanBD-nextBD #positive if next Birth should occur in the 180 days before meanBD and negative if after
        lastIBI[female]<-lastIBI[female]+(s*deviation) #next birth becomes closer to mean BirthDate
        
      } 
    }
  }
  
  birth<-births%%360 #day of year
  birth[birth==0]<-360
  rad_birth<-birth*2*pi/360 #birth in radians
  
  return (as.numeric(r.test(rad_birth)[1])) #Rayleigh statistics
}


######### Results of the simualtion ######### 


#Code for the Figure S3 here :
nruns<-1000

r_random<-vector()
r_synchrony<-vector()

for (i in 1:nruns){
  print(i)
  r_random<-c(r_random,simulation(T,0))
  r_synchrony<-c(r_synchrony,simulation(F,0))
}

t.test(r_random,r_synchrony)

Random=rep('Random',1000)
Asynchrony=rep('Asynchrony',1000)
NEW=data.frame(r_random,Random)
NEW2=data.frame(r_synchrony,Asynchrony)
names(NEW)[1]='R'
names(NEW)[2]='Type'
names(NEW2)[1]='R'
names(NEW2)[2]='Type'
NEW=rbind(NEW,NEW2)
library(ggplot2)
NEW$Type=as.factor(NEW$Type)
NEW$Type=factor(NEW$Type,c('Random','Asynchrony'))
summary(NEW$Type)


ggplot(NEW,aes(x=Type,y=R,fill=Type)) + geom_boxplot() + 
  scale_fill_manual(values=c('#56B4E9','#E69F00')) + 
  scale_x_discrete(labels=c('absent','present')) +
  labs(y='Strength of reproductive seasonality (R)',x='Group reproductive synchrony effect') + 
  theme(legend.position='none',
        axis.title = element_text(size=20),axis.text=element_text(size=18))

# Code for the Figure S4 here:
nruns<-1000
df<-data.frame()
#seq(0.11,1,by=0.01)
for (s in seq(0.11,1,by=0.01)){
  print(s)
  r_randomS<-vector()
  r_synchronyS<-vector()
  for (i in 1:nruns){
    r_randomS<-c(r_randomS,simulation(T,s))
    r_synchronyS<-c(r_synchronyS,simulation(F,s))
  }
  df<-rbind(df,c(s,mean(r_randomS),sd(r_randomS),mean(r_synchronyS),sd(r_synchronyS)))
}
colnames(df)<-c("s","r_random","sd_random","r_synchrony","sd_synchrony")


PLOT=ggplot(df) + geom_line(aes(x=s,y=r_random, color='r_random'), size=1.5) + 
  geom_ribbon(aes(x=s,ymin=r_random - sd_random,ymax=r_random + sd_random, fill='r_random'), 
              alpha=0.3) +
  geom_line(aes(x=s,y=r_synchrony, color='r_synchrony'), size=1.5) + 
  geom_ribbon(aes(x=s,ymin=r_synchrony - sd_synchrony,ymax=r_synchrony + sd_synchrony,fill='r_synchrony'), 
            alpha=0.3)
  
PLOT=PLOT +   scale_color_manual(values = c(
  'r_random' = 'blue','r_synchrony' = 'red'), 
  labels=c('Random models','Synchrony models')) + 
  scale_fill_manual(values = c(
    'r_random' = 'blue','r_synchrony' = 'red'), 
    labels=c('Random models','Synchrony models')) +
  labs(y='Strength of reproductive seasonality (R)',x='Level of reproductive seasonality (s)') +
  theme(axis.text.x = element_text(size=20), axis.title.x =element_text(size=20),
        axis.text.y = element_text(size=20), axis.title.y =element_text(size=20), 
        legend.title = element_blank(), legend.text=element_text(size=20), 
        legend.justification = 'center',legend.position=c(0.8,0.4), 
        legend.background = element_rect(colour = 'transparent', fill='transparent')) + 
  guides(fill=FALSE)
PLOT


  