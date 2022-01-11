library(lme4)
library(ggplot2)
library(car)
library(MuMIn)
library(blme)
library(cowplot)

############################################################################################################
#######                     FIGURE 1                ########################################################
############################################################################################################
tab_NDVI=read.csv2('tab_NDVI.csv')

#we compute the mean NDVI
tab_NDVI$K_NDVI_J=mean(tab_NDVI$mean_month_NDVI_J)
#we compute NDVI_NS here
tab_NDVI$NDVI_Imprev_J=tab_NDVI$mean_month_NDVI_J-tab_NDVI$K_NDVI_J-tab_NDVI$NDVI_APrev_J
#we only keep NDVI from 2004 (start of the study - 1 year)
tab_NDVI=tab_NDVI[which(tab_NDVI$year>2003),]

# PANEL A of Figure 1

B=ggplot(tab_NDVI) + 
  geom_line(aes(x=Time,y=mean_month_NDVI_J, color='mean_month_NDVI_J'),size=1.5, alpha=0.6) +
  geom_line(aes(x=Time,y=K_NDVI_J),size=0.8, color='black', linetype="longdash", alpha=0.5) + 
  geom_line(aes(x=Time,y=NDVI_APrev_J+K_NDVI_J, color='NDVI_APrev_J'),size=0.8, linetype="longdash") +
  geom_line(aes(x=Time,y=NDVI_Imprev_J+K_NDVI_J, color='NDVI_Imprev_J'),size=0.8, linetype="longdash")
B
for (i in (1:16)) {    
  B = B +geom_vline(xintercept =seq(from=13,to=195, by=12),linetype="dashed",color="grey",size=0.2)
}
B
B = B + scale_x_continuous(breaks=c(31,91,151), 
                           labels=c(2005,2010,2015)) + 
  scale_color_manual(values = c(
    'mean_month_NDVI_J' = 'green','NDVI_APrev_J' = 'blue', 'NDVI_Imprev_J' = 'red'), 
    labels=c('NDVI','NDVI_S','NDVI_NS')) 

B = B +  labs(y='Monthly NDVI',x='Year') + 
  theme(axis.text.x = element_text(size=20), axis.title.x =element_text(size=20),
        axis.text.y = element_text(size=20), axis.title.y =element_text(size=20), 
        legend.title = element_blank(), legend.text=element_text(size=20), 
        legend.justification = 'center',legend.position=c(0.8,0.8), 
        legend.background = element_rect(colour = 'transparent', fill='transparent')) 
B


# now code for the panel B of Figure 1
TAB=read.csv2("tab mortality.csv") 
TAB=TAB[which(TAB$Death_18months=="1"|TAB$Death_18months=="0"), ] #195 donnees
TAB$Death_18months=as.factor(TAB$Death_18months)
TAB$Death_18months=factor(TAB$Death_18months)
TAB$Mother=factor(TAB$Mother)
TAB$Relative_Rank=scale(TAB$Relative_Rank, center=TRUE, scale=TRUE)
TAB$Dif_NDVI=scale(TAB$Dif_NDVI)
TAB$Nb_Birth_After1=scale(TAB$Nb_Birth_After1)
TAB$NB_Adult_Females_V2=scale(TAB$NB_Adult_Females_V2)
MOD=glmer(Death_18months ~ sin(Birth_radian+pi/2) + Dif_NDVI + Troop + Relative_Rank + Parity + 
            Sex + NB_Adult_Females_V2 + Nb_Birth_After1 + (1|Mother) , 
          glmerControl(optimizer = 'bobyqa'), data=TAB, family='binomial')
FITTED=fitted(MOD)

TAB=read.csv2("tab mortality", stringsAsFactors = TRUE) 
TAB=TAB[which(TAB$Death_18months=="1"|TAB$Death_18months=="0"), ] #195 donnees
TAB$Death_18months=as.factor(TAB$Death_18months)
TAB$Death_18months=factor(TAB$Death_18months)
TAB=TAB[which(is.na(TAB$Sex)==FALSE),]
TAB=cbind(TAB,FITTED)
C=ggplot(TAB,aes(x=Dif_NDVI,y=FITTED)) + geom_point() + 
  geom_smooth(method='glm',size=2,col="#0072B2",method.args = list(family=binomial)) + 
  labs(x='Non-seasonal NDVI (NDVI_NS)',y='Infant probability to die before weaning') +  
  theme(axis.title.x = element_text(size=20),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),axis.text.y = element_text(size=15)) + 
  scale_x_continuous(breaks=c(-0.05,0,0.05))
C

PLOT_FIN=plot_grid(B,C,nrow=2,ncol=1, labels='AUTO',label_size = 15)
PLOT_FIN
ggsave("FIG1.tiff", units="in", width=8, height=12, dpi=300, compression = 'jpeg')



############################################################################################################
#######                     FIGURE 2                ########################################################
############################################################################################################

### The figure 2 panel A has been created manually on another software. Here you have the code to create panels
# B and C. 

TAB=read.csv2('ForFigure2.csv', stringsAsFactors = TRUE)
TAB$When=factor(TAB$When,levels=c('before','around','after'))
CONC=ggplot(TAB,aes(y=AIC_conception,x=Window,fill=When)) + 
  geom_bar(stat='identity', position=position_dodge(width=0.8),width=0.8) +
  geom_text(aes(label=WindowPrecis), vjust=-0.95, color="black",
            position = position_dodge(0.8), size=10, 
            fontface=ifelse(TAB$WindowPrecis=='[4a]',2,1)) + 
  scale_y_continuous(name='ΔAIC of conception probability models',limits = c(0,18)) +
  scale_fill_manual(values=c('#0072B2','#D55E00','#009E73')) +
  theme(axis.title.y = element_text(size=25, color='black'),axis.text.x = element_blank(),
        axis.text.y = element_text(size=20, color='black'), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = 'none', 
        axis.ticks.length = unit(0.3,'cm'),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
CONC

FPPO=ggplot(TAB,aes(y=AIC_FPPO,x=Window,fill=When)) + 
  geom_bar(stat='identity', position=position_dodge(width=0.8),width=0.8) +
  geom_text(aes(label=WindowPrecis), vjust=-0.95, color="black",
            position = position_dodge(0.8), size=10, 
            fontface=ifelse(TAB$WindowPrecis=='[6a]',2,1)) + #stylé c'est pour mettre en gras uniquement le meilleur modèle
  scale_y_continuous(name='ΔAIC of cycle resumptions probability models', limits=c(0,10), breaks=c(0,5,10)) +
  scale_fill_manual(values=c('#0072B2','#D55E00','#009E73')) +
  theme(axis.title.y = element_text(size=25, color='black'),axis.text.x = element_blank(),
        axis.text.y = element_text(size=20, color='black'), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = 'none',
        axis.ticks.length = unit(0.3,'cm'),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

FPPO
plot_grid(FPPO, CONC ,nrow=2,ncol=1, labels = c('B','C'), label_size = 25, label_colour = 'black')
ggsave("FIG2 right side.tiff", units="in", width=8, height=16, dpi=300, compression = 'jpeg')


############################################################################################################
#######                     FIGURE 3                ########################################################
############################################################################################################


# Panel A (IBI)
TAB=read.csv2('IBI tab.csv') 
TAB$DeathInBw=as.factor(TAB$DeathInBw)
TAB=TAB[which(TAB$DeathInBw=='0'),]
TAB=TAB[which(TAB$DeathFirstYear1=="n"),]
TAB=TAB[which(TAB$Borndead=='n'),] 
TAB$Mother=factor(TAB$Mother)
TAB$Relative_Rank=scale(TAB$Relative_Rank)
TAB$Dif_NDVI=scale(TAB$Dif_NDVI)
TAB$NB_Births_2months_Around=scale(TAB$NB_Births_2months_Around)
TAB$Age_Year_Mother_CARRE=TAB$Age_Year_Mother*TAB$Age_Year_Mother
TAB$Age_Year_Mother=scale(TAB$Age_Year_Mother)
TAB$Age_Year_Mother_CARRE=scale(TAB$Age_Year_Mother_CARRE)
TAB$NB_Adult_Females_V2carre=TAB$NB_Adult_Females_V2*TAB$NB_Adult_Females_V2
TAB$NB_Adult_Females_V2carre=scale(TAB$NB_Adult_Females_V2carre)
TAB$NB_Adult_Females_V2=scale(TAB$NB_Adult_Females_V2)


MOD=lmer(IBIDays ~ sin(Birth_Rad+pi/6) + Dif_NDVI + NB_Adult_Females_V2 + NB_Adult_Females_V2carre + 
           Troop + Relative_Rank + Parity + Sex1 + 
           Age_Year_Mother + Age_Year_Mother_CARRE +
           NB_Births_2months_Around + (NB_Births_2months_Around:Relative_Rank)  + (1|Mother),data=TAB)
FITTED=fitted(MOD) 

TAB=read.csv2('IBI tab.csv') 

TAB$DeathInBw=as.factor(TAB$DeathInBw)
TAB=TAB[which(TAB$DeathInBw=='0'),] #on enleve ces 10 IBI
TAB=TAB[which(TAB$DeathFirstYear1=="n"),]
TAB=TAB[which(TAB$Borndead=='n'),] 

TAB=cbind(TAB,FITTED)
TAB$Rank=factor(TAB$Rank,c('High-Ranking','Mid-Ranking','Low-Ranking'))
IBI=ggplot(TAB, aes(x=Nb_Birth_Before3, y=FITTED,color=Rank)) + geom_point() +
  geom_smooth(method="lm", se=FALSE, size=1.8) + 
  labs(x='Group reproductive synchrony',y='IBI (days)') +  
  scale_colour_colorblind() +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=15),
        legend.position="top",legend.text = element_text(size=15),legend.title=element_blank())

IBI



#Cycle resumption (Panel B)

TAB=read.csv2('Cycle resumption tab.csv')
TAB$ID=as.factor(TAB$ID)
TAB$FPPO=as.factor(TAB$FPPO)
TABFPPO=TAB
TABFPPO$Relative.rank=scale(TABFPPO$Relative.rank, center=TRUE, scale=TRUE)
TABFPPO$NDVI_Seasonality4=scale(TABFPPO$NDVI_Seasonality4)
TABFPPO$Dif_NDVI3=scale(TABFPPO$Dif_NDVI3)
TABFPPO$NB_Adult_Females_V2=as.numeric(scale(TABFPPO$NB_Adult_Females_V2))
TABFPPO$Number_Conc_Update6=scale(TABFPPO$Number_Conc_Update6)

y=TABFPPO$FPPO

MOD=bglmer(y ~ NDVI_Seasonality4 + Dif_NDVI3 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             Number_Conc_Update6 + (Number_Conc_Update6:Relative.rank) + 
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABFPPO, family='binomial')
FITTED=fitted(MOD)

TAB=read.csv2('Cycle resumption tab.csv.csv')
TAB=TAB[which(is.na(TAB$FPPO)==FALSE),]
TAB=TAB[which(is.na(TAB$NB_Adult_Females_V2)==FALSE),]
TAB=TAB[which(is.na(TAB$Number_Conc_Update6)==FALSE),]
tabsansRankNA2=TAB[which(is.na(TAB$Rank)==FALSE),]
tabsansRankNA2$Rank=factor(tabsansRankNA2$Rank,c('High-Ranking','Mid-Ranking','Low-Ranking'))
tabsansRankNA2=cbind(tabsansRankNA2,FITTED)
FPPO=ggplot(tabsansRankNA2, aes(x=Number_Conc_Update6, y=FITTED)) + geom_point(aes(fill=Rank, color=Rank)) + 
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, method.args = list(family=binomial),  aes(fill=Rank, color=Rank),size=1.8) + 
  labs(x='Group reproductive synchrony',y='Probability of cycle resumption') +
  scale_colour_colorblind() +
  theme(axis.title.x = element_blank(),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=15),
        legend.position="none")
FPPO


## Panel C (conception)


TABCONC=read.csv2('Conception table.csv')
TABCONC$Relative.rank=scale(TABCONC$Relative.rank,center=TRUE,scale=TRUE)
TABCONC$Conc=as.factor(TABCONC$Conc)
TABCONC$ID=factor(TABCONC$ID)
TABCONC$Dif_NDVI12=scale(TABCONC$Dif_NDVI12)
TABCONC$NDVI_Seasonality2=scale(TABCONC$NDVI_Seasonality2)
TABCONC$NB_Adult_Females_V2=as.numeric(scale(TABCONC$NB_Adult_Females_V2))
TABCONC$Number_Conc_Update4=scale(TABCONC$Number_Conc_Update4)
MOD=bglmer(Conc ~ NDVI_Seasonality2 + Dif_NDVI12 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             Number_Conc_Update4 + 
             (Number_Conc_Update4:Relative.rank)  +
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABCONC, family='binomial')
Fitted=fitted(MOD)
TAB=read.csv2('Conception table.csv')
TAB=TAB[which(is.na(TAB$Conc)==FALSE),]
TAB=TAB[which(is.na(TAB$NB_Adult_Females_V2)==FALSE),]
TAB=TAB[which(is.na(TAB$Number_Conc_Update4)==FALSE),]
tabsansRankNA2=TAB[which(is.na(TAB$Rank)==FALSE),]
tabsansRankNA2$Rank=factor(tabsansRankNA2$Rank,c('High-Ranking','Mid-Ranking','Low-Ranking'))
tabsansRankNA2=cbind(tabsansRankNA2,Fitted)
CONC=ggplot(tabsansRankNA2, aes(x=Number_Conc_Update4, y=Fitted)) + geom_point(aes(fill=Rank, color=Rank)) + 
  stat_smooth(method="glm", se=FALSE, fullrange=TRUE, method.args = list(family=binomial),  aes(fill=Rank, color=Rank),size=1.8) + 
  scale_y_continuous(limits=c(0,0.6),breaks=c(0,0.2,0.4)) + 
  labs(x='Group reproductive synchrony',y='Probability of conception') + 
  scale_colour_colorblind() +
  theme(axis.title.x = element_text(size=15),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=15),axis.text.y = element_text(size=15),
        legend.position='none')
CONC



plot_grid(IBI,FPPO,CONC,nrow=3,ncol=1, labels = "AUTO", label_size = 25, label_colour = 'black')
ggsave("FIG3.tiff", units="in", width=6, height=12, dpi=300, compression = 'jpeg')


############################################################################################################
#######                   FIGURE S1                ########################################################
############################################################################################################

#it has been created using another software (manually)

############################################################################################################
#######                   FIGURE S2                ########################################################
############################################################################################################


TABCONC=read.csv2('Conception table.csv')
TABCONC$Month=as.factor(TABCONC$Month)
TABCONC$Relative.rank=scale(TABCONC$Relative.rank,center=TRUE,scale=TRUE)
TABCONC$Conc=as.factor(TABCONC$Conc)
TABCONC$ID=factor(TABCONC$ID)
TABCONC$Dif_NDVI12=scale(TABCONC$Dif_NDVI12)
TABCONC$NDVI_Seasonality2=scale(TABCONC$NDVI_Seasonality2)
TABCONC$NB_Adult_Females_V2=as.numeric(scale(TABCONC$NB_Adult_Females_V2))
TABCONC$Number_Conc_Update4=scale(TABCONC$Number_Conc_Update4)

MOD=bglmer(Conc ~ NDVI_Seasonality2 + Dif_NDVI12 + NB_Adult_Females_V2 + Troop + Relative.rank + Parity +  
             Number_Conc_Update4 + 
             (Number_Conc_Update4:Relative.rank)  +
             (1|ID), glmerControl(optimizer='bobyqa'), data=TABCONC, family='binomial')
Fitted=fitted(MOD)

TABCONC=read.csv2('/Users/julesdezeure/ownCloud/Analyses_Chapter1/Post_Revision_AmNat/Analyses_Post_Revision_AmNat/TABModelCONC.csv')
TABCONC=TABCONC[which(is.na(TABCONC$Rank)==FALSE),]
TABCONC$Rank=factor(TABCONC$Rank,c('High-Ranking','Mid-Ranking','Low-Ranking'))
TABCONC=TABCONC[which(is.na(TABCONC$Conc)==FALSE),]
TABCONC=cbind(TABCONC,Fitted)

p=ggplot(TABCONC, aes(x=NDVI_Seasonality2, y=Fitted)) + geom_point() + 
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE,method.args = list(family=binomial), size=2,col='#009E73') + 
  labs(x='Seasonal NDVI over the past 2 months',y='Probability of conception') +
  theme(axis.title.x = element_text(size=20),axis.text.x = element_text(size=15),
        axis.title.y = element_text(size=20),axis.text.y = element_text(size=15))
p
