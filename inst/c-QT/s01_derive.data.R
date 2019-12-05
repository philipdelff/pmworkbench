## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 19.11.2019
## Reviwer:

## Script no 1
## Deriving analysis data 

## houskeeping
#rm(list = ls())

## Import packages
library(knitr)
## Model fitting
library(lme4)  
## Estimation of model parameter p-values and CI
library(lmerTest)
library(pbkrtest)
## Estimation of ddQTcF by contrasts
library(lsmeans)
## Handeling of model outputs
library(broom)
## Graphics and data wrangling
library(tidyverse)
## Graphics
library(ggthemes)
## Graphics
library(cowplot)
## Graphics
library(ggrepel)


## Import dataset
#raw.dat<-read.csv("//samba-hpc/QCP_MODELING_OLD/CVMD/AZD5718/PMX/Activity_01_SMAD/Data/QT/D7550C00001NM_QT.csv")
raw.dat<-read_csv("SourceData/D7550C00001NM_QT.csv")

## Labels for graphics
conc.label <- expression(paste("Concentrations (umol/l)"))


DELTAQTcF.label <- expression(paste(Delta,"QTcF (ms)"))
DELTADELTA.label <- expression(paste(Delta,Delta,"QTcF (ms)"))
TIME.label <- "Nominal hours"

## AZ col
AZcol<-rep(c("#FFAB00", "#830051", "#003865", "#68D2DF", "#3C1053", "#C4D600", "#3F4445"),20)


## Set theme for plots
theme_set(theme_bw(base_size = 9))
theme_replace(plot.title   = element_text(size = rel(1.3), 
                                          face   = "bold",
                                          hjust  = 0,
                                          margin = margin(10, 0, 10, 0)),
              legend.position = "bottom",
              legend.title    = element_blank(),
              panel.grid      = element_blank())

## My summary funciton for using dplyr::summarize_(). 
my.sum.fun<-funs(mean   = mean(. ,na.rm=T),
                 median = median(. , na.rm=T),
                 min    = min(. ,na.rm=T),
                 max    = max(. ,na.rm=T),
                 sd     = sd(. ,na.rm=T), 
                 n      = sum(!is.na(.)),
                 se     = sd(. ,na.rm=T)/sqrt(sum(!is.na(.))),
                 LCL    = mean(. ,na.rm=T)+qnorm(0.05)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))),
                 UCL    = mean(. ,na.rm=T)+qnorm(0.95)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))))


## Cohort 1-8 = SAD, 9-12=MAD  
qtpk<-raw.dat %>% 
  #filter(MDV==0) %>% 
  filter(EVID==0) %>% 
  mutate(DV=as.numeric(as.character(DV)),
         QT=as.numeric(as.character(QT)),
         RR=as.numeric(as.character(RR)),
         DHR=as.numeric(as.character(DHR)),
         BQTCF=as.numeric(as.character(BQTCF)),
         TRTID=as.factor(TRTID),
         QTCF=as.numeric(as.character(QTCF)),
         DQTCF=as.numeric(as.character(DQTCF)),
         DOSE=as.factor(DOSE), 
         ##Optional: get Bezett's corrected QT
         QTCB=QT/((RR/1000)^0.5), 
         ## Imput zero concentrations in the placebo group
         DV=ifelse(TRTID==0,0,DV),
         ## Imput zero concentration at nominal time 0
         DV=ifelse(NOMTIME==0 & is.na(DV), 0, DV))



## Calculate mean placebo DQTCF
Placebo.QTCF<-qtpk %>% 
  filter(DOSE==0) %>% 
  group_by(NOMTIME, PART) %>% 
  summarise(Placebo.dQTcF.mean=mean(DQTCF, na.rm=T))

## Merge mean placebo and qtpk
qtpk<-qtpk %>% left_join(., Placebo.QTCF)


## Calculate DDQTCF
qtpk<-qtpk %>% mutate(DDQTCF=DQTCF-Placebo.dQTcF.mean)

## Calculate populaiton mean baseline QTCF and ad to qtpk dataset
popMean.BQTCF<-qtpk %>% 
  distinct(PATIENT, .keep_all=T) %>%
  summarise(popMean.BQTCF=mean(BQTCF))


qtpk$QTcF.mB<-popMean.BQTCF[[1]]


## Make informative variable names
qtpk<-qtpk %>% 
  mutate(DOSE_TRT_FOOD=paste(TRT, DOSE, "mg", FASTED),
         ACTIVE=ifelse(TRTID==0,"Placebo", "Active"),
         MAD_SAD_DAY=ifelse(PART=="MAD" & (DAY==1 | DAY==2)," MAD Day 1 dose",
                     ifelse(PART=="MAD" & DAY>=10  ,"MAD Day 10 dose",
                     ifelse(PART=="MAD" & DAY==9    ,"MAD Day 9 dose",
                     ifelse(PART=="SAD"            ,"SAD Day 1 dose",  "Not used"))))) %>% 
  select(PATIENT, PART, TRT, DOSE, NOMTIME, DV, TRTID, QT, QTCF, DDQTCF,DQTCF, BQTCF, QTCB,  
         DHR, HR,RR, DHR, QTcF.mB, DOSE_TRT_FOOD, ACTIVE, DAY, MAD_SAD_DAY, FASTED) %>% 
  filter(MAD_SAD_DAY!="Not used" & MAD_SAD_DAY!="MAD Day 9 dose") 


write_csv(qtpk, "DerivedData/qtpk.csv")

qtpk %>%  group_by(MAD_SAD_DAY, DOSE_TRT_FOOD) %>%
  summarise(Subjects_n=n_distinct(PATIENT),
            PK_n=sum(!is.na(DV)), 
            QT_n=sum(!is.na(QT))) %>%
  write_csv(., "DerivedData/Observations.csv")

  
