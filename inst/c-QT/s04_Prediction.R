## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 03.02.2017
## Reviwer:

## Script no 4 prediciton table using the prespecified model

## houskeeping
#rm(list = ls())

## Import packages-------------------------------------------

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


## sourse fitting script------------------------------------------------------------------------

source("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/scripts/s03_20170302_Modeling.R")


##Calculation of Cmax, and Tmax-----------------------------------------------------------
Cmax.df<-qtpk %>% 
  # Removes missing concentrations and exlude placebo
  filter(!is.na(DV) ,TRTID!=0) %>%                 
  group_by(PATIENT, MAD_SAD_DAY ,DOSE_TRT_FOOD ) %>%    
  ## Get Cmax and Tmax for each TREATMENT and ID
  mutate(Cmax.i = max(DV),                 
         Tmax   = NOMTIME[which.max(DV)]) %>%   
  ungroup() %>%   
  # Only select the first observation per ID and TREATMENT
  distinct(PATIENT,TRTID, DOSE_TRT_FOOD ,MAD_SAD_DAY, .keep_all=T) %>%  
  # Summary statistics by treatment
  group_by(DOSE_TRT_FOOD, MAD_SAD_DAY) %>%                        
  summarise(Cmax     = mean(Cmax.i, na.rm=T),        
            SD.Cmax  = sd(Cmax.i, na.rm = T),
            GMEAN    = exp(mean(log(Cmax.i))),
            GCV      = 100 * sqrt(exp(var(log(Cmax.i))) - 1),
            Cmax_i   = max(Cmax.i, na.rm=T), 
            Tmax.med = median(Tmax),
            Max.Tmax = max(Tmax, na.rm=T),
            Min.Tmax = min(Tmax, na.rm=T),
            n.id     = n_distinct(PATIENT)) %>%
  arrange(Cmax) %>% 
  mutate_each(funs(round(.,4)), -DOSE_TRT_FOOD, -MAD_SAD_DAY)

## Predict effect at Cmax----------------------------------------------------------------------------------------------
## Specify the concentrations of interest 
## Concentrations should be ordered from min to max

Cmax<-Cmax.df %>% ungroup() %>%  slice(which.max(GMEAN)) %>%  select(GMEAN)


Cmax.ref.grid <- ref.grid(Pre_spec_model,               # Model 
                          at=list(DV=c(0,Cmax$GMEAN)),  # Concentration to make prediction
                          TRTID=c(0,1))                 # Placebo and active


## Generate LS-Means estimated of baseline adjusted QTcF for active and placebo
## by CONC, predictions are averaged over the levels of TIME. 
Cmax.lsm<-lsmeans::lsmeans(Cmax.ref.grid, c("DV", "TRTID"))


## Get DQTCF estimates
Cmax.dQTcF<-summary(Cmax.lsm, level=0.9) %>% filter(TRTID==1, DV!=0)

dQTcF.tab<-Cmax.dQTcF %>% 
  dplyr::select(-SE, -df, -TRTID) %>%
  rename(`Concentration (µmol/l)`= DV,
         `Lower 90% CI`=lower.CL,
         `Upper 90% CI`= upper.CL, 
         Estimate=lsmean) %>%
  round(.,2) %>%
  mutate(Treatment="FLAPi (1200 mg SAD)",
         Exposure=c("Geometric Mean Cmax"), 
         QTcF="Baseline corrected") %>%
  filter(Exposure=="Geometric Mean Cmax") %>%
  dplyr::select(Treatment, QTcF, Exposure,  everything())

## Get  DDQTCF estimates
Cmax.ddQTcF.tmp<-summary(contrast(Cmax.lsm, method = "trt.vs.ctrl1"),
                         infer = c(TRUE, FALSE), 
                         level = .90, adjust = "none") # Do not adjust for multiplicity

## Select only relevant contrasts
Cmax.ddQTcF<-Cmax.ddQTcF.tmp[(nrow(Cmax.ddQTcF.tmp)/2):nrow(Cmax.ddQTcF.tmp)+1,]

## Formating of ddQTcF output table
Cmax.ddQTcF<- Cmax.ddQTcF %>% separate(contrast, sep=",", into="CONC", extra="drop") %>%
  mutate(CONC=as.numeric(CONC)) %>%
  filter(CONC!=0) %>%  #remove zero conc
  dplyr::select(-SE, -df) %>%
  rename(`Concentration (µmol/l)`= CONC,
         `Lower 90% CI`=lower.CL,
         `Upper 90% CI`= upper.CL, 
         Estimate=estimate) %>%
  round(.,2) %>%
  mutate(Treatment="FLAPi (1200 mg SAD)",
         Exposure=c("Geometric Mean Cmax"), 
         QTcF="Baseline and placebo corrected") %>%
  filter(Exposure=="Geometric Mean Cmax") %>%
  dplyr::select(Treatment, QTcF, Exposure,  everything())


Prediciton.table<-rbind(dQTcF.tab, Cmax.ddQTcF) 

write.csv(Prediciton.table, "C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Prediciton.table.csv")
