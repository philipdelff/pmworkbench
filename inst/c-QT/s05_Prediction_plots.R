## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 03.07.2017
## Reviwer:

## Script no 5 prediciton plots using the prespecified model

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

## sourse Prediciton script------------------------------------------------------------------------

source("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/scripts/s04_20170306_Prediction.R")

## Settings for plots------------------------------------------------------------------------------
## Labels for graphics
conc.label <- expression(paste("Concentrations (µmol/l)"))
DELTAQTcF.label <- expression(paste(Delta,"QTcF (ms)"))
DELTADELTA.label <- expression(paste(Delta,Delta,"QTcF (ms)"))
TIME.label <- "Nominal hours"

## AZ col
AZcol<-rep(c("#FFAB00", "#830051", "#003865", "#68D2DF", "#3C1053", "#C4D600", "#3F4445"),20)


## Set theme for plots
theme_set(theme_bw(base_size = 12))
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

## Wraper for titles 
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

## Data for binplot--------------------------------------------------------------------------------

## Bin plot CONC-DDQT
## Summary by concentrion decile data set 
Bin.qtpk<- qtpk %>%
  filter(!is.na(DV) & !is.na(QT) & TRTID==1) %>%
  group_by(Decile=ntile(DV, 10), TRTID) %>% 
  summarise_each(my.sum.fun, DV, DDQTCF, DQTCF) %>% 
  mutate(TRTID=as.factor(TRTID))


## Create data for plotting bin width
plot.bins<-expand.grid(cutpoints=c(Bin.qtpk$DV_min, max(qtpk$DV, na.rm=T)),
                       y.plot=min(qtpk$DDQTCF, na.rm=T)*1.2)


## Prediction based on lsmeans DQTCF and DDQTCF----------------------------------------------------------------------
## 
## Concentrations should be ordered from min to max


Cmax.ref.grid <- ref.grid(Pre_spec_model,                     # Model 
                          at=list(DV=seq(from=0, 
                                         to=max(exp.resp$DV, na.rm=T), 
                                         length.out = 100)),  # Concentration to make prediction
                          TRTID=c(0,1))                       # Placebo and active


## Generate LS-Means estimated of baseline adjusted QTcF for active and placebo
## by CONC, predictions are averaged over the levels of TIME. 
Cmax.lsm<-lsmeans::lsmeans(Cmax.ref.grid, c("DV", "TRTID"))


## Get DQTCF estimates
Cmax.dQTcF<-summary(Cmax.lsm, level=0.9) %>% filter(TRTID==1)

## Get  DDQTCF estimates
Cmax.ddQTcF.tmp<-summary(contrast(Cmax.lsm, method = "trt.vs.ctrl1"),
                         infer = c(TRUE, FALSE), 
                         level = .90, adjust = "none") # Do not adjust for multiplicity

## Select only relevant contrasts
Cmax.ddQTcF<-Cmax.ddQTcF.tmp[(nrow(Cmax.ddQTcF.tmp)/2):nrow(Cmax.ddQTcF.tmp)+1,]
Cmax.ddQTcF<- Cmax.ddQTcF %>% 
  separate(contrast, sep=",", into="DV", extra="drop") %>%
  mutate(DV=as.numeric(DV), 
         lsmean=estimate)

## Plotting DQTCF-----------------------------------------------------------------------

Pred.plot.dQTcF<-ggplot()+
  geom_ribbon(data=Cmax.dQTcF, 
              aes(x=DV, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=Cmax.dQTcF,  aes(x=DV, y=lsmean), size=1, col="black")+
  geom_hline(aes(yintercept=0), linetype=2)+
  labs(x=conc.label, y=DELTAQTcF.label)+
  scale_color_manual(values=AZcol[2:4])+
  theme(legend.position="none")


Pred.plot.dQTcF.bin<-Pred.plot.dQTcF+
   geom_pointrange(data=Bin.qtpk, aes(x=DV_median,
                                     ymax=DQTCF_UCL,
                                     ymin= DQTCF_LCL,
                                     y=DQTCF_mean,
                                     col=TRTID)) +
  geom_point(data=plot.bins, aes(y=y.plot, x=cutpoints), shape="|", size=2 )+
  geom_line(data=plot.bins, aes(y=y.plot, x=cutpoints), size=0.25)

Pred.plot.dQTcF.bin.obs<-Pred.plot.dQTcF.bin+
  geom_point(data=qtpk, aes(x=DV, y=DQTCF, col=as.factor(TRTID)), alpha=0.2)+
  labs(title=wrapper("Non-significant exposure-\u0394QTcF relationship",75),
       subtitle=(wrapper("The black line represents predictions from the prespecified C-QTcF model.The shaded area represent the 90% CI of the prediction. The points and and bars repesent \u0394QTcF mean and 90% CI at the meadian concentration in a bin", 120)))

Pred.plot.dQTcF.bin.obs

ggsave("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Pred.plot.dQTcF.bin.obs.png",
       width = 9, height = 6)
  



## Plotting DDQTCF-----------------------------------------------------------------------

Pred.plot.ddQTcF<-ggplot()+
  geom_ribbon(data=Cmax.ddQTcF, 
              aes(x=DV, ymin= lower.CL, ymax=upper.CL), fill="black", alpha=0.1)+
  geom_line(data=Cmax.ddQTcF,  aes(x=DV, y=lsmean), size=1, col="black")+
 
  geom_hline(aes(yintercept=10), linetype=2)+
  geom_hline(aes(yintercept=0), linetype=1)+
  labs(x=conc.label, y=DELTADELTA.label)+
  scale_color_manual(values=AZcol[2:4])+
  theme(legend.position="none")

Pred.plot.ddQTcF.bins<-Pred.plot.ddQTcF+
  geom_point(data=qtpk, aes(x=DV, y=DDQTCF, col=as.factor(TRTID)), alpha=0.2)+
  geom_pointrange(data=Bin.qtpk, aes(x=DV_median,
                                     ymax=DDQTCF_UCL,
                                     ymin= DDQTCF_LCL,
                                     y=DDQTCF_mean,
                                     col=TRTID)) +
  geom_point(data=plot.bins, aes(y=y.plot, x=cutpoints), shape="|", size=2 )+
  geom_line(data=plot.bins, aes(y=y.plot, x=cutpoints), size=0.25)


## Get Cmax data for label
label.df<-Prediciton.table %>%  filter(QTcF=="Baseline and placebo corrected")



Pred.plot.ddQTcF.text<-Pred.plot.ddQTcF+
  geom_segment(data=label.df,  aes(x=6.88, xend=6.88, y=-Inf, yend=5.27), col="red",  size=1)+
  geom_segment(data=label.df,  aes(x=6.88, xend=0, y=5.27, yend=5.27), col="red",
               arrow = arrow(length = unit(0.5, "cm")), size=1)+
  geom_label(data=label.df, aes(x=0, y=6), label="Upper 90% CI at supratherapetic Cmax:5.27 ms", 
             fill="red", col="white", hjust=0)+
  labs(title=wrapper("The upper 90% CI of \u0394\u0394QTcF is estimated <10 ms at supratherapetic Cmax",75),
       subtitle=wrapper("The black line represents predictions from the prespecified C-QTcF model. The shaded area
represent the 90% CI of the prediction, Supratherapetic 1200 mg Cmax: 6.8 \u00B5mol/L",120))

Pred.plot.ddQTcF.text
ggsave("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Pred.plot.ddQTcF.text.png",
       width = 9, height = 6)


