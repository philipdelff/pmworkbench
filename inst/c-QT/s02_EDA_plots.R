## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 19.11.2019
## Reviwer:

## Script no 2 graphical analysis




## Import packages-------------------------------------------

## houskeeping
#rm(list = ls())

## Data wrangling
library(tidyverse)
## Auxiliary package for ggplot2 (all plots)
library(cowplot)
## Auxiliary package for ggplot2 (Hysteresis plots)
library(ggrepel)


## Caption 
date_time<-Sys.time()
script_name<-"s02_20170302_EDA_plots.R"
data_source <- "data source: D7550C00001NM_QT.csv"
caption <- paste0(script_name, "\n", data_source, "\n", date_time)

## Settings for plots-----------------------------------------------------
## Labels for graphics
conc.label <- expression(paste("Concentrations (umol/l)"))
DELTAQTcF.label <- expression(paste(Delta,"QTcF (ms)"))
DELTADELTA.label <- expression(paste(Delta,Delta,"QTcF (ms)"))
TIME.label <- "Nominal hours"

## AZ colors
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


## Wraper for titles 
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

## Import dataset-----------------------------------------------
qtpk<-read_csv("DerivedData/qtpk.csv")

## Change order of factors
qtpk$ACTIVE<-factor(qtpk$ACTIVE, levels = c("Placebo", "Active"))

#View(qtpk)



## Matric plot to verify if HR correction is adequate --------------------------------
mplot.data<- qtpk %>% 
  dplyr::select(PATIENT, RR, QT, QTCF, QTCB, PART, ACTIVE) %>%
  gather(value = value, key = type, -RR, -PATIENT, -ACTIVE, -PART)


## Generate a matrix plot
m.plot <- ggplot(data=mplot.data, aes(RR, value, col=ACTIVE))+
  geom_point(alpha=0.2)+
  facet_grid(PART~type, scales = "free_y")+
  geom_smooth(aes(RR, value, fill=ACTIVE),
              size=1, method="lm", se=T)+
  labs(title="QTcF is adequatly corrected for HR", 
       x="RR interval", 
       y="QT interval (ms)", 
       subtitle="", 
       caption=caption)+
  scale_color_manual(values=AZcol)+
  scale_fill_manual(values=AZcol)

m.plot

ggsave("Results/mplot.png", width = 9, height = 6)

## Concentration dHR and dQT by time figure ------------------------------------------

## prepare summary data by time and treatmetn and part 

Sum.qtpk<-qtpk %>%
  # Make long dataset
  dplyr::select(PATIENT, NOMTIME, DHR, DQTCF, DV, ACTIVE, DOSE_TRT_FOOD, PART, MAD_SAD_DAY, FASTED) %>%
  gather(key = KEY, value = VALUE, -PATIENT, -ACTIVE, -NOMTIME, -DOSE_TRT_FOOD, -PART, -MAD_SAD_DAY, -FASTED) %>%
  ## Group by these factors
  group_by(KEY, PART,MAD_SAD_DAY, DOSE_TRT_FOOD,  ACTIVE, NOMTIME) %>%
  ## Summarize HR, PK, and QT by factors using my simmary function 
  summarise_at(vars(VALUE), my.sum.fun) %>%
  ## Show mean and Sd for PK and CI for HR and QT
  mutate(P.UCL = ifelse(KEY  == "DV", mean+sd, UCL),
         P.LCL = ifelse(KEY  == "DV", mean-sd, LCL)) %>%
  ungroup() %>%
  ## Labels for facets in plot
  mutate(KEY = ifelse(KEY == "DV", "Concentration (umol/l)",
               ifelse(KEY == "DHR", "\u0394Heart rate (bpm)",
                      "\u0394QTcF (ms)")))

## All data not colelcted for all time points
Sum.qtpk<-Sum.qtpk %>% filter(!is.na(mean)) 

## Make figure by time figure
By.Time.Figure<-ggplot(data=Sum.qtpk, aes(x=NOMTIME,y=mean, col=ACTIVE))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=NOMTIME,y=mean, linetype=ACTIVE, group=DOSE_TRT_FOOD))+
  labs(title="Concentration, DQTCF, and HR versus nominal time",
       x="Time (h)", y="",
       subtitle=wrapper("Placebo (SAD or MAD) versus all doses. Concentrations are shown as mean \u00b1SD, HR and BP are shown as mean and 90% CI",120), 
       caption=caption)+
  facet_grid(KEY~MAD_SAD_DAY, scales = "free_y")+
  scale_color_manual(values=AZcol)+
  scale_fill_manual(values=AZcol) 



ggsave("Results/Time_Plot1.png",By.Time.Figure, width = 9, height = 6)

By.Time.Figure2<-ggplot(data=Sum.qtpk, aes(x=NOMTIME,y=mean, col=DOSE_TRT_FOOD))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=NOMTIME,y=mean, linetype=DOSE_TRT_FOOD, group=DOSE_TRT_FOOD))+
  labs(title="Concentration, DQTCF, and HR versus nominal time",
       x="Time (h)", y="",
       subtitle=wrapper("All doses and Placebo (SAD or MAD). Concentrations are shown as mean \u00b1SD, HR and BP are shown as mean and 90% CI",120), 
       caption=caption)+
  facet_grid(KEY~MAD_SAD_DAY, scales = "free_y")+
  scale_color_manual(values=AZcol)+
  scale_fill_manual(values=AZcol)



ggsave("Results/Time_Plot2.png", By.Time.Figure2, width = 9, height = 6)

## only keep highest dose
Sum.qtpk %>% 
  filter(((ACTIVE=="Placebo"| DOSE_TRT_FOOD=="AZD5718 oral suspension (amorphous) 600 mg FASTED") & PART=="MAD") | 
           ((ACTIVE=="Placebo"| DOSE_TRT_FOOD=="AZD5718 oral suspension (amorphous) 1200 mg FASTED") & PART=="SAD")) %>% 
  ggplot(data=., aes(x=NOMTIME,y=mean, col=ACTIVE))+
  geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
  geom_point()+
  geom_line(aes(x=NOMTIME,y=mean, linetype=DOSE_TRT_FOOD, group=DOSE_TRT_FOOD))+
  labs(title="Concentration, DQTCF, and HR versus nominal time",
       x="Time (h)", y="",
       subtitle=wrapper("Top dose versus Placebo (SAD or MAD). Concentrations are shown as mean \u00b1SD, HR and BP are shown as mean and 90% CI",120), 
       caption=caption)+
  facet_grid(KEY~MAD_SAD_DAY, scales = "free_y")+
  scale_color_manual(values=AZcol)+
  scale_fill_manual(values=AZcol) -> By.Time.Figure3


ggsave("Results/Time_Plot3.png", By.Time.Figure3, width = 9, height = 6)


## Hysteresis plot -----------------------------------------------------------------------------------------------
hyst.qtpk<- qtpk %>%
  ## grouping factors 
  group_by(DOSE_TRT_FOOD, ACTIVE, MAD_SAD_DAY, NOMTIME) %>% 
  filter(!is.na(QT)) %>% 
  ## summarize ddQT and PK by factors using my summary function 
  summarise_at(vars(DV, DDQTCF), my.sum.fun)



## Hysteresis plot SAD part
hyst.qtpk %>%  filter(MAD_SAD_DAY=="SAD Day 1 dose" & ACTIVE!="Placebo") %>% 
ggplot(data=., aes(x=DV_mean, y=DDQTCF_mean, col=DOSE_TRT_FOOD))+
  geom_path(alpha=0.7)+
  geom_pointrange(aes(ymin=DDQTCF_LCL, ymax=DDQTCF_UCL), alpha=0.5)+
  geom_text_repel(aes(label=paste0(as.factor(NOMTIME), " h")))+
  labs(x=conc.label, y=DELTADELTA.label,
       title="No apparent evidence of hysteresis in SAD",
       subtitle="Mean Placebo corrected and baseline adjusted QTcF versus mean concentration at each nominal time",
       caption=caption)+
  scale_color_manual(values=AZcol)+
  theme(legend.position="none")+
  facet_wrap(~DOSE_TRT_FOOD, scale="free") -> Hysteresis.plot.SAD1

Hysteresis.plot.SAD1

ggsave("Results/Hysteresis.plot.SAD.png",Hysteresis.plot.SAD1, width = 10, height = 6)

## Hysteresis plot MAD part
hyst.qtpk %>%  filter(MAD_SAD_DAY=="MAD Day 1 dose" & ACTIVE!="Placebo") %>% 
ggplot(data=., aes(x=DV_mean, y=DDQTCF_mean, col=DOSE_TRT_FOOD))+
  geom_path(alpha=0.7)+
  geom_pointrange(aes(ymin=DDQTCF_LCL, ymax=DDQTCF_UCL), alpha=0.5)+
  geom_text_repel(aes(label=paste0(as.factor(NOMTIME), " h")))+
  labs(x=conc.label, y=DELTADELTA.label,
       title="No apparent evidence of hysteresis in MAD day 1",
       subtitle="Mean Placebo corrected and baseline adjusted QTcF versus mean concentration at each nominal time",
       caption=caption)+
  scale_color_manual(values=AZcol)+
  theme(legend.position="none")+
  facet_wrap(~DOSE_TRT_FOOD, labeller = "label_both", scale="free") ->Hysteresis.plot.MAD1


ggsave("Results/Hysteresis.plot.MAD.png",Hysteresis.plot.MAD1, width = 8, height = 6)

## Hysteresis plot MAD part
hyst.qtpk %>%  filter(MAD_SAD_DAY=="MAD Day 10 dose" & ACTIVE!="Placebo") %>% 
  ggplot(data=., aes(x=DV_mean, y=DDQTCF_mean, col=DOSE_TRT_FOOD))+
  geom_path(alpha=0.7)+
  geom_pointrange(aes(ymin=DDQTCF_LCL, ymax=DDQTCF_UCL), alpha=0.5)+
  geom_text_repel(aes(label=paste0(as.factor(NOMTIME), " h")))+
  labs(x=conc.label, y=DELTADELTA.label,
       title="No apparent evidence of hysteresis in MAD day 10",
       subtitle="Mean Placebo corrected and baseline adjusted QTcF versus mean concentration at each nominal time",
       caption=caption)+
  scale_color_manual(values=AZcol)+
  theme(legend.position="none")+
  facet_wrap(~DOSE_TRT_FOOD, labeller = "label_both", scale="free") ->Hysteresis.plot.MAD10


ggsave("Results/Hysteresis.plot.MAD_10.png",Hysteresis.plot.MAD10, width = 8, height = 6)


## Bin plot CONC-DDQT --------------------------------------------------------------------------

## create bins (n=10) and summarize PK, DDQTCF, and DQTCF using my summary funciton 
Bin.qtpk<- qtpk %>%
  ## remove missing DV, and QT and Placebo
  filter(!is.na(DV) & !is.na(QT) & ACTIVE=="Active") %>%
  ## Create Deciles and use that as a factor for summary function 
  group_by(Decile=ntile(DV, 10), ACTIVE) %>% 
  summarise_at(vars(DV, DDQTCF, DQTCF), my.sum.fun)

## Create data for plotting bin width
## Shows the cutpoints of the bins (min and max)
plot.bins<-expand.grid(cutpoints=c(Bin.qtpk$DV_min, max(qtpk$DV, na.rm=T)),
                       y.plot=min(qtpk$DDQTCF, na.rm=T)*1.2)

## Code for creating Exposure response bin plot
Explor.plot<-ggplot()+
  geom_point(data=qtpk, aes(x=DV, y=DDQTCF, col=ACTIVE), alpha=0.2)+
  geom_smooth(data=qtpk, method="lm",se=F, col="black", size=0.5, linetype=2, aes(x=DV, y=DDQTCF))+
  geom_smooth(data=qtpk, method="loess",se=F, col="#f44336", aes(x=DV, y=DDQTCF))+
  geom_pointrange(data=Bin.qtpk, aes(x=DV_median, 
                                     ymax=DDQTCF_mean+DDQTCF_sd, 
                                     ymin=DDQTCF_mean-DDQTCF_sd,
                                     y=DDQTCF_mean, 
                                     col=ACTIVE)) +
  geom_point(data=plot.bins, aes(y=y.plot, x=cutpoints), shape="|", size=2 )+
  geom_line(data=plot.bins, aes(y=y.plot, x=cutpoints), size=0.25)+
  labs(x=conc.label, y=DELTADELTA.label,
       title="No major non-linear trend is seen",
       subtitle="Exposure QTcF bin plot",
       caption="")+
  scale_color_manual(values=AZcol)+
  theme(legend.position="none")

Explor.plot2<-Explor.plot+scale_x_log10()+labs(caption=caption)

plot_grid(Explor.plot, Explor.plot2, align = "h")
ggsave("Results/Exploratory_DQTCF.png",Explor.plot2, width = 10, height = 6)
