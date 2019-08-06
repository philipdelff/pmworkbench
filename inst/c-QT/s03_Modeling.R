## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 03.02.2017
## Reviwer:

## Script no 3 Fitting prespecified model

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

## Labels for graphics
conc.label <- expression(paste("Concentrations (µmol/l)"))


DELTAQTcF.label <- expression(paste(Delta,"QTcF (ms)"))
DELTADELTA.label <- expression(paste(Delta,Delta,"QTcF (ms)"))
TIME.label <- "Nominal hours"

## AZ col
AZcol<-rep(c("#FFAB00", "#830051", "#003865", "#68D2DF", "#3C1053", "#C4D600", "#3F4445"),20)

## Import dataset-----------------------------------------------
qtpk<-read_csv("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/DerivedData/qtpk.csv")






## Summerize analysis dataset --------------------------------------------------------------------------------
Missing.data<-qtpk %>%  group_by(DOSE_TRT_FOOD, MAD_SAD_DAY) %>%
  summarise(N.id=n_distinct(PATIENT), 
            N.PK.obs=sum(!is.na(DV)),
            N.PK.missing=sum(is.na(DV)),
            N.QT.obs=sum(!is.na(DQTCF)),
            N.QT.missing=sum(is.na(DQTCF))) 
 

write.csv(Missing.data,"C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/DerivedData/Misisng.data.csv")
      

## Data prep----------------------------------------------------------------------------------------------------
exp.resp<- qtpk %>%
  dplyr::select(PATIENT, NOMTIME, DQTCF, DOSE, DOSE_TRT_FOOD, MAD_SAD_DAY,
                DV, TRTID, QTcF.mB, 
                BQTCF, ACTIVE, PART) %>%
  ## Time as factor 
  mutate(NOMTIME=as.factor(NOMTIME), 
         TRTID=as.factor(TRTID)) %>% 
  filter(!is.na(DV)) %>% 
  filter(!is.na(DQTCF))

## Fitting prespecified model ----------------------------------------------------------------------------------
Pre_spec_model <- lmerTest::lmer(DQTCF ~ NOMTIME + TRTID + MAD_SAD_DAY + I(BQTCF-QTcF.mB) + DV + (DV|PATIENT),
                                 data = exp.resp)

summary(Pre_spec_model)
## Extraction of parameter estimates 
Par.tab.temp<-as.data.frame(coef(summary(Pre_spec_model, 
                                         ddf="Kenward-Roger"))) 

## Formating the table and calculating CI
Par.tab<-Par.tab.temp %>%
  rownames_to_column(var="Fixed effect parameter") %>%
  rename(p  = `Pr(>|t|)`, 
         se = `Std. Error`) %>%
  mutate(`Relative standard error (%)`=se/Estimate*100, 
         `Lower 95% CI` = Estimate + qt(0.025, df=df) * se,
         `Upper 95% CI` = Estimate + qt(0.975, df=df) * se) %>%
  mutate(p=ifelse(p<0.001, "<0.001", signif(p,3))) %>%
  dplyr::select(`Fixed effect parameter`, Estimate, `Lower 95% CI`,
                `Upper 95% CI`, `Relative standard error (%)`, `p`) %>% 
  mutate_each(funs(round(.,2)), -`Fixed effect parameter`, -p) 


## P value for slope not significant
Par.tab

## Save fixed effect table
write.csv(Par.tab, "C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Fixed_effects_table.csv", 
          row.names=F)


## Estimation of confidence intervals of random effect using profiling---------------------------------------
## Random effects are selected by specifying param="theta_"

# cc <- confint.merMod(Pre_spec_model, parm = "theta_", method = "profile", oldNames = FALSE)
# vc<-as.data.frame(VarCorr(Pre_spec_model),order="lower.tri")
# cc.vc<-cbind(vc,cc )
# 
# Random.effect.table<- cc.vc %>%
#   rownames_to_column(., var="Random effect parameter") %>%
#   dplyr::select(`Random effect parameter`, sdcor, `2.5 %`, `97.5 %`) %>%
#   rename(`Estimate (sd/corr)` = sdcor,
#          `Lower 95% CI` = `2.5 %`,
#          `Upper 95% CI` = `97.5 %`) %>%
#   mutate_each(funs(round(.,2)), -`Random effect parameter`)
# 
# 
# Random.effect.table
# # 
# ## Save random effect table
# write.csv(Random.effect.table, "C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Random_effects_table.csv",
#           row.names=F)


## Goodnes of fit plots------------------------------------------------------------------------------------------
## GOOF plot funciton 
GOF_Pre_Spec = function(Pre_spec_model, lg=F){
  ## Function to produce GOF plots. Takes only lme4 model object as argument
  ## The funciton assums that the model has the variables NOMTIME, DV, I(QTcF.B - QTcF.mB),
  ## TRTID and DQTCF
  ## Function depends on Broom and ggplot2 packages. 
  
  GOF.dat<-augment(Pre_spec_model) %>%
    mutate(Sd.rediduals=residuals(Pre_spec_model, scaled = T))
  
  IPRED.DV<-ggplot(data=GOF.dat, aes(x=.fitted, y=DQTCF))+
    geom_point(alpha=0.2)+
    geom_abline(intercept=0, slope=1)+
    geom_smooth(aes(x=.fitted, y=DQTCF),
                method="loess",size=1, se=T)+
    labs(x="Model-predicted value (ms)",
         y="Observed Value (ms)")+
    #coord_cartesian(ylim = c(, 140), xlim = c(-40, 140))+
    theme(legend.position="none")
  
  QQ.plot<-ggplot(GOF.dat, aes(sample = Sd.rediduals)) +
    stat_qq(alpha=0.2) +
    #coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 5))+
    geom_abline(slope=1)+
    labs(x="Theoretical Quantiles",
         y="Standardized Residuals") +
    scale_color_gdocs()+
    theme(legend.position="none")
  
  RES.CONC<-ggplot(GOF.dat, aes(x=DV, y= Sd.rediduals)) +
    geom_point(alpha=0.2) +
    geom_abline(slope=0) +
    labs(x=conc.label,
         y="Standardized Residuals") +
    geom_smooth()+
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    scale_color_gdocs() +
    theme(legend.position="none")
  
  RES.Base.QTcF<-ggplot(GOF.dat, aes(x=`I.BQTCF...QTcF.mB.`, y= Sd.rediduals)) + 
    geom_point(alpha=0.2) +
    geom_abline(slope=0) +
    labs(x="Centered Baseline QTcF (ms)",
         y="Standardized Residuals") +
    geom_smooth() +
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    scale_color_gdocs() +
    theme(legend.position="none")
  
  RES.TIME<-ggplot(GOF.dat, aes(x=NOMTIME, y= Sd.rediduals)) +
    geom_boxplot(notch = TRUE, outlier.color="white") +
    geom_abline(slope=0, col="#f44336") +
    labs(x="Nominal Time (hr)",
         y="Standardized Residuals") +
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    scale_color_gdocs() +
    theme(legend.position="none")
  
  RES.TIME<-ggplot(GOF.dat, aes(x=NOMTIME, y= Sd.rediduals)) +
    geom_boxplot(notch = TRUE, outlier.color="white")+
    geom_abline(slope=0, col="#f44336")+
    labs(x="Nominal Time (hr)",
         y="Standardized Residuals") +
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    #facet_wrap(~ACTIVE, labeller = "label_both") +
    theme(legend.position="none")
  
  RES.ACTIVE<-ggplot(GOF.dat, aes(x=TRTID, y= Sd.rediduals)) + 
    geom_boxplot(notch = TRUE, outlier.color="white") +
    geom_abline(slope=0, col="#f44336") +
    labs(x="Active treatment",
         y="Standardized Residuals") +
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    scale_color_gdocs() +
    theme(legend.position="none")
  
  RES.MAD_SAD<-ggplot(GOF.dat, aes(x=MAD_SAD_DAY, y= Sd.rediduals)) + 
    geom_boxplot(notch = TRUE, outlier.color="white") +
    geom_abline(slope=0, col="#f44336") +
    labs(x="Active treatment",
         y="Standardized Residuals") +
    geom_abline(intercept=1.96,slope=0,linetype="dashed") +
    geom_abline(intercept=-1.96,slope=0,linetype="dashed") +
    scale_color_gdocs() +
    theme(legend.position="none")
  
  plot_grid(IPRED.DV, QQ.plot, RES.CONC, RES.Base.QTcF, RES.TIME, RES.ACTIVE, RES.MAD_SAD, 
            ncol=2, align = "hv", labels="auto")
  
} 

GOF_Pre_Spec(Pre_spec_model)

ggsave("C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/Results/Prespecified_model_GOF.png",
       width = 8, height = 8)
