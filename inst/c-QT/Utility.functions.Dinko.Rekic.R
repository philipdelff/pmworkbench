## Utility functions for C-QT analysis
## Dinko Rekic
## 2017 03 03

## Data wrangling
library(tidyverse)
## Auxiliary package for ggplot2 (all plots)
library(cowplot)
## Auxiliary package for ggplot2 (Hysteresis plots)
library(ggrepel)


## Labels for graphics---------------------------------------------------------------------
conc.label <- expression(paste("Concentrations (µmol/l)"))
DELTAQTcF.label <- expression(paste(Delta,"QTcF (ms)"))
DELTADELTA.label <- expression(paste(Delta,Delta,"QTcF (ms)"))
TIME.label <- "Nominal hours"

## AZ colors-------------------------------------------------------------------------------
AZcol<-rep(c("#FFAB00", "#830051", "#003865", "#68D2DF", "#3C1053", "#C4D600", "#3F4445"),20)

## Set theme for plots---------------------------------------------------------------------
theme_set(theme_bw(base_size = 9))
theme_replace(plot.title   = element_text(size = rel(1.3), 
                                          face   = "bold",
                                          hjust  = 0,
                                          margin = margin(10, 0, 10, 0)),
              legend.position = "bottom",
              legend.title    = element_blank(),
              panel.grid      = element_blank())


## My summary funciton for using dplyr::summarize_().-------------------------------------- 
my.sum.fun<-funs(mean   = mean(. ,na.rm=T),
                 median = median(. , na.rm=T),
                 min    = min(. ,na.rm=T),
                 max    = max(. ,na.rm=T),
                 sd     = sd(. ,na.rm=T), 
                 n      = sum(!is.na(.)),
                 se     = sd(. ,na.rm=T)/sqrt(sum(!is.na(.))),
                 LCL    = mean(. ,na.rm=T)+qnorm(0.05)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))),
                 UCL    = mean(. ,na.rm=T)+qnorm(0.95)*(sd(. ,na.rm=T)/sqrt(sum(!is.na(.)))))


## Funciton to make dataset for matric plot-------------------------------------------- 
make.matrix.data<-function(qtpk)
{
  ## Variable needed to make a plot
  var.names<-c("PATIENT", "RR", "QT", "QTCF", "QTCB", "ACTIVE")  
  ## If PART not defined by user then make var PART==1
  if(prod(var.names %in% names(qtpk))==0){
    ## If varible not presentent retur error message
    stop("Needed variables are missing")

  } else {
    ## create long dataset
    matrix.data<- qtpk %>% 
      dplyr::select(PATIENT, RR, QT, QTCF, QTCB, PART, ACTIVE) %>%
      gather(value = value, key = type, -RR, -PATIENT, -ACTIVE, -PART)
  }
  
  return(matrix.data)
  
}
  

## Funciton to make the matrix plot----------------------------------------------------
plot.matrix<-function(matrix.data){
  matrix.plot. <- ggplot(data=matrix.data, aes(RR, value, col=ACTIVE))+
      geom_point(alpha=0.2)+
      facet_grid(PART~type, scales = "free_y")+
      geom_smooth(aes(RR, value, fill=ACTIVE),
                  size=1, method="lm", se=T)+
      labs(title="RR vs QT, QTcB and QTcF", 
           x="RR interval", 
           y="QT interval (ms)", 
           subtitle="", 
           caption="")+
      scale_color_manual(values=AZcol)+
      scale_fill_manual(values=AZcol)
    
    return(matrix.plot.)
    
}

## Function to make dataset for by time plot------------------------------------------

make.by.time.data<-function(qtpk){
  ## Variable needed to make a plot
  var.names<-c("PATIENT", "NOMTIME", "DHR", "DV", "DOSE_TRT", "ACTIVE")
  ## If PART not defined by user then make var PART==1
  if(prod(var.names %in% names(qtpk))==0){
    ## If varible not presentent retur error message
    stop("Needed variables are missing" )

  } else {
    by.time.dat<-qtpk %>% 
      # Make long dataset
      dplyr::select(PATIENT, NOMTIME, DHR, DQTCF, DV, ACTIVE, DOSE_TRT, PART) %>%
      gather(key = KEY, value = VALUE, -PATIENT, -ACTIVE, -NOMTIME, -DOSE_TRT, -PART) %>%
      ## Group by these factors
      group_by(KEY, PART, DOSE_TRT, ACTIVE, NOMTIME) %>%
      ## Summarize HR, PK, and QT by factors using my simmary function 
      summarise_each(my.sum.fun, VALUE) %>%
      ## Show mean and Sd for PK and CI for HR and QT
      mutate(P.UCL = ifelse(KEY  == "DV", mean+sd, UCL),
             P.LCL = ifelse(KEY  == "DV", mean-sd, LCL)) %>%
      ungroup() %>%
      ## Labels for facets in plot
      mutate(KEY = ifelse(KEY == "DV", "Concentrations (µmol/l), mean \u00b1SD",
                          ifelse(KEY == "DHR", "\u0394Heart rate (bpm), mean and 90% CI",
                                 "\u0394QTcF (ms), mean and 90% CI")))
    
    by.time.dat<-by.time.dat %>% filter(!is.na(mean)) ## All data not colelcted for all time points
    return(by.time.dat)
  }
}
  
## Funciton to make by time plot-------------------------------------------------------------

plot.by.time<-function(by.time.dat){
  ## Code to make figure by time figure
  By.Time.Figure<-ggplot(data=by.time.dat, aes(x=NOMTIME,y=mean, col=ACTIVE))+
    geom_linerange(aes(ymin=P.LCL, ymax=P.UCL))+
    geom_point()+
    geom_line(aes(x=NOMTIME,y=mean, linetype=ACTIVE, group=DOSE_TRT))+
    labs(title="Concentration, DQTCF, and HR versus nominal time",
         x="Time (h)", y="",
         subtitle="Placebo (SAD or MAD) versus all doses", 
         caption="data source:D7550C00001NM_QT.csv")+
    facet_grid(KEY~PART, scales = "free_y")+
    scale_color_manual(values=AZcol)+
    scale_fill_manual(values=AZcol)
  
  return(By.Time.Figure)
}

## Function to make hysteresis data

make.hyst.dat<-function(qtpk){
  ## Variable needed to make a plot
  var.names<-c("NOMTIME", "DDQTCF", "DV", "DOSE_TRT", "ACTIVE")
  ## If PART not defined by user then make var PART==1
  if(prod(var.names %in% names(qtpk))==0){
    ## If varible not presentent retur error message
    stop("Needed variables are missing;" )  
  } else {
    
    hyst.dat<- qtpk %>%
      ## grouping factors 
      group_by(DOSE_TRT, ACTIVE, PART, NOMTIME) %>% 
      filter(!is.na(QT)) %>% 
      ## summarize ddQT and PK by factors using my summary function 
      summarise_each(my.sum.fun, DV, DDQTCF)
  }
    
    return(hyst.dat)
}

# 

# 
qtpk %>% make.matrix.data(.) %>% plot.matrix(.)

qtpk %>% make.by.time.data(.) %>% plot.by.time(.)
