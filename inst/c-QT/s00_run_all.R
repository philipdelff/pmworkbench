## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 19.11.2019
## Reviwer: 

## Script no 0, runn all other scripts

## houskeeping
rm(list = ls())

## Path------------------------------------------------------------------------
path<-"scripts/"

## Source----------------------------------------------------------------------

source(paste0(path, "s01_20170302_derive.data.R"))

source(paste0(path, "s02_20170302_EDA_plots.R"))

source(paste0(path, "s03_20170302_Modeling.R"))

source(paste0(path, "s04_20170306_Prediction.R"))

source(paste0(path, "s05_20170306_Prediction_plots.R"))
