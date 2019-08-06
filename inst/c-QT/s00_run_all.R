## EDA of QT SAD and MAD 
## Author: Dinko Rekic
## Date: 03.07.2017
## Reviwer: 

## Script no 0, runn all other scripts

## houskeeping
rm(list = ls())

## Path------------------------------------------------------------------------
### PLEASE EDIT.
path<-"C:/Users/knhc208/OneDrive - AZCollaboration/Active Projects/AZ13702997_FLAP/cqt_20170203_sad_mad/scripts/"

## Source----------------------------------------------------------------------

source(paste0(path, "s01_derive.data.R"))

source(paste0(path, "s02_EDA_plots.R"))

source(paste0(path, "s03_Modeling.R"))

source(paste0(path, "s04_Prediction.R"))

source(paste0(path, "s05_Prediction_plots.R"))
