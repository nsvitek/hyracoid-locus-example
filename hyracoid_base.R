#This script is intended to be used after linear measurements are developed

#two options for analyses: analyze upper molars or lower molars?
arcade <- "lowers"
# arcade <- "uppers"

#read in data -----------
library(dplyr) #organize
library(ggplot2) #plot
library(readxl) #read procavia data
library(reshape2)
library(caret)
library(khroma) #for  Paul Tol color options
library(MASS) #for LDA in post-review analyses

# locateData <- "C:/Users/nsvit/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position"
locateData<-"D:/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position"

locateScripts<-"C:/scripts/hyracoid-locus-example"

literature_format_script<-"literature_format_linear_data.R"
#settings for plotting -------------------------------
single.column.width<-3.27
double.column.width<-6.61
scale_n<-as.character(colour("sunset")(10))
scale_locus<-as.character(colour("muted")(3))

#read measurement spreadsheets
# UPPERS ---------
if(arcade=="uppers"){
  data.procavia<-read_excel(paste(locateData,"modern-procavia-metadata-uppers.xlsx",sep="/"))
  data.comparative<-read_excel(paste(locateData,"hyracoidea_fayum_locus_comparisons.xlsx",sep="/"),sheet="Measurements-Upper")
  linear_format_script<-"procavia_format_linear_data_uppers.R"
  comparative_format_script<-"fayum_format_linear_data_uppers.R"
  
  #Analyses are the same for upper and lowers, so set an appropriate directory for results
  setwd(locateData)
  setwd("output_uppers")
}

# #LOWERS ----------
if(arcade=="lowers"){
  data.procavia<-read_excel(paste(locateData,"modern-procavia-metadata.xlsx",sep="/"))
  data.comparative<-read_excel(paste(locateData,"hyracoidea_fayum_locus_comparisons.xlsx",sep="/"),sheet="Measurements-Lower")
  linear_format_script<-"procavia_format_linear_data_lowers.R"
  comparative_format_script<-"fayum_format_linear_data_lowers.R"
  
  #Analyses are the same for upper and lowers, so set an appropriate directory for results
  setwd(locateData)
  setwd("output_lowers")
}

# data cleanup ------
#source the scripts to read, filter, format, and merge data
source(paste(locateScripts,"hyracoid_merge_data.R",sep="/"))

# analytical steps -------
#then source the shared analytical steps 
source(paste(locateScripts,"hyracoid_analyses.R",sep="/"))
