library(reshape2)
library(ggplot2)
library(readxl)
library(dplyr)
setwd("D:/Dropbox/Documents/research/Turkana/hyracoidea/records")
raw.measurements<-"../Topernawi_hyracoid_published_tooth_sizes.csv"

#estimating body masses of Topernawi species based on tooth dimensions
#Schwartz et al. 1995 on exactly this topic prefer 
#the Janis 1990 equations for hyracoids-perissodactyls for upper & lower m2 lengths
# equations as functions uppers ----
#m2/ length model
#log10(mass.kg) = log10(length.cm) * 2.90 + 1.209
# log(Y) = log(X) * b + log(a) + log(r)

estimate.from.uppers<-function(length.cm){
  log.mass<-log10(length.cm) * 2.90 + 1.209 
  return(log.mass)
}

#detransformed, corrected estimate should be:
#multiply the detransformed value (Y) by exp(var/2), where var is error variance (some constant)
#var = RMS = SEE^2, SEE = 0.476 for M2/

# exp(1) # = e
# exp(log(10)^2/2 * SEE^2) #CF detransformed using the antilog of natural log
# cf<-10^(log(10)^2/2 * log10(exp(1)) * SEE^2) #CF detransformed base 10 from base10 data transform

#but that equation using the published %SEE doesn't work, even with the following conversion from
#van valkenburgh
# #comes from log10(%SEE+100) - 2
# #so for example, if log standard error is 0.220
# LSE<-.220
# 10^(LSE+2)-100 #a 66% %SEE
# log10(166)-2
#check from Schwartz et al. 1995 Table 2
# (log(1000.4/968.4) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.11
# (log(187.1/181.1) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.11
# (log(9.7/9.4) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.11

estimate.from.uppers.corrected<-function(length.cm){
  # SEE.P = 0.476 #from Janis 1990 #I can't get this number to work.
  # SEE<-log10(100+(SEE.P*100))-2
  #working backwards from raw/corrected detransformed masses instead gives
  SEE<-0.111
  cf<-(log(10)^2/2 * log10(exp(1)) * SEE^2) #Correction factor, Table 2, Smith 1993
  log.mass<-log10(length.cm) * 2.90 + 1.209
  corrected.mass<-10^log.mass * 10^cf
  return(corrected.mass)
}

# #test, using a quick .csv made from measurements published in Schwartz et al. 1995
# data.schwartz<-read.csv("temporary_measurements.csv")
# #log mass estimate
# estimate.from.uppers(data.schwartz$M2.mm/10) 
# #detransformed estimate without correction
# 10^estimate.from.uppers(data.schwartz$M2.mm/10) 
# # 9.4 181.1 968.4
# estimate.from.uppers.corrected(data.schwartz$M2.mm/10)
# # 9.7, 187.1, 1000.4

# equations as functions lowers -----
#lowers: intercept 1.216, slope 3.010, %SEE .274
#m/2 length model
#log10(mass.kg) = log10(length.cm) * 3.010 + 1.216

estimate.from.lowers<-function(length.cm){
  log.mass<-log10(length.cm) * 3.010 + 1.216
  return(log.mass)
}

#check
#log mass estimate
# estimate.from.lowers(data.schwartz$m2.mm/10) 
# #detransformed estimate without correction
# 10^estimate.from.lowers(data.schwartz$m2.mm/10) 
# # 9.4 153.4 968.8

#detransformed, corrected estimate should be:
#9.2 155.4 1001.7
#check from Schwartz et al. 1995 Table 2
# (log(1001.7/988.8) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.069
# (log(155.4/153.4) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.069
# (log(9.4/9.2) * 2 / log(10)^2) %>% sqrt() #an SEE of 0.09

estimate.from.lowers.corrected<-function(length.cm){
  # SEE<-0.274 #from Janis 1990 #I can't get this number to work.
  #working backwards from raw/corrected detransformed masses instead gives
  SEE<-0.069
  cf<-(log(10)^2/2 * log10(exp(1)) * SEE^2) #Correction factor, Table 2, Smith 1993
  log.mass<-log10(length.cm) *  3.010 + 1.216
  corrected.mass<-10^log.mass * 10^cf
  return(corrected.mass)
}

estimate.from.lowers.corrected(data.schwartz$m2.mm/10) 

# Read in Topernawi Data, Format ------
#read in the raw data
data.raw<-read.csv(raw.measurements)
#filter down to upper and lower m2s of hyracoids with lengths
colnames(data.raw)

data.filter<-data.raw %>% 
  filter(Locus %in% c("M2","m2"))
  
data.filter$Length.MD.mm <- as.numeric(data.filter$L)
data.filter<-data.filter %>% filter(!is.na(Length.MD.mm))

# Apply Models ------
data.filter$estimated.log.mass<-0

#m/2 length model
#log10(mass.kg) = log10(length.mm) * 3.01 - 1.79

lowers<-which(data.filter$Locus=="m2")

data.filter$estimated.log.mass[lowers]<-estimate.from.lowers(data.filter$Length.MD.mm[lowers]/10)

uppers<-which(data.filter$Locus=="M2")
data.filter$estimated.log.mass[uppers]<-estimate.from.uppers(data.filter$Length.MD.mm[uppers]/10)

#correction:
data.filter$mass.kg<-10^(data.filter$estimated.log.mass)
data.filter$mass.kg.corrected<-0
data.filter$mass.kg.corrected[lowers]<-estimate.from.lowers.corrected(data.filter$Length.MD.mm[lowers]/10)

#uppers
data.filter$mass.kg.corrected[uppers]<-estimate.from.uppers.corrected(data.filter$Length.MD.mm[uppers]/10)

cbind(data.filter$mass.kg.corrected,data.filter$Species)

# # Graph -----
# setwd(working.directory)
ggplot(data = data.filter, aes(x = Species,y=mass.kg.corrected)) +
  geom_point(aes(shape = Locus)) + theme_minimal() + 
  xlab(NULL) + ylab("Corrected Mass (kg)") +labs(shape="Tooth \nPosition") +
  theme(axis.text.x = element_text(angle=315,hjust=-0.05),legend.position="top")
  
ggsave("mass_plot.pdf",width=3.539,height=4,units="in")

#error bars -------
#%SEE: if  it's 66%, then 68% (1 SD) of values expected to fall within +/- 66% of predicted value
data.averaged<-data.filter %>% group_by(Species, Locus) %>%
  mutate(mean.mass=mean(mass.kg.corrected))
data.averaged$P.SEE<-0.476 #value for upper molars, Janis 1990 appendix
data.averaged$P.SEE[which(data.averaged$Locus=="m2")]<-0.274 #value for lower molars
#note: using average reconstructed masses and %SEEs gives ranges that are much smaller than
#the variation between individual specimens. The resulting error bars are absurdly narrow
#in comparison to between-specimen ranges.
#Therefore, error bars not plotted because not informative about uncertainty. 
#Individual points are a better individual estimate of uncertainty. 

# save a table
write.csv(data.averaged,"body_size_results.csv")

data.averaged %>% summarise(mean(mass.kg.corrected)) %>% write.csv("body_size_averages.csv")
# ggplot(data = data.averaged, aes(x = Species,y=mass.kg.corrected)) +
#   geom_point(aes(shape = Locus),size=3) +
#   # geom_errorbar(aes(ymin=(mean.mass-2*P.SEE), ymax=(mean.mass+2*P.SEE))) +
#   theme_minimal() + 
#   xlab(NULL) + ylab("Corrected Mass (kg)") +labs(shape="Tooth \nPosition") +
#   theme(text = element_text(size=25), axis.text.x = element_text(angle=315,hjust=-0.05)) 
# 
# data.averaged %>% mutate(min.error=mean.mass-2*P.SEE)
