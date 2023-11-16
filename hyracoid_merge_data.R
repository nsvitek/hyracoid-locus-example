#####
# This script is sourced by hyracoid_base.R
#####
#MODERN PROCAVIA FORMATTING -----------
#get measurement data formatted for plotting (in a separate, arcade-specific script)
source(paste(locateScripts,linear_format_script,sep="/"))

#COMPARATIVE FAYUM MEASUREMENTS FORMATTING ------------
#remove undescribed species and rejected measurements
data.comp.filter<-data.comparative %>% 
  filter(c(cat_num!="12471"&cat_num!="18667"&cat_num!="20509"))

#reformat data to make averaging easier
data.comp.filter.long<-melt(data.comp.filter, id=c("coll","cat_num","genus","species", "Side","Position"))

#fix quick data type issue preventing math
data.comp.filter.long$value<-as.numeric(data.comp.filter.long$value)

#average left and right sides of same specimen
asymmetry<-data.comp.filter.long %>% group_by(coll,cat_num,genus,species,Position,variable) %>% 
  summarise(avg = mean(value))
data.comp.avg<- dcast(asymmetry, coll + cat_num + genus+ species + Position ~ variable)

#get measurement data formatted for plotting (in a separate, arcade-specific script)
source(paste(locateScripts,comparative_format_script,sep="/"))

#COMPARATIVE LITERATURE ----
source(paste(locateScripts,literature_format_script,sep="/"))

if(arcade=="uppers"){ #minor formatting change to make late code work
  data.lit$Position<-tolower(data.lit$Position)
}
#merge all datasets ----------
data.jaws.all<-rbind(select(data2merge,Sample,Position,length,genus,species), 
                 select(data.c2merge,Sample,Position,length,genus,species),
                 select(data.lit,Sample,Position,length,genus,species))

# calculate length relative to m1 -------
#absolute length is not going to work well with statistics. 
#what I'm really interested in is the m/1 being shorter than m/2, etc.
#so translate to relative length.
#doing this the slow, painful way.
data.jaws.all$m1.length<-NA
for (i in 1:nrow(data.jaws.all)){
  wanted.value<-mean(data.jaws.all$length[which(data.jaws.all$species==data.jaws.all$species[i] & 
                                         data.jaws.all$Position=="m1")])
  if (length(wanted.value)==0){ next }
  data.jaws.all$m1.length[i]<-wanted.value
}
data.jaws.all$rel.length<-data.jaws.all$length/data.jaws.all$m1.length

#ditch NA's
data.jaws.all<-data.jaws.all[!is.na(data.jaws.all$rel.length),]
#also ditch NaN's
data.jaws.all<-data.jaws.all[!is.nan(data.jaws.all$rel.length),]




