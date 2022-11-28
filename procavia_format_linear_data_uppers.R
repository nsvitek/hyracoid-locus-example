#MODERN PROCAVIA -----------
#take average of triplicates to reduce error --------
data.procavia$length<-apply(cbind(data.procavia$length.1,data.procavia$length.2,data.procavia$length.3),1,mean)
data.procavia$width.para<-apply(cbind(data.procavia$width.para.1,data.procavia$width.para.2,data.procavia$width.para.3),1,mean)
data.procavia$width.meta<-apply(cbind(data.procavia$width.meta.1,data.procavia$width.meta.2,data.procavia$width.meta.3),1,mean)

#trim dataset to necessary columns
data<-dplyr::select(data.procavia,Sample,Side,Position,length,width.para,width.meta)

#reformat data to make averaging easier
data.long<-melt(data, id=c("Sample", "Side","Position"))

#average left and right sides of same specimen
asymmetry<-data.long %>% group_by(Sample,Position,variable) %>% summarise(avg = mean(value))
data.avg<- dcast(asymmetry, Sample + Position ~ variable)

#calculate ratios ------
data.avg$rel.meta<-data.avg$width.meta/data.avg$length #compare to rel.widths
data.avg$rel.widths<-data.avg$width.para/data.avg$width.meta #compare to rel.tal.width
data.avg$proportion<-data.avg$width.para/data.avg$length #compare to relative length

#absolute length is not going to work well with statistics. 
#what I'm really interested in is the m/1 being shorter than m/2, etc.
#so translate to relative length.
#doing this the slow, painful way.
data.avg$m1.length<-NA
for (i in 1:nrow(data.avg)){
  wanted.value<-data.avg$length[which(data.avg$Sample==data.avg$Sample[i] &
                                        data.avg$Position=="m1")]
  if (length(wanted.value)==0){ next }
  data.avg$m1.length[i]<-wanted.value
}
data.avg$rel.length<-data.avg$length/data.avg$m1.length

#plot ratio data ------
#make a plot similar to Fig. 6 of Vitek and Chen 2022: grouped series of scatterplots + boxplots
#to give initial sense of utility of various metrics. 
linear2plot<-data.avg %>% 
  select(Position,proportion,rel.length,rel.meta,rel.widths) %>%
  melt(id="Position")

#Make facet label names
measure.labels<-c("paraloph : length","length : m1 length","metloph : length","paraloph : metaloph")
names(measure.labels)<-c("proportion","rel.length","rel.meta","rel.widths")
