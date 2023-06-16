#in this case, measures of interest
data.comp.filter$rel.widths<-as.numeric(data.comp.filter$width.para)/as.numeric(data.comp.filter$width.meta)
data.comp.filter$rel.meta<-data.comp.filter$width.meta/data.comp.filter$length
data.comp.filter$proportion<-data.comp.filter$width.para/data.comp.filter$length

#standardize between the two sets
data.avg$genus<-"Procavia"
data.avg$species<-"capensis"
data.comp.filter$length<-as.numeric(data.comp.filter$length)

#merge the two datasets
data2merge<-data.avg %>% dplyr::select(Sample, Position, length, 
                                       rel.widths, rel.meta, proportion,
                                       genus, species)
data.comp.filter$Sample<-paste(data.comp.filter$coll,data.comp.filter$cat_num,sep="_")

data.c2merge<-data.comp.filter %>% 
  dplyr::select(Sample, Position, length,rel.widths, 
                rel.meta, proportion,genus, species) %>%
  melt(id=c("Sample","Position","genus","species")) %>%
  group_by(Sample,Position,variable,genus,species) %>% summarise(avg = mean(value)) %>% 
  dcast(Sample + Position + genus + species ~ variable)

data.jaws<-rbind(data2merge,data.c2merge)

#absolute length is not going to work well with statistics. 
#what I'm really interested in is the m/1 being shorter than m/2, etc.
#so translate to relative length.
#doing this the slow, painful way.
data.jaws$m1.length<-NA
for (i in 1:nrow(data.jaws)){
  wanted.value<-data.jaws$length[which(data.jaws$Sample==data.jaws$Sample[i] & 
                                         data.jaws$Position=="m1")]
  if (length(wanted.value)==0){ next }
  data.jaws$m1.length[i]<-wanted.value
}
data.jaws$rel.length<-data.jaws$length/data.jaws$m1.length

# check
# cbind(data.jaws$Sample, data.jaws$rel.length)[data.jaws$Position=="m1",]

linear2plot2<-data.jaws %>% 
  select(species,Position,rel.length,rel.widths,rel.meta, proportion) %>%
  melt(id=c("species","Position"))

#Make facet label names
measure.labels<-c("paraloph : length","length : m1 length","metloph : length","paraloph : metaloph")
names(measure.labels)<-c("proportion","rel.length","rel.meta","rel.widths")

