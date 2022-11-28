#in this case, measures of interest are length and relative width. 
data.comp.filter$rel.widths<-as.numeric(data.comp.filter$trigonid.width)/as.numeric(data.comp.filter$talonid.width)
#standardize between the two sets
data.avg$genus<-"Procavia"
data.avg$species<-"capensis"
data.comp.filter$length<-as.numeric(data.comp.filter$length)

#merge the two datasets
data2merge<-data.avg %>% dplyr::select(Sample, Position, length, rel.widths, genus, species)
data.comp.filter$Sample<-paste(data.comp.filter$coll,data.comp.filter$cat_num,sep="_")
data.c2merge<-data.comp.filter %>% 
  dplyr::select(Sample, Position, length,rel.widths,genus, species) %>%
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

linear2plot2<-data.jaws %>% 
  select(species,Position,rel.length,rel.widths) %>%
  melt(id=c("species","Position"))

#Make facet label names
measure.labels<-c("length : M1 length", "trigonid width : talonid width")
names(measure.labels)<-c("rel.length","rel.widths")

#statistics: ANOVA ------
sink("ANOVA_results.txt")

print("Trigonid Width : Talonid Width")
aov(as.numeric(factor(Position)) ~ rel.widths, data = data.jaws) %>% summary

print("Length : M1 Length")
aov(as.numeric(factor(Position)) ~ rel.length, data = data.jaws) %>% summary

print("Test for Interaction Between Ratios")
aov(as.numeric(factor(Position)) ~ rel.widths * rel.length, data = data.jaws) %>% summary

sink()

#make data object for LDA
model_data<-data.jaws %>% dplyr::select(Sample, Position, rel.length, rel.widths,) %>% 
  na.omit()