#COMPARATIVE LITERATURE ----
file.name<-paste(locateData,"published-tooth-sizes.xlsx",sep="/")
data.literature<-read_excel(file.name,sheet="LiteratureMeasurements")
data.genus<-read_excel(file.name,sheet="GenusKey")

#add in genus names -----
data.literature$genus<-data.genus$Genus[match(data.literature$species,data.genus$Species)]

#remove genera == unclear
data.literature<-data.literature[-which(data.literature$genus=="unclear"),]

#use Synonym key to rename species as needed
have.synonyms<-which(!is.na(data.genus$Synonym))
rename.guide<-data.genus$Species[have.synonyms]
need.renaming<-which(data.literature$species %in% rename.guide)
new.names.index<-have.synonyms[match(data.literature$species[need.renaming],rename.guide)]
new.names<-data.genus$Synonym[new.names.index]
data.literature$species[need.renaming]<-new.names

# format ------
data.lit.long<-dplyr::select(data.literature,-c(publication,page,Notes)) %>%
  melt(id=c("institution", "number","genus","species","side")) %>% filter(!is.na(value))
data.lit.long$measure<-"length"
data.lit.long$measure[grep(".* W$", data.lit.long$variable)]<-"width"
data.lit.long$position<-gsub("(.*) [L|W]$", "\\1", data.lit.long$variable)
data.lit.long$Sample<-paste(data.lit.long$institution,data.lit.long$number)


if(arcade=="lowers"){
    data.lit.filter<-filter(data.lit.long,position %in% c("m1","m2","m3"))
}
if(arcade=="uppers"){
  data.lit.filter<-filter(data.lit.long,position %in% c("M1","M2","M3"))
  
}

data.lit<-dcast(data.lit.filter,Sample + genus + species + position ~ measure,
                    fun.aggregate=mean,na.rm=TRUE)

colnames(data.lit)[which(colnames(data.lit)=="position")]<-"Position"

#may someday want species average values
# data.lit.avg.long<-data.lit %>% select(-variable) %>% group_by(genus,species,value, measure,position) %>% 
#   summarise(value = mean(value, na.rm=TRUE))