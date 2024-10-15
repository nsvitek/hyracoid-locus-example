#analyzing length and width measurements of hyracoid molars
#to compare:
#     (1) coefficients of variation of individual molar positions
#     (2) relative sizes of upper vs. lower tooth positions (ex: m1 vs M1)
#     (3) relative sizes of metameres (ex: m1 vs. m2)

#read in data -----------
library(dplyr) #organize
library(ggplot2) #plot
library(readxl) #read data
library(reshape2)
library(ggthemes)

locateData<-"D:/Dropbox/Documents/research/Turkana/hyracoidea/Vitek_Princehouse_2024_hyracoid_tooth_position"
data.procavia.hi<-read_excel(paste(locateData,"modern-procavia-uppers.xlsx",sep="/"))
data.procavia.lo<-read_excel(paste(locateData,"modern-procavia-lowers.xlsx",sep="/"))
data.comparative.hi<-read_excel(paste(locateData,"hyracoidea_fayum_locus_comparisons.xlsx",sep="/"), sheet = "Measurements-Upper")
data.comparative.lo<-read_excel(paste(locateData,"hyracoidea_fayum_locus_comparisons.xlsx",sep="/"), sheet = "Measurements-Lower")
data.topernawi<-read_excel(paste(locateData,"../../Topernawi_hyracoid_published_tooth_sizes.csv",sep="/"))
data.literature<-read_excel(paste(locateData,"published-tooth-sizes.xlsx",sep="/"))
data.genus<-read_excel(paste(locateData,"published-tooth-sizes.xlsx",sep="/"),sheet="GenusKey")

#plot settings
column.width<-3.539
column.width.twothirds<-4.8
double.columns.width<-7.166
locus.palette<-c("#EECC66","#997700","#000000")

#set the working directory
setwd("D:/Dropbox/Documents/research/Turkana/hyracoidea/manuscript-hyracoid-Topernawi-description")
# format PROCAVIA measurements -----------
#formatting to put together the two spreadsheets
data.procavia.hi$Position<-toupper(data.procavia.hi$Position)

colnames(data.procavia.hi)<-gsub("\\.para","",colnames(data.procavia.hi))
colnames(data.procavia.lo)<-gsub("trigonid\\.","",colnames(data.procavia.lo))

data.procavia.lo<-select(data.procavia.lo,-c(talonid.width.1,talonid.width.2,talonid.width.3))
data.procavia.hi<-select(data.procavia.hi,-c(width.meta.1,width.meta.2,width.meta.3))

data.procavia<-rbind(data.procavia.lo,data.procavia.hi)

#take average of triplicates to reduce error
data.procavia$length<-apply(cbind(data.procavia$length.1,data.procavia$length.2,data.procavia$length.3),1,mean)
data.procavia$width<-apply(cbind(data.procavia$width.1,data.procavia$width.2,data.procavia$width.3),1,mean)

data.procavia.slim<-select(data.procavia,Sample,Side,Position,length,width)

data.procavia.long<-melt(data.procavia.slim, id=c("Sample", "Side","Position"))

#average left and right sides of same specimen
asymmetry<-data.procavia.long %>% group_by(Sample,Position,variable) %>% summarise(value = mean(value))

#take species-level mean, standard deviation, and sample size
data.procavia.avg<-asymmetry %>% group_by(Position, variable) 

#add in taxonomic data
data.procavia.avg$species<-"capensis"

colnames(data.procavia.avg)[which(colnames(data.procavia.avg)=="variable")]<-"measure"
colnames(data.procavia.avg)[which(colnames(data.procavia.avg)=="position")]<-"Position"

# format FAYUM measurements ------------
#formatting to put together the two spreadsheets
data.comparative.hi$Position<-toupper(data.comparative.hi$Position)

data.comparative.lo<-select(data.comparative.lo,-c(talonid.width,trigonid.length))
data.comparative.hi<-select(data.comparative.hi,-c(width.meta))

colnames(data.comparative.hi)<-gsub("\\.para","",colnames(data.comparative.hi))
colnames(data.comparative.lo)<-gsub("trigonid\\.","",colnames(data.comparative.lo))

data.comparative<-rbind(data.comparative.hi,data.comparative.lo)
#remove undescribed species and rejected measurements
data.comp.filter<-filter(data.comparative,c(cat_num!="12471"&cat_num!="18667"&cat_num!="20509"))

#fix some numeric/text designations (due to typed NA's?)
data.comp.filter$length<-as.numeric(data.comp.filter$length)
data.comp.filter$width<-as.numeric(data.comp.filter$width)

#take averages by taxon and tooth locus, slim down dataset to necessary elements
data.comp.avg<-data.comp.filter %>% select(cat_num, species, Position,length,width) %>% 
  melt(id.vars = c("species","Position","cat_num")) %>%
  group_by(species,Position,variable)

colnames(data.comp.avg)[which(colnames(data.comp.avg)=="variable")]<-"measure"
colnames(data.comp.avg)[which(colnames(data.comp.avg)=="cat_num")]<-"Sample"


# format TOPERNAWI measurements ---------
#fix some numeric/text designations (due to typed NA's?)
data.topernawi$width<-as.numeric(data.topernawi$W_M)
data.topernawi$length<-as.numeric(data.topernawi$L)

#choose filters
hyracoid.ids<-c("Nengohyrax josephi", "Abdahyrax ekuwomi", "Geniohyus ewoi", "Thyrohyrax lokutani", "Thyrohyrax ekaii")
hyracoid.molars<-c("m1", "m2","m3","M1","M2","M3")

#apply filters to Topernawi catalog
topernawi.hyracoid<-data.topernawi %>% filter(Locus %in% hyracoid.molars)

#take averages by taxon and tooth locus, slim down dataset to necessary elements
topernawi.avg<-topernawi.hyracoid %>% 
  melt(id.vars = c("Specimen","Species","Locus")) %>%
  group_by(Species,Locus,variable)

colnames(topernawi.avg)[which(colnames(topernawi.avg)=="Species")]<-"species"
colnames(topernawi.avg)[which(colnames(topernawi.avg)=="Locus")]<-"Position"
colnames(topernawi.avg)[which(colnames(topernawi.avg)=="variable")]<-"measure"
colnames(topernawi.avg)[which(colnames(topernawi.avg)=="SpecimenNumber")]<-"Sample"

# format LITERATURE measurements ---------

data.lit.long<-select(data.literature,-c(publication,page,Notes)) %>%
  melt(id=c("institution", "number","species","side")) %>% filter(!is.na(value))
data.lit.long$measure<-"length"
data.lit.long$measure[grep(".* W$", data.lit.long$variable)]<-"width"
data.lit.long$position<-gsub("(.*) [L|W]$", "\\1", data.lit.long$variable)
data.lit.long$Sample<-paste(data.lit.long$institution,data.lit.long$number)

data.lit.avg.long<-data.lit.long %>% select(-c(variable,side,institution, number)) %>% 
  group_by(species, measure,position)
colnames(data.lit.avg.long)[which(colnames(data.lit.avg.long)=="position")]<-"Position"

# merge the datasets ----------
data.all<-rbind(data.procavia.avg, data.comp.avg, topernawi.avg, data.lit.avg.long)

data.all$genus<-data.genus$Genus[match(data.all$species,data.genus$Species)]

#remove genera == unclear
data.all<-data.all[-which(data.all$genus=="unclear"),]

#use Synonym key to rename species as needed
have.synonyms<-which(!is.na(data.genus$Synonym))
rename.guide<-data.genus$Species[have.synonyms]
need.renaming<-which(data.all$species %in% rename.guide)
new.names.index<-have.synonyms[match(data.all$species[need.renaming],rename.guide)]
new.names<-data.genus$Synonym[new.names.index]
data.all$species[need.renaming]<-new.names
data.all$name<-paste(data.all$genus, data.all$species)
for(i in 1:nrow(data.all)){
  if(grepl("NA .*",data.all$name[i])){data.all$name[i]<-data.all$species[i]}
}

colnames(data.all)[which(colnames(data.all)=="name")]<-"Species"

data.jaws<-data.all %>% 
  group_by(Species, measure,Position) %>% 
  summarise(mean = mean(value, na.rm=TRUE), sigma = sd(value, na.rm=TRUE), N = n()) %>%
  filter(Position %in% c("m1","m2","m3","M1","M2","M3")) 

data.jaws$species

# COEFFICIENT OF VARIATION ------
#calculate coefficient of variation
data.jaws$CV<-round(data.jaws$sigma / data.jaws$mean *100, 3)
data.jaws$Topernawi<-data.jaws$Species %in% hyracoid.ids

#plot data

ggplot(data = data.jaws, aes(x = Position, y = CV,alpha= Topernawi, 
                             color = Topernawi, shape = Topernawi)) + 
  geom_jitter() + 
  scale_alpha_manual(values=c(0.5,1)) +
  scale_color_manual(values=c("gray50","black")) +
  facet_wrap(ncol=2, "measure", scales="free_y", strip.position="top") +
  theme_minimal()
ggsave("context_CV.pdf",device = cairo_pdf, width = double.columns.width,
       height = column.width,units="in",dpi=600)

#report coefficients of variation in a table
data.jaws %>% select(-c(Topernawi,mean, sigma)) %>%
  filter(!is.na(CV)) %>%
  melt(id.vars=c("Position","measure","Species")) %>%
  dcast(Species + measure ~ Position + variable) %>%
  group_by(Species) %>% arrange(.by_group = TRUE)  %>% 
write.csv("SITableX-hyracoid-CV.csv", na="")

# UPPER VS. LOWER PROPORTIONS ------
#add upper/lower designation columns
data.jaws$Arcade<-"upper"
data.jaws$Arcade[which(data.jaws$Position %in% c("m1","m2","m3"))]<-"lower"
data.jaws$Locus<-"first"
data.jaws$Locus[which(grepl("[mM]2",data.jaws$Position))]<-"second"
data.jaws$Locus[which(grepl("[mM]3",data.jaws$Position))]<-"third"

#reformat for plotting
data.arcade<-data.jaws %>% select(-c(sigma,N,CV,Position)) %>%
  dcast(Species + measure + Topernawi + Locus ~ Arcade, value.var = "mean")
#plot
ggplot(data = data.arcade, aes(x = lower, y = upper,alpha= Topernawi, 
                               color = Topernawi, shape = Topernawi)) + 
  geom_point() + 
  scale_alpha_manual(values=c(0.5,1)) +
  scale_color_manual(values=c("gray50","black")) +
  facet_wrap(measure ~ Locus, nrow=2,ncol=3, scales="free", strip.position="top") +
  theme_minimal()
ggsave("context_arcade_proportions.pdf",device = cairo_pdf, width = double.columns.width,
       height = column.width*2,units="in",dpi=600)

# RCA ------
#reformat to facilitate relative crown area (RCA) calculations

data.rca.format<-data.jaws %>% select(-c(Position,sigma,N,CV)) %>%
  dcast(Species + Topernawi + Arcade + Locus ~ measure,value.var = "mean",na.rm=TRUE) %>%
  mutate(area = length * width) %>% select(-c(length,width)) %>%
  dcast(Species + Topernawi + Arcade ~ Locus,value.var = "area",na.rm=TRUE)
data.rca<-data.rca.format[complete.cases(data.rca.format),]

#calculate RCA values
data.rca$m2.m1 <- data.rca$second / data.rca$first
data.rca$m3.m1 <- data.rca$third / data.rca$first

#plot
ggplot(data = data.rca, aes(x = m2.m1, y = m3.m1,alpha= Topernawi, color = Topernawi,
                            shape=Topernawi)) + 
  geom_point() + 
  scale_alpha_manual(values=c(0.5,1)) +
  scale_color_manual(values=c("gray50","black")) +
  facet_wrap(~ Arcade, ncol=2, scales="free", strip.position="top") +
  theme_minimal()
ggsave("context_rca.pdf",device = cairo_pdf, width = double.columns.width,
       height = column.width,units="in",dpi=600)
