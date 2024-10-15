library(dplyr) #organize
library(ggplot2) #plot
library(readxl) #read  data
library(reshape2)
library(ggthemes)

locateData<-"D:/Dropbox/Documents/research/Turkana/hyracoidea"
data.topernawi<-read.csv(paste(locateData,"Topernawi_hyracoid_published_tooth_sizes.csv",sep="/"))

locus.palette<-c("#EECC66","#997700","#000000")

# Basic formatting transformations of strings and integers ------

#fix some numeric/text designations (due to typed NA's?)
data.topernawi$width<-as.numeric(data.topernawi$W_M)
data.topernawi$width.distal<-as.numeric(data.topernawi$W_D)
data.topernawi$length<-as.numeric(data.topernawi$L)

# format TOPERNAWI measurements ---------
#choose filters
hyracoid.ids<-c("Nengohyrax josephi","Abdahyrax ekuwomi","Geniohyus ewoi",
                "Thyrohyrax lokutani","Thyrohyrax ekaii")
hyracoid.molars<-c("m1","m2","m3","M1","M2","M3")

#apply filters to Topernawi catalog
topernawi.hyracoid<-data.topernawi %>% filter(Locus %in% hyracoid.molars)

#figure of ratios for tooth IDs ------
topernawi.hyracoid$relative.width<-topernawi.hyracoid$width/topernawi.hyracoid$width.distal
topernawi.hyracoid$arcade<-"upper"
topernawi.hyracoid$arcade[which(topernawi.hyracoid$Locus %in% c("m1","m2","m3"))]<-"lower"

ggplot(data=filter(topernawi.hyracoid,arcade=="lower"), 
       aes(x=relative.width, y=length, color=Locus,shape=Species)) + 
  geom_point()  +
  theme_minimal() +
  scale_color_manual(name="Species",values=locus.palette) +
  scale_shape_discrete(name="Tooth Position") +
  ylab("Length (mm)") + xlab("Relative Trigonid:Talonid Width") +
  theme(legend.text= element_text(face="italic")) +
  scale_y_continuous(trans='log10')
ggsave("FigX-LocusMeasures.pdf",width=5,height=4,units="in")

ggplot(data=filter(topernawi.hyracoid,arcade=="upper"), 
       aes(x=as.numeric(width), y=as.numeric(length), color=Locus,shape=Species)) + 
  geom_point()  +
  theme_minimal() +
  scale_color_manual(name="Species",values=locus.palette) +
  scale_shape_discrete(name="Tooth Position") +
  ylab("Length (mm)") + xlab("Width (mm)") +
  theme(legend.text= element_text(face="italic")) +
  scale_y_continuous(trans='log10')
ggsave("FigSIX-LocusMeasuresUppers.pdf",width=5,height=4,units="in")

#calculate RCA ratios ------
#slim down dataset
topernawi.hyracoid<-select(topernawi.hyracoid,-L,-W_M,-W_D)

#calculate some additional measurements
topernawi.hyracoid$area<-topernawi.hyracoid$length * topernawi.hyracoid$width
topernawi.hyracoid$proportion<-topernawi.hyracoid$length / topernawi.hyracoid$width

#create a second format to facilitate relative crown area (RCA) calculations
# take averages by taxon and tooth locus
topernawi.avg<-topernawi.hyracoid %>% group_by(Species,Locus) %>%
  summarise(area = mean(area, na.rm=TRUE))

colnames(topernawi.avg)[which(colnames(topernawi.avg)=="Species")]<-"species"
colnames(topernawi.avg)[which(colnames(topernawi.avg)=="Locus")]<-"Position"
data.jaws.spp<-melt(topernawi.avg, id=c("species","Position")) %>% 
  dcast(species ~ Position + variable,value.var = "value",fun.aggregate=mean,na.rm=TRUE)

#calculate RCA values
data.jaws.spp$m2.m1<- data.jaws.spp$m2_area / data.jaws.spp$m1_area

data.jaws.spp$m1.m2<- 1/ data.jaws.spp$m2.m1

data.jaws.spp$m3.m1<- data.jaws.spp$m3_area / data.jaws.spp$m1_area

data.jaws.spp$m1.m3<- 1/ data.jaws.spp$m3.m1

data.jaws.spp$M2.M1<- data.jaws.spp$M2_area / data.jaws.spp$M1_area
data.jaws.spp$M3.M1<- data.jaws.spp$M3_area / data.jaws.spp$M1_area

# TAXONOMIC TRAITS -------
# #look at Cooper et al. 2014 matrix character 52
# ggplot(data=data.jaws.spp, aes(x=m1.m2)) + 
#   geom_histogram() +
#   geom_vline(aes(xintercept = 1), linetype = "dashed", color = "red") +
#   geom_vline(aes(xintercept = 0.65), linetype = "dashed", color = "red")
# 
# #similar character for m1:m3 as 52?
# ggplot(data=data.jaws.spp, aes(x=m1.m3)) + 
#   geom_histogram() +
#   geom_vline(aes(xintercept = 0.7), linetype = "dashed", color = "red")

#RCA ------
#reformat yet again to get upper and lower proportions in a single column.
uppers<-data.jaws.spp %>% select(species,M2.M1,M3.M1)
lowers<-data.jaws.spp %>% select(species,m2.m1,m3.m1)
names(uppers)<-c("species","m2.m1","m3.m1")
uppers$arcade<-"uppers"
lowers$arcade<-"lower"
data.rca<-rbind(uppers,lowers)
# write.csv(data.jaws.spp,"RCA_results.csv")

#which plots the same information in a different format
#make the two polygons of "forbidden space" in Kavanagh model
y.lims<-range(data.rca$m3.m1,na.rm=T)
x.lims<-range(data.rca$m2.m1,na.rm=T)
poly.right<-data.frame(id=rep("one",4),
                       x=c(1,1,y.lims[2]+0.1,y.lims[2]+0.1),
                       y=c(0,1,y.lims[2]+0.1,0))
poly.left<-data.frame(id=rep("two",4),
                      x=c(x.lims[1]-0.1,1,1,x.lims[1]-0.1),
                      y=c(x.lims[1]-0.1,1,y.lims[2]+0.1,y.lims[2]+0.1))

ggplot(data=data.rca, aes(x=m2.m1, y=m3.m1)) +   #,color=species
  # geom_polygon(data=poly.right,aes(x=x,y=y))+
  # geom_polygon(data=poly.left,aes(x=x,y=y))+
  geom_point(size=5,aes(shape=species,color=arcade)) + 
  scale_color_manual(name="Dental Arcade",values=locus.palette[c(1,3)]) +
  coord_cartesian(xlim = x.lims, ylim = y.lims) +
  theme_minimal() + theme(legend.position="top")
ggsave("Topernawi-hyracoid-RCA.pdf",width=3.539,height=4,units="in")