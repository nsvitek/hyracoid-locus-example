#####
# This script is sourced by hyracoid_base.R
#####
#MODERN PROCAVIA -----------
#get measurement data formatted for plotting (in a separate, arcade-specific script)
source(paste(locateScripts,linear_format_script,sep="/"))

#plot
ggplot(linear2plot, aes(x=Position, y = value)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(aes(color=as.factor(Position)))  + 
  # xlab(c(expression(M[1]),expression(M[2]),expression(M[3]))) 
  theme_minimal() + ylab("Ratio") +
  # facet_grid() + 
  facet_wrap(vars(variable), scales="free_y", strip.position="top",
             labeller = labeller(variable = measure.labels))+
  theme(axis.text.x = element_text(face="italic"),
        strip.text.y.right = element_text(angle=0))+
  scale_x_discrete(labels=c(expression(M[1]),expression(M[2]),expression(M[3]))) +
  scale_color_manual(name="Locus",values=scale_locus,
                     labels=c(expression(M[1]),expression(M[2]),expression(M[3])))

ggsave("linear_boxplots.pdf", device = cairo_pdf, width = double.column.width,
       height = double.column.width,units="in",dpi=600)

#COMPARATIVE FAYUM MEASUREMENTS ------------
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

##START HERE. WORK ON THE UPPERS SCRIPT
#get measurement data formatted for plotting (in a separate, arcade-specific script)
source(paste(locateScripts,comparative_format_script,sep="/"))

# plot larger sample scatterplot --------
ggplot(linear2plot2, aes(x=species, y = value, color=Position)) +
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitterdodge()) + 
  theme_minimal() + ylab("Ratio") + xlab("Species") +
  facet_wrap(nrow=2, vars(variable), scales="free_y", strip.position="top",
             labeller = labeller(variable = measure.labels))+
  theme(axis.text.x = element_text(face="italic"),
        strip.text.y.right = element_text(angle=0))+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_color_manual(name="Locus",values=scale_locus[c(2,3,1)],
                     labels=c(expression(M[1]),expression(M[2]),expression(M[3])))
ggsave("comparative_ratios.pdf",device = cairo_pdf, width = double.column.width,
       height = double.column.width,units="in",dpi=600)

# #Reporting: Sample Table ----
#pre-review version
# data.jaws %>% group_by(genus,species) %>% summarise(n = length(unique(Sample))) %>%
#   write.csv("hyracoid_sampling.csv")

#post-review version
n.to.write<-data.jaws %>% group_by(genus,species,Position) %>% 
  summarise(n = length(unique(Sample))) #%>% 
dcast(n.to.write, genus + species ~ Position, value.var="n") %>% 
  write.csv("hyracoid_sampling_revised.csv")

#statistics: LDA --------
#note, only ~20% of teeth are m3 in dataset
length(which(model_data$Position=="m3"))/nrow(model_data)
length(which(model_data$Position=="m2"))/nrow(model_data)

#Linear discriminant analyses held in another script
source(paste(locateScripts,"hyracoid_lda.R",sep="/"))
