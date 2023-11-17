# REPORTING: Sample Table ----
n.to.write<-count(data.jaws.all,genus,species,Position) %>% 
  dcast(genus + species ~ Position, value.var="n")
 
write.csv(n.to.write,"hyracoid_sampling_revised.csv")

#keep a dataset suitable for assessing standard deviation ------
#work with data.jaws & n.to.write
well.sampled.species<-n.to.write$species[which(n.to.write[,3]>=8)]

data.large<-data.jaws.all[which(data.jaws.all$species %in% well.sampled.species),]

stdevs.broad<-data.large %>% group_by(species,Position) %>% 
  summarise(stdev=sd(rel.length))
write.csv(stdevs.broad,"relative_length_standard_deviations.csv")

# ANALYSES ---------
# visualize data: scatterplot --------
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
# model variation --------
#use Procavia capensis (highest sample size among species) to model standard deviation of traits
stdevs<-linear2plot %>% group_by(Position, variable) %>% 
  summarise(stdev=sd(value, na.rm=TRUE))
write.csv(stdevs,"ratio_standard_deviations.csv")

# ARCADE-SPECIFIC: lowers-------
if(arcade=="lowers"){
  source(paste(locateScripts,"hyracoid_analyses_lowers.R",sep="/"))
}

if(arcade=="uppers"){
  source(paste(locateScripts,"hyracoid_analyses_uppers.R",sep="/"))
}

# ARCADE-GENERAL: include literature measurements ------
sink("ANOVA_results_literature.txt")

print("Length : M1 Length")
aov(as.numeric(factor(Position)) ~ rel.length * species, data = data.jaws.all) %>% summary

sink()

mean.ratios<-data.jaws.all %>% group_by(Position,genus,species) %>% 
  summarise(l = round(mean(rel.length, na.rm=TRUE),3)) %>% 
  dcast(genus + species ~ Position)

#is there directionality down tooth row? ------
directionality.table<-cbind(mean.ratios$m1 < mean.ratios$m2,
                            mean.ratios$m2 < mean.ratios$m3)

write.csv(cbind(mean.ratios,directionality.table),"ratios_increase_down_row_literature.csv")

#are successive mean values outside of prior locus 2*sd? -------

m1.shorter<-((mean.ratios$m1-mean.ratios$m2)/stdevs$stdev[6])  %>% pnorm() %>% round(4)
m2.shorter<-((mean.ratios$m2-mean.ratios$m3)/stdevs$stdev[10])  %>% pnorm() %>% round(4)

overlap.table<-cbind(mean.ratios,m1.shorter,m2.shorter)

#add correction for multiple tests:
overlap.table$m1.bh<-m1.shorter %>% p.adjust(method = "BH")
overlap.table$m2.bh<-m2.shorter %>% p.adjust(method = "BH")

write.csv(overlap.table,"ratios_differ_literature.csv")
