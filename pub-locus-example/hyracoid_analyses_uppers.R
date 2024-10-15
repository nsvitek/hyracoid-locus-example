# ANOVAs ----
#statistics: ANOVA ------
sink("ANOVA_results.txt")

print("Paraloph Width : Length")
aov(as.numeric(factor(Position)) ~ proportion * species, data = data.jaws) %>% summary

print("Metaloph Width : Length")
aov(as.numeric(factor(Position)) ~ rel.meta * species, data = data.jaws) %>% summary


print("Paraloph Width : Metaloph Width")
aov(as.numeric(factor(Position)) ~ rel.widths * species, data = data.jaws) %>% summary

print("Length : M1 Length")
aov(as.numeric(factor(Position)) ~ rel.length * species, data = data.jaws) %>% summary

print("Test for Interaction Between Ratios /n relative widths vs relative lengths")
aov(as.numeric(factor(Position)) ~ rel.widths * rel.length * species, data = data.jaws) %>% summary

print("relative widths vs. metalophid width")
aov(as.numeric(factor(Position)) ~ rel.widths * rel.meta * species, data = data.jaws) %>% summary

sink()

#is locus distinctiveness shared across species? ----

#table of mean per-locus values, 
mean.ratios<-data.jaws %>% group_by(Position,genus,species) %>% 
  summarise(rw = round(mean(rel.widths, na.rm=TRUE),3), 
            l = round(mean(rel.length, na.rm=TRUE),3),
            ml = round(mean(rel.meta, na.rm=TRUE),3)) %>% 
  melt(id.vars = c("Position","genus","species")) %>% 
  dcast(genus + species ~ Position + variable)

#is there directionality down tooth row? ------
directionality.table<-cbind(mean.ratios$m1_rw < mean.ratios$m2_rw,
                            mean.ratios$m2_rw < mean.ratios$m3_rw,
                            mean.ratios$m1_l < mean.ratios$m2_l,
                            mean.ratios$m2_l < mean.ratios$m3_l,
                            mean.ratios$m1_ml < mean.ratios$m2_ml,
                            mean.ratios$m2_ml < mean.ratios$m3_ml)
colnames(directionality.table)<-c("m1.para.narrower","m2.para.narrower",
                                  "m1.shorter","m2.shorter",
                                  "m1.meta.narrower","m2.meta.narrower")
write.csv(cbind(mean.ratios,directionality.table),"ratios_increase_down_row.csv")
#more mixed.  

#are successive mean values outside of prior locus 2*sd? -------
#given that we are working from lowers with hypothesis that directionality is consistent,
#do a one-tailed instead of a two-tailed test
#z transform then pnorm()
p.vals.m1.rw<-((mean.ratios$m1_rw-mean.ratios$m2_rw)/stdevs$stdev[8])  %>% pnorm() %>% round(4)
p.vals.m2.rw<-((mean.ratios$m2_rw-mean.ratios$m3_rw)/stdevs$stdev[12])  %>% pnorm() %>% round(4)
p.vals.m1.l<-((mean.ratios$m1_l-mean.ratios$m2_l)/stdevs$stdev[6])  %>% pnorm() %>% round(4)
p.vals.m2.l<-((mean.ratios$m2_l-mean.ratios$m3_l)/stdevs$stdev[10])  %>% pnorm() %>% round(4)
p.vals.m1.ml<-((mean.ratios$m1_ml-mean.ratios$m2_ml)/stdevs$stdev[7])  %>% pnorm() %>% round(4)
p.vals.m2.ml<-((mean.ratios$m2_ml-mean.ratios$m3_ml)/stdevs$stdev[11])  %>% pnorm() %>% round(4)

overlap.table<-cbind(mean.ratios,p.vals.m1.rw,
                     p.vals.m2.rw,
                     p.vals.m1.l,
                     p.vals.m2.l, p.vals.m1.ml, p.vals.m2.ml)

overlap.table$m1.rw.bh<-p.vals.m1.rw %>% p.adjust(method = "BH")
overlap.table$m2.rw.bh<-p.vals.m2.rw %>% p.adjust(method = "BH")
overlap.table$m1.l.bh<-p.vals.m1.l %>% p.adjust(method = "BH")
overlap.table$m2.l.bh<-p.vals.m2.l %>% p.adjust(method = "BH")
overlap.table$m1.ml.bh<-p.vals.m1.ml %>% p.adjust(method = "BH")
overlap.table$m2.ml.bh<-p.vals.m2.ml %>% p.adjust(method = "BH")

write.csv(overlap.table,"ratios_differ.csv")
