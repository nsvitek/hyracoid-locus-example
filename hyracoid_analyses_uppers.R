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

overlap.table<-cbind(mean.ratios,p.vals.m1.rw,
                     p.vals.m2.rw,
                     p.vals.m1.l,
                     p.vals.m2.l)

write.csv(overlap.table,"ratios_differ.csv")
