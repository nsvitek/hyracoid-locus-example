remodel<-100 #number of times to reassign specimens and build new models

#build objects to keep track of misidentified specimens
misclassify<-vector(mode="list",length=remodel)

#empty array to hold multiple confusion matrices
#Note (because I can't figure out how to set this info in dimnames): rows are predicted, 
#cols are reference
total.confusion<-array(dim=c(3,3,remodel))
dimnames(total.confusion)[1:2]<-list(c("m1","m2","m3"),c("m1","m2","m3"))

#empty array to hold model metrics
empty.col<-rep(NA,remodel)
total.eval<-data.frame(N_train = empty.col, N_test = empty.col, 
                       Accuracy = empty.col,Kappa = empty.col)

start.seed<-1001 #make results more repeatable by setting the seed
i<-1 #start counter

# create models, evaluation data --------------------
for(j in 1:remodel){
  set.seed(start.seed)
  
  #because teeth within specimens are not independent from another, need to randomly sample at the level of specimens
  train.specimen<-factor(model_data$Sample) %>% levels() %>% 
    as.data.frame(.) %>% sample_frac(0.75)
  test.specimen<-model_data$Sample[!(model_data$Sample %in% train.specimen$.)]
  
  #specify which images go in the training set
  inTrain<-which(model_data$Sample %in% train.specimen$.) 
  
  #make the 3 datasets, removing the specimen column
  training <- model_data[ inTrain,-1]
  testing  <- model_data[-inTrain,-1]
  
  total.eval$N_train[i]<-nrow(training)
  total.eval$N_test[i]<-nrow(testing)
  
  #Training using lda method
  ldaFit <- train(
    Position ~ .,
    data = training,
    method = "lda",
    trControl = trainControl(method="cv") #, #ctrl,
    #metric = performance_metric
  )
  
  ldaClasses <- predict(ldaFit, newdata = testing)
  ldaProbs <- predict(ldaFit, newdata = testing, type = "prob")
  ldaConfusion<-confusionMatrix(ldaClasses, factor(testing$Position))
  
  ldaTogether<-cbind(obs=factor(testing$Position),pred=ldaClasses,ldaProbs)
  
  #Results LDA
  total.confusion[,,i]<-ldaConfusion$table #the actual confusion matrix
  total.eval$Accuracy[i]<-ldaConfusion$overall[1]
  total.eval$Kappa[i]<-ldaConfusion$overall[2]
  
  #check which ones gives wrong predictions
  wrong_pred <- which(ldaClasses!=factor(testing$Position))
  test.specimen %in% model_data$Sample
  misclassify[[i]]<-test.specimen[wrong_pred]
  
  #jump up
  start.seed<-start.seed+i*100
  i<-j+1
}

write.csv(total.confusion,"LDA_confusion.csv")
write.csv(total.eval,"LDA_evaluation.csv")

sink("LDA_misclassify.txt")
print(misclassify)
sink()

# illustrate evaluation metrics --------
#create a faceted plot like Puschel et al. 2018: each facet a metric, three rows for each model,
#and a simple line plot: geom_pointrange()
all.eval.long<-melt(total.eval,id.vars=c("N_train","N_test"))
head(all.eval.long)

ggplot(data=all.eval.long, aes(x=value, y = variable)) + 
  geom_line() +
  stat_summary(fun.data = "mean_cl_normal") +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("LDA_eval.pdf", device = cairo_pdf, width = double.column.width, 
       height = double.column.width,units="in",dpi=600)

#Create a table for reporting mean and sd values
total.eval %>% summarise_all(list(mean,sd), na.rm=T) %>% write.csv("model_evaluation.csv")

# compare confusion matrices -------
mean.confusion<-apply(total.confusion, c(1,2), mean) 
write.csv(mean.confusion, "confusion_summary_mean.csv")
sd.confusion<-apply(total.confusion, c(1,2), sd) 
write.csv(sd.confusion, "confusion_summary_sd.csv")

confusion.table<-data.frame(predicted.class=factor(c("m1","m2","m3","m1","m2","m3","m1","m2","m3")),
                            reference.class=factor(c("m1","m1","m1","m2","m2","m2","m3","m3","m3")),
                            mean.confusion=as.vector(mean.confusion) %>% round(1),
                            sd.confusion=as.vector(sd.confusion) %>% round(1))

confusion.table <- confusion.table %>%
  mutate(Classification = ifelse(predicted.class == reference.class, "correct", "incorrect")) %>%
  group_by(reference.class) %>%
  mutate(prop = round(mean.confusion/sum(mean.confusion),2))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups 
#(see dplyr code above as well as original confusion matrix for comparison)
confuse.lda.raw<-ggplot(data = confusion.table, mapping = aes(x = reference.class,
                  y = predicted.class, fill = Classification)) +
  geom_tile(alpha = confusion.table$prop) +
  geom_text(aes(label = paste(mean.confusion,"\n(",sd.confusion,")",sep="")),
            vjust = .5, fontface  = "bold", alpha = .5) +
  scale_fill_manual(values = c(correct = scale_n[3], incorrect = scale_n[10])) +
  theme_minimal() + ylab("Predicted Locus") +
  theme(axis.title.x=element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="italic"),legend.position = "none",
        axis.text.y = element_text(face="italic", angle=90),
        axis.title.y=element_text(size=8))+
  ylim(rev(levels(confusion.table$reference.class)))

confuse.lda.prop<-ggplot(data = confusion.table, mapping = aes(x = reference.class,
                   y = predicted.class, fill = Classification)) +
  geom_tile(alpha = confusion.table$prop) +
  geom_text(aes(label = confusion.table$prop),
            vjust = .5, fontface  = "bold", alpha = .5) +
  scale_fill_manual(values = c(correct = scale_n[3], incorrect = scale_n[10])) +
  theme_minimal() + ylab("Predicted Locus") +
  theme(axis.title.x=element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="italic"),legend.position = "none",
        axis.text.y = element_text(face="italic", angle=90),
        axis.title.y=element_text(size=8))+
  ylim(rev(levels(confusion.table$reference.class)))

ggsave("confusion_tables_raw.pdf", confuse.lda.raw,
       device = cairo_pdf, width = single.column.width,
       height=single.column.width, units="in",dpi=600)

ggsave("confusion_tables_prop.pdf", confuse.lda.prop,
       device = cairo_pdf, width = single.column.width,
       height=single.column.width, units="in",dpi=600)

ggsave("confusion_legend.pdf", confuse.lda.prop + theme(legend.position = "right"), device = cairo_pdf, width = single.column.width,
       height=single.column.width, units="in",dpi=600)

# Examine misclassified images -----
#calculate baseline frequencies of species, not specimens 
count.species<-as.factor(data.jaws$species) %>% summary()
freq.species<-count.species/nrow(data.jaws)
plot.freq.spp<-data.frame(Frequency=freq.species,
                           Species=names(freq.species))

# #[was originally developed for specimens]
# count.specimens<-as.factor(data.jaws$Sample) %>% summary()
# freq.specimens<-count.specimens/nrow(data.jaws)
# 
# plot.freq.spec<-data.frame(Frequency=freq.specimens,
#                            Specimen=names(freq.specimens))
# lda.misclassify.all<-lda.freq[rep(row.names(lda.freq), lda.freq$Frequency), ]
# lda.freq.specimen<-lda.misclassify.all$Specimen %>% as.factor %>% summary()/nrow(lda.misclassify.all)
# plot.freq.spec$MF.lda<-0.0
# plot.freq.spec$MF.lda[match(names(lda.freq.specimen),rownames(plot.freq.spec))]<-lda.freq.specimen
# plot.freq.spec$RMF.lda<-plot.freq.spec$MF.lda/plot.freq.spec$Frequency

#LDA
lda.freq<-misclassify %>% unlist %>% as.factor() %>% summary() %>% as.data.frame
lda.freq$Specimen<-lda.freq %>% rownames
colnames(lda.freq)[1]<-"Frequency"
lda.freq$species<-NA
for (i in 1:nrow(lda.freq)){
  lda.freq$species[i]<-data.jaws$species[which(data.jaws$Sample==lda.freq$Specimen[i])[1]]
}

#reformat
lda.freq.spp<-lda.freq %>% group_by(species) %>% summarise(sum(Frequency)) %>%
  as.data.frame
row.names(lda.freq.spp)<-lda.freq.spp$species
lda.misclassify.all.spp<-lda.freq.spp[rep(lda.freq.spp$species,
                                          lda.freq.spp$`sum(Frequency)`),]
#calculate a per-specimen # of times misclassified per total instances of misclassified specimens
lda.freq.spp<-lda.misclassify.all.spp$species %>% as.factor %>% summary()/nrow(lda.misclassify.all.spp)

plot.freq.spp$MF.lda<-0.0 #start with assumption that Misclassification Frequency=0
#unless shown otherwise
plot.freq.spp$MF.lda[match(names(lda.freq.spp),rownames(plot.freq.spp))]<-lda.freq.spp
plot.freq.spp$RMF.lda<-plot.freq.spp$MF.lda/plot.freq.spp$Frequency

# plot misclassification frequencies -------
plot.lda.spec<-ggplot(plot.freq.spp, aes(x=Species, y=RMF.lda)) +
  geom_hline(yintercept=1,lty=2) +
  geom_segment( aes(x=Species, xend=Species, y=0, yend=RMF.lda), color=scale_n[8]) +
  geom_point(size=3,color=scale_n[10],alpha=1) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),axis.title.x=element_blank(),
    axis.ticks.y = element_blank(),axis.text.y=element_text(size=8))
ggsave("misclassification_freq.pdf", plot.lda.spec, device = cairo_pdf, width = double.column.width,
       height = 4,units="in",dpi=600)
