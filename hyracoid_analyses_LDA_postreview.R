#Hyracoid locus code after revision: Request for LDA
library(gridExtra)

# Preparation ---------
#as of first draft, this code works after data are merged, 
#before the rest of the analyses are necessarily called in `hyracoid_analyses.R`
#Technically, most LDA intended for multivariate, most tests in paper are univariate, so work up
# Functions for univariate LDA visualization --------
plot.locus.overlap<-function(data.jaws, species.name){
  for.curves<-data.jaws %>% filter(species == species.name) %>% group_by(Position) %>% 
    summarise(Mean=mean(length), Variance = var(length))
  #Use mean and variance of ratios for 2 loci, following p. 208 of Bookstein 2018 textbook
  #modifying from here: https://stackoverflow.com/questions/74415495/plot-normal-distribution-when-only-mean-and-standar-deviation-exists-in-ggplot2
  made.curves<-for.curves %>% 
    mutate(min_range = min(Mean - 4 * sqrt(Variance)),
           max_range = max(Mean + 4 * sqrt(Variance))) %>%
    mutate(group = row_number()) %>%
    rowwise() %>%
    reframe(x = seq(min_range, max_range, length.out = 500),
            y = dnorm(x, Mean, sqrt(Variance)),
            is_2sd = x > Mean - 2 * sqrt(Variance) & 
              x < Mean + 2 * sqrt(Variance),
            group = group)
  
  range.1<-range(made.curves$x[which(made.curves$group==1 & made.curves$is_2sd==TRUE)])
  range.2<-range(made.curves$x[which(made.curves$group==2 & made.curves$is_2sd==TRUE)])
  range.3<-range(made.curves$x[which(made.curves$group==3 & made.curves$is_2sd==TRUE)])
  
  #make empty objects to later fill with an indicator for what values to color in
  made.curves$fill_1<-made.curves$fill_2a<-made.curves$fill_2b<-made.curves$fill_3<-FALSE
  
  #find the region of group 1 that's within 2sds of group 2
  made.curves$fill_1[which(made.curves$x[which(made.curves$group==1)] > min(range.2))]<-TRUE
  
  #find the region of group 2 that's within 2sds of groups 3
  made.curves$fill_2b[which(made.curves$group==2 & (made.curves$x > range.3[1]) &
                             made.curves$group==2 & (made.curves$x < range.3[2]) )]<-TRUE
  #same, but group 2 vs. group 1
  made.curves$fill_2a[which(made.curves$group==2 & (made.curves$x > range.1[1]) &
                             made.curves$group==2 & (made.curves$x < range.1[2]))]<-TRUE
  
  #find the region of group 3 that's within 2sds of group 2
  made.curves$fill_3[which(made.curves$group==3 & (made.curves$x < range.2[2]))]<-TRUE
  
  #make a plot that colors ranges of "hard to classify" values
  output.plot<-ggplot(made.curves,aes(x, y)) +
    geom_area(data = . %>% filter(fill_1), fill='red', alpha=0.2, position = 'identity') +
    geom_area(data = . %>% filter(fill_2a), fill='gold', alpha=0.2, position = 'identity') +
    geom_area(data = . %>% filter(fill_2b), fill='gold', alpha=0.2, position = 'identity') +
    geom_area(data = . %>% filter(fill_3), fill='blue', alpha=0.2, position = 'identity') +
    geom_line(aes(group = group)) +
    theme_minimal() +  #base_size = 16
    labs(x = "", y = "Proportional Frequency") + #x = "Length",
    scale_x_discrete(labels = NULL, breaks = NULL)
  return(output.plot)
}

#density plot function (for consistency across species)
plot.density.2d<-function(data.jaws, species.name){
  plot.me<-ggplot(filter(data.jaws, species==species.name),
                  aes(x = length, y = rel.widths, color = Position, shape = Position)) +
    geom_density_2d(alpha=0.2) +
    geom_point() + theme_minimal() + theme(legend.position="none") +
    labs(x = "Length", y = "Relative Widths")
  geom_point()
  return(plot.me)
}
# Procavia capensis ------- 
# well-sampled empirical: #Procavia has both ratios??
#for just size, Saghatherium antiquum, Prohyrax hendeyi, maybe Afrohyrax championi
#look at raw histogram to get a sense of the data
# ggplot(filter(data.jaws, species=="capensis"),
#        aes(x = length, fill = Position)) +
#   geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity',binwidth=0.25)
#note, in case you're thinking it: that large m2 and relatively large m1? That's UMZC_H4970J and 
#yes, that's absolutely an m1 and m2, they are ajacent to an erupting p4: https://www.morphosource.org/concern/parent/0000S5237/media/000022983
# data.jaws %>% filter(species == "capensis", Sample =="UMZC_H5081B") 

#plot length discrimination ability
capensis.curve.plot<-plot.locus.overlap(data.jaws, "capensis")

#now bring in the second dimension (discriminator)
capensis.density.plot<-plot.density.2d(data.jaws, "capensis")

# Saghatherium bowni -----
#plot length discrimination
bowni.curve.plot<-plot.locus.overlap(data.jaws, "bowni")

#now bring in the second dimension (discriminator)
bowni.density.plot<-plot.density.2d(data.jaws, "bowni")

# Thyrohyax domorictus ----
#plot length discrimination
domorictus.curve.plot<-plot.locus.overlap(data.jaws, "domorictus")

#now bring in the second dimension (discriminator)
domorictus.density.plot<-plot.density.2d(data.jaws, "domorictus")

# composite plot ------
composite.density<-grid.arrange(capensis.curve.plot,bowni.curve.plot, domorictus.curve.plot,
             capensis.density.plot,bowni.density.plot, domorictus.density.plot,
             ncol=3, nrow = 2, heights = c(2,3))
ggsave("visual_LDA.pdf",composite.density,device = cairo_pdf, width = double.column.width,
       height = double.column.width,units="in",dpi=600)

#LDA capensis ------
#note that using relative length vs. length makes a difference in accuracy, but
#relative length implies that you know that two teeth belong together in a jaw.
#Therefore, for LDA the fairer test of discriminatory ability uses absolute length, which allows for
#teeth to be isolated. 

#calculate the linear discriminant model
LinDisc<-lda(Position ~ length + rel.widths, data=filter(data.jaws, species == "capensis"),CV=TRUE)
#evaluate the predictive ability for the model on the data used to build it
#not ideal, but sample size doesn't allow for better option
conf.matrix<-table(filter(data.jaws, species == "capensis")$Position, LinDisc$class)
compare.lda<-data.jaws %>% filter(species == "capensis") %>% dplyr::select(Position) %>% 
  cbind(LinDisc$class)

#% accuracy
length(which(compare.lda$Position==compare.lda$`LinDisc$class`))/nrow(compare.lda)

#write confusion matrix and accuracy to files. 

# LDA bowni ------
LinDisc<-lda(Position ~ length + rel.widths, data=filter(data.jaws, species == "bowni"))

PredPos<-predict(LinDisc,filter(data.jaws, species == "bowni"))$class

compare.lda<-data.jaws %>% filter(species == "bowni") %>% dplyr::select(Position) %>% cbind(PredPos)
compare.lda$Correct<-compare.lda$Position==compare.lda$PredPos

#% correct classification
length(which(compare.lda$Position==compare.lda$PredPos))/nrow(compare.lda)

# LDA Thyrohyrax domorictus ----
#LDA
LinDisc<-lda(Position ~ length + rel.widths, data=filter(data.jaws, species == "bowni"))

PredPos<-predict(LinDisc,filter(data.jaws, species == "domorictus"))$class

compare.lda<-data.jaws %>% filter(species == "domorictus") %>% dplyr::select(Position) %>% cbind(PredPos)
compare.lda$Correct<-compare.lda$Position==compare.lda$PredPos

#% correct classification
length(which(compare.lda$Position==compare.lda$PredPos))/nrow(compare.lda)

# Meroehyrax kyongi as applied example ---------

data.test<-read_xlsx("../test_case_meroehyrax.xlsx") %>% 
  select(-c(institution, species, publication, page, Notes)) %>% melt(id.var="number", na.rm=TRUE)

data.test$Position<-gsub("(m[123]).*","\\1",data.test$variable)
data.test$measure<-gsub("m[123] (.*)","\\1",data.test$variable)

data.tf<-data.test %>% select(-variable) %>% dcast(formula = number + Position ~ measure)

data.tf$rel.widths<-data.tf$WM/data.tf$WD

ggplot(data.tf,  aes(x = L, y = rel.widths, color = Position)) +
  # geom_density_2d() +
  geom_point() + geom_text(aes(label = number),vjust=-1.)



#for uppers: Procavia capensis, Saghatherium bowni, Thyrohyrax meyeri -----