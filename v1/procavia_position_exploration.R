#Script for exploratory analysis of high-resolution morphometric data
#Using an example of hyracoid lower molars
#Relies on scripts published as part of Vitek et al. 2017, Ecology & Evolution
#(Github Repository: observer-free-morphotype-characterization)
#Analyzes output from an auto3dgm alignment produced by PuenteAlignment-Slurm code
#(see Github repository PuenteAlignment-Slurm for more details on alignment process)
###############################################################################
locateScripts <- "C://scripts/observer-free-morphotype-characterization/"

locateData <- "C:/Users/nsvit/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position"
# locateData <- "D:/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position"

shape.folder<-"procavia-65k-1024-211119-qc" #lower teeth
# shape.folder<-"procavia-u56k-1024-220816" #upper teeth

metadata.file<-"modern-procavia-metadata.xlsx"
# metadata.file<-"modern-procavia-metadata-uppers.xlsx"

# Load Dependencies ------------------------------------------------------------------
setwd(locateScripts)
library(dplyr) #makes the %>% pipes work in functions. 
library(geomorph) #for reading in shape data
library(readxl) #for reading in metadata
library(scales) #for rescale for colors, not cutting anymore
library(khroma) #for Paul Tol color schemes
library(ggplot2)
source("rearrange3ddat.R")
source("flipPC.R")
source("anderson.R")
source("figpiece3d.R")
source("find_repeatablePCs.R")
source("PCheat.R")
source("sensitivity_utils.R")

# Figure Settings --------------------------------------------------------------
#Note that PeerJ only has one column, regardless of figure width
page.width<-6.5 #width in inches for PeerJ figures. 

#Set color scheme for PC plots using Paul Tol color scheme
color_scheme<-colour("high contrast")(3)

#set up a palette for heat map
rainbow1<-colour("sunset")

#set up a palette for tooth loci
scale_locus<-colour("discrete rainbow")(3)

# Read in data -----------------------------------------------------------------
setwd(locateData) #this setup assumes that your metadata and auto3dgm output folder are in the same place
taxa<-read_excel(metadata.file) #read in metadata

#navigate to data, then read it in and make it something that other functions can use
setwd(shape.folder)
molars<-read.morphologika("morphologika_unscaled_high.txt")  %>% 
  preprocess(.) 

#link datasets. There are two extra specimens in alignment not studied in the end.
#need to remove those extra shapes so that shape data and metadata line up
choose.match<-substr(dimnames(molars$scaled)[[3]],1,17) %in% substr(taxa$Filename,1,17)
length(dimnames(molars$scaled)[[3]])

molars$cs<-molars$centroid[choose.match]
molars$m2d<-molars$m2d[choose.match,]
molars$n<-nrow(molars$m2d)
molars$scaled<-molars$scaled[,,choose.match]

# Take Left-Right Averages -----------------------------------------------------
#method needs to work with the shape data as well.
#approach is to build an index of which rows in taxa & molars$m2d are left-right pairs
#doing this the slow, painful way.
taxa$SpecimenLocus<-paste(taxa$Sample,taxa$Position,sep="_") %>% factor

#make a placeholder matrix of correct dimensions
molars.avg<-matrix(data=NA, nrow=length(levels(taxa$SpecimenLocus)),
                   ncol=ncol(molars$m2d))
#just a placeholder matrix of correct dimensions, won't need the SpecimenLocus column downstream
taxa.avg<-taxa[1:length(levels(taxa$SpecimenLocus)),] %>% select(-c(Side,SpecimenLocus,Filename)) 

#average left and right sides of same specimen
for (i in 1:length(levels(taxa$SpecimenLocus))){
  desired.rows<-which(taxa$SpecimenLocus==levels(taxa$SpecimenLocus)[i])
  if (length(desired.rows) > 1){
    molars.avg[i,]<-colMeans(molars$m2d[desired.rows,])
    } else {
      molars.avg[i,]<-molars$m2d[i,]
    } #in case there is only 1 side, no need to take avg
  
  #can't take average of first two columns, sample and position (characters)
  taxa.avg[i,c(1:2)]<-taxa[desired.rows[1],c(1, 3)]
  if (length(desired.rows) > 1){
  taxa.avg[i,c(3:ncol(taxa.avg))]<-apply(taxa[desired.rows,5:(ncol(taxa)-1)],2,mean) %>% 
    as.list()
  } else {
    taxa.avg[i,c(3:ncol(taxa.avg))]<-taxa[desired.rows,5:(ncol(taxa)-1)]
  } #in case there is only 1 side, no need to take avg
}

# PCA --------------------------------------------------------------------------
#principal components analysis. m2d is just a 2D form of your 3D data
PCA<-prcomp(molars.avg,scale.=FALSE)

# plot initial PCA -------------------------------------------------------------
#get percent variation explained by PCs, for plotting
PC1.percent<-summary(PCA)$importance[2,1] %>% round(3) * 100 
PC2.percent<-summary(PCA)$importance[2,2] %>% round(3) * 100

#put PCs and metadata together for ease of ggplotting
PCA2Plot<-cbind(PCA$x,taxa.avg)

#plot
ggplot(PCA2Plot, aes(x = PC1, y= PC2, color = Position)) +
  geom_point(size=1) +
  scale_color_manual(values=unname(color_scheme),
                     labels=c(expression(M[1]),expression(M[2]),expression(M[3]))) +
  xlab(paste("PC 1 (",PC1.percent,"%)", sep="")) +
  ylab(paste("PC 2 (",PC2.percent,"%)", sep="")) +
  theme_classic() + theme(legend.position="top") #, text=element_text(family = "mono")
ggsave("PCA_with_wear.eps", width=page.width/2, height=page.width/2, units="in")

#calculate PC1 max and min shapes.
PC1.shapes<-pcdif(PCA,mshp(molars$m2d), pcs=1)

#find colors for each point corresponding to heat map of differences
differences_PC1<-shpdif(PC1.shapes$pc1$max,PC1.shapes$pc1$min,rainbow1,alter="none",outlier=FALSE)

#plot
open3d()
plot3d(PC1.shapes$pc1$max,axes=F,col=differences_PC1,size=10,xlab="",ylab="",zlab="")
writePLY("PC1_max.ply",format="ascii",pointRadius=0.005)
close3d()

open3d()
plot3d(PC1.shapes$pc1$min,axes=F,col=differences_PC1,size=10,xlab="",ylab="",zlab="")
writePLY("PC1_min.ply",format="ascii",pointRadius=0.005)
close3d()

#calculate PC2 max and min shapes
#calculate PC1 max and min shapes.
PC2.shapes<-pcdif(PCA,mshp(molars$m2d), pcs=2)

#find colors for each point corresponding to heat map of differences
differences_PC2<-shpdif(PC2.shapes$pc2$max,PC2.shapes$pc2$min,rainbow1,alter="none",outlier=FALSE)

#plot
open3d()
plot3d(PC2.shapes$pc2$max,axes=F,col=differences_PC2,size=10,xlab="",ylab="",zlab="")
writePLY("PC2_max.ply",format="ascii",pointRadius=0.005)
close3d()


open3d()
plot3d(PC2.shapes$pc2$min,axes=F,col=differences_PC2,size=10,xlab="",ylab="",zlab="")
writePLY("PC2_min.ply",format="ascii",pointRadius=0.005)
close3d()


# Procrustes ANOVA ------
gdf.hyrax<-geomorph.data.frame(coords = molars$m2d[choose.match,],
                               locus = taxa$Position) 

shape.anova<-procD.lm(coords ~ locus, data = gdf.hyrax, print.progress = FALSE)
summary(shape.anova)

#code copied from geomorph help examples to look at pairwise statistics between groups
PW <- pairwise(shape.anova, groups = gdf.hyrax$locus, covariate = NULL)
PW.sum<-summary(PW, test.type = "dist", confidence = 0.95, stat.table = TRUE)

#write results for reporting in manuscript
write.csv(summary(shape.anova)$table, file = "../output_statistics_all/ProcANOVA_overall.csv")
write.csv(PW.sum$summary.table, file = "../output_statistics_all/ProcANOVA_pairwise.csv")

# calculate mean shapes --------
mean.m1<-mshp(molars.avg[which(taxa.avg$Position=="m1"),])
mean.m2<-mshp(molars.avg[which(taxa.avg$Position=="m2"),])
mean.m3<-mshp(molars.avg[which(taxa.avg$Position=="m3"),])

# heat maps --------
differences_means12<-shpdif(mean.m1$meanshape,mean.m2$meanshape,rainbow1,alter="none",outlier=FALSE) #note $meanshape
differences_means23<-shpdif(mean.m2$meanshape,mean.m3$meanshape,rainbow1,alter="none",outlier=FALSE) #note $meanshape

open3d()
plot3d(mean.m1$meanshape,axes=F,col=differences_means12,size=10,xlab="",ylab="",zlab="")
writePLY("m1.vs.m2.ply",format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(mean.m2$meanshape,axes=F,col=differences_means12,size=10,xlab="",ylab="",zlab="")
writePLY("m2.vs.m1.ply",format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(mean.m2$meanshape,axes=F,col=differences_means23,size=10,xlab="",ylab="",zlab="")
writePLY("m2.vs.m3.ply",format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(mean.m3$meanshape,axes=F,col=differences_means23,size=10,xlab="",ylab="",zlab="")
writePLY("m3.vs.m2.ply",format="ascii",pointRadius=0.005)
rgl.close()