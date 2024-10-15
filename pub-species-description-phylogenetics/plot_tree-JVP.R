#phylogenetic trees for Topernawi hyracoid systematic work

# install.packages("BiocManager", repos = "https://cloud.r-project.org")
# install.packages("treeio")
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)

# set path, objects used no matter what ------
getwd()
setwd("D://Dropbox/Documents/research/Turkana/hyracoidea/phylogeny/hyracoid_mrbayes_v4/")

#Define which taxa should be highlighted on the tree (based on sampling)
topernawi.taxa<-c("Nengohyrax_josephi","Abdahyrax_philipi","Geniohyus_ewoi",
                  "Thyrohyrax_lokutani", "Thyrohyrax_ekaii") 


#try tidyverse summary tree option -----
#note that these MrBayes-based objects no longer works with base plot. 
bayesian.tree.mkv<-read.mrbayes("mkv/halfcompat/hyracoidea_mrbayes.nex.con.tre")
bayesian.tree.mk<-read.mrbayes("mk/halfcompat/hyracoidea_mrbayes.nex.con.tre")

# #then need to edit these trees so they are readable nexus-formatted files
# #these objects can be plotted in base plot
# mcc.mkv<-read.nexus("mkv/hyracoidea_mrbayes.MCC.tre")
# map.mkv<-read.nexus("mkv/hyracoidea_mrbayes.MAP.tre")
# mcc.mk<-read.nexus("mk/hyracoidea_mrbayes.MCC.tre")
# map.mk<-read.nexus("mk/hyracoidea_mrbayes.MAP.tre")

# #these two trees are identical. Start with the mkv tree as the baseline "tree"
# tree<-bayesian.tree.mkv

# edit majority rule trees for plotting -----
#need to substitute "_" for " "
bayesian.tree.mkv@phylo$tip.label<-gsub("(.*)_(.*)","\\1 \\2",bayesian.tree.mkv@phylo$tip.label)
bayesian.tree.mk@phylo$tip.label<-gsub("(.*)_(.*)","\\1 \\2",bayesian.tree.mk@phylo$tip.label)

#need to add BPP as node labels, improve color scheme
bayesian.tree.mkv@data$prob<-as.numeric(bayesian.tree.mkv@data$prob)
bayesian.tree.mk@data$prob<-as.numeric(bayesian.tree.mk@data$prob)

#need to fix timescale label. Solution from Guangchuang Yu
#https://github.com/YuLab-SMU/ggtree/issues/87

#tip ranges would be nice. Use the geom_range option below

backward.tree<-ggplot(bayesian.tree.mkv,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
  theme_tree2() + #hexpand(0.02,direction=1) + #scale_x_continuous(labels = abs)
  geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
  geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2) #+
  scale_color_continuous(low="gray",high="black") #+
    #improve the color scheme
time.scale<-seq(from=-70,to=0,by=10)
revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
ggsave("mkv-halfcompat.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 


#need to fix timescale label. Solution from Guangchuang Yu
#https://github.com/YuLab-SMU/ggtree/issues/87

#tip ranges would be nice. Use the geom_range option below

backward.tree<-ggplot(bayesian.tree.mk,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
  theme_tree2() + #hexpand(0.02,direction=1) + #scale_x_continuous(labels = abs)
  geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
  geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2) #+
scale_color_continuous(low="gray",high="black") #+
#improve the color scheme
time.scale<-seq(from=-70,to=0,by=10)
revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
ggsave("mk-halfcompat.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 


#there is also a ggdensitree() option
# make plot of TNT tree --------
library(ape)
parsimony.tree<-ape::read.tree(text="(Ocepeia_daouiensis ,(Phosphatherium_escuilliei ,(Eritherium_azzouzorum ,(Seggeurius_amourensis ,(Microhyrax_lavocati ,(Namahyrax_corvus ,(Dimaitherium_patnaiki ,((Thyrohyrax_litholagus ,((((Afrohyrax_championi ,(Titanohyrax_angustidens ,(Antilohyrax_pectidens ,(Titanohyrax_andrewsi ,Titanohyrax_sp._nov. )))),(Geniohyus_ewoi ,((Bunohyrax_fajumensis ,Pachyhyrax_crassidentatus ),(((Bunohyrax_major ,(Geniohyus_diphycus ,Geniohyus_mirus )),(Megalohyrax_eocaenus ,Megalohyrax_sp._nov. )),(Nengohyrax_josephi ,Abdahyrax_philipi ))))),(Thyrohyrax_lokutani ,Thyrohyrax_ekaii )),(Thyrohyrax_meyeri ,(Thyrohyrax_pygmaeus ,(Thyrohyrax_domorictus ,(Procavia_capensis ,Prohyrax_hendeyi )))))),(Saghatherium_bowni ,(Saghatherium_antiquum ,Selenohyrax_chatrathi ))))))))));
")
jpeg("SI_Fig1_TNT_Tree.jpg",width=6,height=8,units="in",res=300)
par(mar=c(0.2,0.2,0.2,0.2))
plot(parsimony.tree, label.offset=0.3, edge.width =2,cex=1)
dev.off()
