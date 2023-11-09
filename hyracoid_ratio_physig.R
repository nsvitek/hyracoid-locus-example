#phylogenetic tree for Topernawi hyracoid tooth locus work
#base plot version (cannot auto-get node values)
library(ape)
library(phytools)
library(readxl) #for tip date spreadsheet
library(RRphylo) #for dates


# set path, objects used no matter what ------
getwd()
setwd("D:/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position")

cooper.tree<-read.newick("hyracoid_cooper.tre")
#alternative, if no file:
# cooper.newick<-"((Ocepeia_daouiensis,Phosphatherium_escuilliei),Eritherium_azzouzorum,((((((Afrohyrax_championi,(Antilohyrax_pectidens,((Titanohyrax_andrewsi,Titanohyrax_sp._nov.),Titanohyrax_angustidens))),((Saghatherium_antiquum,Saghatherium_bowni),Selenohyrax_chatrathi)),((((Procavia_capensis,Prohyrax_hendeyi),Thyrohyrax_domorictus),Thyrohyrax_meyeri),Thyrohyrax_pygmaeus)),(((((Bunohyrax_fajumensis,Pachyhyrax_crassidentatus),Bunohyrax_major),(Geniohyus_diphycus,Geniohyus_mirus)),(Megalohyrax_eocaenus,Megalohyrax_sp._nov.)),Thyrohyrax_litholagus)),(Dimaitherium_patnaiki,Namahyrax_corvus)),(Microhyrax_lavocati,Seggeurius_amourensis)));"
# cooper.tree<-read.tree(text=cooper.newick)
plot(cooper.tree)

#add branch length, computed.
cooper.tree<-compute.brlen(cooper.tree,method="Grafen")

#add tip dates to scale tree -------
dates.raw<-read_xlsx("../phylogeny/hyracoid-tip-dates.xlsx")
tipAges<-dates.raw$min
names(tipAges)<-dates.raw$tip.label
tipAges<-tipAges[which(names(tipAges) %in% cooper.tree$tip.label)]
nodeAges<-70.1
names(nodeAges)<-getMRCA(cooper.tree,c("Ocepeia_daouiensis","Procavia_capensis")) #root
cooper.scaled<-scaleTree(cooper.tree,node.ages=nodeAges)
cooper.scaled<-scaleTree(cooper.scaled,tip.ages=tipAges)

plot(cooper.scaled)
axisPhylo()

#also fix the polytomies ------
cooper.fixed<-fix.poly(cooper.scaled,type="resolve")
plot(cooper.fixed)

?fix.poly
#map locus-ratio traits on tree ----
# test for phylogenetic signal/conservatism ------

# #plotting ------
# 
# backward.tree<-ggplot(bayesian.tree.mkv,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
#   theme_tree2() + #hexpand(0.02,direction=1) + #scale_x_continuous(labels = abs)
#   geom_range(range='height_0.95HPD',alpha=0.4,size=2) +
#   geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2) #+
#   scale_color_continuous(low="gray",high="black") #+
#     #improve the color scheme
# time.scale<-seq(from=-70,to=0,by=10)
# revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
# ggsave("mkv-allcompat.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 