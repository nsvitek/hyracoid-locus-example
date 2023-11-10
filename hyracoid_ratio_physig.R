#phylogenetic tree for Topernawi hyracoid tooth locus work
#base plot version (cannot auto-get node values)
library(ape) #used throughouts
library(phytools) #I think maybe also used?
library(readxl) #for tip date spreadsheet
library(RRphylo) #for applying dates to tree, resolving polytomies

library(picante) #for blomberg's k
library(geiger) #for pagel's lambda
library(caper) #for pagel's lambda

# set path, objects used no matter what ------
getwd()
# setwd("D:/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position")
setwd("C:/Users/nsvit/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position")

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

#map locus-ratio traits on tree ----
#read in relative length values from saved output from literature
ratios.raw<-read.csv("output_lowers/ratios_differ_literature.csv")

#create a variable to match tree tip name formatting
ratios.raw$tip.name<-paste(ratios.raw$genus,ratios.raw$species,sep="_")

keep.these.tips<-which(cooper.fixed$tip.label %in% ratios.raw$tip.name)
drop.these.tips<-cooper.fixed$tip.label[-keep.these.tips]
cooper.trim<-drop.tip(cooper.fixed,drop.these.tips)
cooper.trim<-drop.tip(cooper.trim,c("Thyrohyrax_litholagus","Thyrohyrax_pygmaeus",
                                    "Thyrohyrax_meyeri",
                                    "Thyrohyrax_domorictus","Prohyrax_hendeyi",
                                   "Procavia_capensis"))
cooper.trim$node.label<-NA

plot(cooper.trim)
?drop.clade()

#pull the m2 ratio into a separate formatted object, following tutorials
m2.ratio<-ratios.raw$m2
names(m2.ratio)<-ratios.raw$tip.name
m2.ratio<-m2.ratio[complete.cases(m2.ratio)] #inferred names causing issues with caper downstream


#pull the m3 ratio into a separate formatted object, following tutorials
m3.ratio<-ratios.raw$m3
names(m3.ratio)<-ratios.raw$tip.name
m3.ratio<-m3.ratio[complete.cases(m3.ratio)]

#Blomberg et al's K
Kcalc(m2.ratio[cooper.fixed$tip.label], cooper.fixed)
phylosignal(m2.ratio[cooper.fixed$tip.label], cooper.fixed, reps = 1000)

Kcalc(m2.ratio[cooper.trim$tip.label], cooper.trim)
phylosignal(m2.ratio[cooper.trim$tip.label], cooper.trim, reps = 1000)

#Pagel's lambda, geiger
lambda_GL_ML<-fitContinuous(phy= cooper.trim, dat = m2.ratio,
                            model = "lambda")

#comparison to what you'd get if lambda were 1
lambda_GL_L1<-fitContinuous(phy = cooper.trim, dat = m2.ratio,
                            model = "BM")
LLR1 <- -2*(lambda_GL_L1$opt$lnL - lambda_GL_ML$opt$lnL)
pchisq(LLR1, df=1,lower.tail = FALSE)# p =1


#suggestion that OU model better fits conservatism #https://onlinelibrary.wiley.com/doi/10.1111/j.1420-9101.2010.02144.x
lambda_GL_L2<-fitContinuous(phy = cooper.trim, dat = m2.ratio,
                            model = "OU")
LLR2<- -2*(lambda_GL_L2$opt$lnL - lambda_GL_L1$opt$lnL)
#can't tell difference between BM and OU model?
pchisq(LLR2, df=1,lower.tail = FALSE)# p =1

#Pagel's lambda, caper
for.caper<-comparative.data(phy = cooper.trim, data = ratios.raw, 
                            names.col = tip.name,
                          vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)


model.lambda<-pgls(m2 ~ 1, data = for.caper, lambda='ML')

summary(model.lambda)

#https://wiki.duke.edu/pages/viewpage.action?pageId=131172281
#now with the m3

Kcalc(m3.ratio[cooper.fixed$tip.label], cooper.fixed)
phylosignal(m3.ratio[cooper.fixed$tip.label], cooper.fixed, reps = 1000)

Kcalc(m3.ratio[cooper.trim$tip.label], cooper.trim)
phylosignal(m3.ratio[cooper.trim$tip.label], cooper.trim, reps = 1000)

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