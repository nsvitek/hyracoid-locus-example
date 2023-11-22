#phylogenetic tree for Topernawi hyracoid tooth locus work

#intended for use after data are merged, 
#before the rest of the analyses are necessarily called in `hyracoid_analyses.R`

library(ape) #used throughouts
library(phytools) #I think maybe also used?
library(readxl) #for tip date spreadsheet
library(RRphylo) #for applying dates to tree, resolving polytomies

library(picante) #for blomberg's k

# set path, objects used no matter what ------
getwd()
# setwd("D:/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position")
# setwd("C:/Users/nsvit/Dropbox/Documents/research/Turkana/hyracoidea/hyracoid_tooth_position")

cooper.tree<-read.newick("../hyracoid_cooper.tre")
#alternative, if no file:
# cooper.newick<-"((Ocepeia_daouiensis,Phosphatherium_escuilliei),Eritherium_azzouzorum,((((((Afrohyrax_championi,(Antilohyrax_pectidens,((Titanohyrax_andrewsi,Titanohyrax_sp._nov.),Titanohyrax_angustidens))),((Saghatherium_antiquum,Saghatherium_bowni),Selenohyrax_chatrathi)),((((Procavia_capensis,Prohyrax_hendeyi),Thyrohyrax_domorictus),Thyrohyrax_meyeri),Thyrohyrax_pygmaeus)),(((((Bunohyrax_fajumensis,Pachyhyrax_crassidentatus),Bunohyrax_major),(Geniohyus_diphycus,Geniohyus_mirus)),(Megalohyrax_eocaenus,Megalohyrax_sp._nov.)),Thyrohyrax_litholagus)),(Dimaitherium_patnaiki,Namahyrax_corvus)),(Microhyrax_lavocati,Seggeurius_amourensis)));"
# cooper.tree<-read.tree(text=cooper.newick)
plot(cooper.tree)

#add branch length, computed.
cooper.tree<-compute.brlen(cooper.tree,method="Grafen")

#add tip dates to scale tree -------
dates.raw<-read_xlsx("../../phylogeny/hyracoid-tip-dates.xlsx")
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
ratios.raw<-data.jaws.all %>% group_by(Position,genus,species) %>% 
  summarise(l = round(mean(rel.length, na.rm=TRUE),3)) %>% 
  dcast(genus + species ~ Position)

#create a variable to match tree tip name formatting
ratios.raw$tip.name<-paste(ratios.raw$genus,ratios.raw$species,sep="_")

keep.these.tips<-which(cooper.fixed$tip.label %in% ratios.raw$tip.name)
drop.these.tips<-cooper.fixed$tip.label[-keep.these.tips]
cooper.trim<-drop.tip(cooper.fixed,drop.these.tips)
# cooper.trim<-drop.tip(cooper.trim,c("Seggeurius_amourensis"))
# cooper.trim<-drop.tip(cooper.trim,c("Thyrohyrax_litholagus","Thyrohyrax_pygmaeus",
#                                     "Thyrohyrax_meyeri",
#                                     "Thyrohyrax_domorictus","Prohyrax_hendeyi",
#                                    "Procavia_capensis"))
cooper.trim$node.label<-NULL #remove node labels entirely
plot(cooper.trim)

#Blomberg et al's K ------------
#pull the m2 ratio into a separate formatted object, following tutorials
m2.ratio<-ratios.raw$m2
names(m2.ratio)<-ratios.raw$tip.name
m2.ratio<-m2.ratio[complete.cases(m2.ratio)] #inferred names causing issues with caper downstream

#pull the m3 ratio into a separate formatted object, following tutorials
m3.ratio<-ratios.raw$m3

names(m3.ratio)<-ratios.raw$tip.name
m3.ratio<-m3.ratio[complete.cases(m3.ratio)]

print("m2, all taxa")
Kcalc(m2.ratio[cooper.fixed$tip.label], cooper.fixed)
phylosignal(m2.ratio[cooper.fixed$tip.label], cooper.fixed, reps = 1000)

print("m2 without Thyrohyrax/Procavaia")
Kcalc(m2.ratio[cooper.trim$tip.label], cooper.trim)
phylosignal(m2.ratio[cooper.trim$tip.label], cooper.trim, reps = 1000)

print("m3, all taxa")
Kcalc(m3.ratio[cooper.fixed$tip.label], cooper.fixed)
phylosignal(m3.ratio[cooper.fixed$tip.label], cooper.fixed, reps = 1000)

print("m3 without Thyrohyrax/Procavaia")
Kcalc(m3.ratio[cooper.trim$tip.label], cooper.trim)
phylosignal(m3.ratio[cooper.trim$tip.label], cooper.trim, reps = 1000)

# K with relative widths -------
ratios.width<-data.jaws %>% group_by(Position,genus,species) %>% 
  summarise(l = round(mean(rel.widths, na.rm=TRUE),3)) %>% 
  dcast(genus + species ~ Position)

#create a variable to match tree tip name formatting
ratios.width$tip.name<-paste(ratios.width$genus,ratios.width$species,sep="_")

keep.these.tips<-which(cooper.fixed$tip.label %in% ratios.width$tip.name)
drop.these.tips<-cooper.fixed$tip.label[-keep.these.tips]
cooper.trim<-drop.tip(cooper.fixed,drop.these.tips)
cooper.trim$node.label<-NULL #remove node labels entirely
plot(cooper.trim)

m1.w<-ratios.width$m1
names(m1.w)<-ratios.width$tip.name
print("m1, relative widths")
phylosignal(m1.w[cooper.trim$tip.label], cooper.trim, reps = 1000)

m2.w<-ratios.width$m2
names(m2.w)<-ratios.width$tip.name
print("m2, relative widths")
phylosignal(m2.w[cooper.trim$tip.label], cooper.trim, reps = 1000)

m3.w<-ratios.width$m3
names(m3.w)<-ratios.width$tip.name
print("m3, relative widths")
phylosignal(m3.w[cooper.trim$tip.label], cooper.trim, reps = 1000)

# #categorical D statistic, using library caper: underpowered -------
# ratios.bin.lit<-read.csv("ratios_differ_literature.csv")
# ratios.bin.lit$m2.binary<-0
# ratios.bin.lit$m2.binary[which(ratios.bin.lit$m1.shorter<0.05)]<-1
# 
# library(caper)
# # ratios.raw$m2.binary<-factor(ratios.raw$m2.binary)
# cooper.trim$node.label<-NULL #remove node labels entirely
# for.caper<-comparative.data(phy = cooper.trim, data = ratios.bin.lit, 
#                             names.col = tip.name,
#                             vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
# #NOTE: if you have node.labels in tree, won't work. Remove them. 
# phylo.d(for.caper,binvar=m2.binary,permut=10000)

# delta statistic, Borges, ape--------
ratios.bin.lit<-read.csv("ratios_differ_literature.csv")
ratios.bin.lit$m2.binary<-0
ratios.bin.lit$m2.binary[which(ratios.bin.lit$m1.shorter<0.05)]<-1

ratios.bin.lit$tip.name<-paste(ratios.bin.lit$genus,ratios.bin.lit$species,sep="_")

keep.these.tips<-which(cooper.fixed$tip.label %in% ratios.bin.lit$tip.name)
drop.these.tips<-cooper.fixed$tip.label[-keep.these.tips]
cooper.trim<-drop.tip(cooper.fixed,drop.these.tips)

source("../code.R")
tree<-cooper.trim
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
trait<-ratios.bin.lit$m2.binary
names(trait)<-ratios.bin.lit$tip.name
trait<-trait[cooper.trim$tip.label]
trait[1]<-1

#like D (and K also, turns out), delta can't handle an invariant trait
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100)


# simulation to confirm rtrait works -----
tree.rando<-rtree(100)
brown.rando<-rtrait(tree.rando, R = matrix(c(0, lambda.0, lambda.0, 0), 2),nstates=2)
plot(tree.rando,tip.color=c("red","blue")[brown.rando])
deltaA.rando <- delta(brown.rando,tree.rando,lambda.0,0.0589,10000,10,100)

for (i in 1:100){
  rtrait <- sample(brown.rando)
  random_delta[i] <- delta(rtrait,tree.rando,lambda.0,0.0589,10000,10,100)
} #this way of estimating p seems useless if you have an invariant or little-variant trait:
# shuffling won't change the value, or amt of information because it'll be same ancestral state either way
p_value <- sum(random_delta>deltaA)/length(random_delta) #p value is 1
#I'm not sure this is working...
#shouldn't this example simulated with rtrait definitely have phylogenetic signal?
boxplot(c(random_delta,deltaA.rando))
abline(h=deltaA.rando,col="red")
#oh, okay, so 1 is equivalent to a p values of <0.00001 or somethign small. So you could consider it 
#significant if p < 0.025 or p > 0.0975
#importantly, this test confirms that rtrait is doing what we expect: simulating a trait
#with signal that fits the tree structure. Now, apply to the smaller (less powerful) hyracoid tree:

# simulation to look at impact of increasing trait inertia on delta ------
nreps<-100
sim.delta<-list()
lambda.0<-0.1

for (i in 1:nreps){
  brown.traits.0<-rtrait(tree, R = matrix(c(0, lambda.0, lambda.0, 0), 2),nstates=2)
  # plot(tree,tip.color=c("red","blue")[brown.traits.0])
  deltaA.0 <- delta(brown.traits.0,tree,lambda.0,0.0589,10000,10,100)
  
  # simulate the character becoming more and more intertial (conserved) from root to tip
  brown.traits.sim<-brown.traits.0
  delta.vals<-deltaA.0
  for(tip in 1:(length(brown.traits.sim)-1)){
    brown.traits.sim[tip]<-2
    deltaA.sim <- delta(brown.traits.sim,tree,lambda.0,0.0589,10000,10,100)
    delta.vals<-c(delta.vals,deltaA.sim)
    #if the next step would create invariant character distribution, top and move to next
    if(length(which(brown.traits.sim==1))==1) {break}
  }
  sim.delta[[i]]<-delta.vals
}

ending.values<-sapply(sim.delta, tail, 1) %>% range
ending.values2<-sapply(sim.delta, tail, 2) %>% range

sim.table<-sapply(sim.delta, '[', seq(max(sapply(sim.delta, length)))) %>% as.data.frame
sim.table$tip<-seq(1:nrow(sim.table))
sim.long<-sim.table %>% melt(id.vars="tip",variable.name = "replicate")

write.csv(sim.table,"simulated_rel_length_delta.csv")


ggplot(sim.long, aes(x = tip, y = log(value), group = replicate)) + 
  geom_line(alpha = 0.2) +
  geom_hline(yintercept=log(ending.values),linetype="dashed",color="red")


# #plotting ------
# 
# backward.tree<-ggplot(bayesian.tree.mkv,aes(x,y)) + geom_tree() + geom_tiplab(fontface=3) + #could color branches by BPP with ",color=prob" in aes
#   theme_tree2() + #hexpand(0.02,direction=1) + #scale_x_continuous(labels = abs)
#   geom_range(range='height_0.95HPD',alpha=0.5,size=2) +
#   geom_nodelab(aes(x=branch,label=round(prob,2)),vjust=-0.5,size=2) #+
#   scale_color_continuous(low="gray",high="black") #+
#     #improve the color scheme
# time.scale<-seq(from=-70,to=0,by=10)
# revts(backward.tree) +scale_x_continuous(breaks=time.scale,labels=time.scale)
# ggsave("mkv-allcompat.pdf", width=122,units="mm") #FigX-Tree.pdf #will need some hand editing, so fix size. 