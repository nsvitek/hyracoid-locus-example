#phylogenetic tree for Topernawi hyracoid tooth locus work

#intended for use after data are merged, 
#before the rest of the analyses are necessarily called in `hyracoid_analyses.R`

library(ape) #used throughouts
library(phytools) #I think maybe also used?
library(readxl) #for tip date spreadsheet
library(RRphylo) #for applying dates to tree, resolving polytomies

library(picante) #for Blomberg's k

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
cooper.trim<-drop.tip(cooper.trim,c("Thyrohyrax_litholagus","Thyrohyrax_pygmaeus",
                                    "Thyrohyrax_meyeri",
                                    "Thyrohyrax_domorictus","Prohyrax_hendeyi",
                                   "Procavia_capensis"))
cooper.trim$node.label<-NULL #remove node labels entirely
plot(cooper.trim)

#Blomberg et al's K with relative lengths (uppers and lowers) ------------
#pull the m2 ratio into a separate formatted object, following tutorials
m2.ratio<-ratios.raw$m2
names(m2.ratio)<-ratios.raw$tip.name
m2.ratio<-m2.ratio[complete.cases(m2.ratio)] #inferred names causing issues with caper downstream

#pull the m3 ratio into a separate formatted object, following tutorials
m3.ratio<-ratios.raw$m3

names(m3.ratio)<-ratios.raw$tip.name
m3.ratio<-m3.ratio[complete.cases(m3.ratio)]

sink("results_K.txt")
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

# K with relative widths (upper and lower) -------
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
m1.w<-m1.w[complete.cases(m1.w)]
print("m1, relative widths")
phylosignal(m1.w[cooper.trim$tip.label], cooper.trim, reps = 1000)


m2.w<-ratios.width$m2
names(m2.w)<-ratios.width$tip.name
m2.w<-m2.w[complete.cases(m2.w)]
print("m2, relative widths")
phylosignal(m2.w[cooper.trim$tip.label], cooper.trim, reps = 1000)

m3.w<-ratios.width$m3
names(m3.w)<-ratios.width$tip.name
m3.w<-m3.w[complete.cases(m3.w)]
print("m3, relative widths")
phylosignal(m3.w[cooper.trim$tip.label], cooper.trim, reps = 1000)
sink()

# delta statistic, Borges, ape--------
ratios.bin.lit<-read.csv("ratios_differ_literature.csv")
ratios.bin.lit$tip.name<-paste(ratios.bin.lit$genus,ratios.bin.lit$species,sep="_")
ratios.bin.lit$trait.binary<-0

trait.choice<-"m1.shorter"
source(paste(locateScripts,"retention_borges_delta.R",sep="/"))

trait.choice<-"m2.shorter"
source(paste(locateScripts,"retention_borges_delta.R",sep="/"))

