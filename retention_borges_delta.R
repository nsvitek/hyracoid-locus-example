column.index<-which(colnames(ratios.bin.lit)==trait.choice)
ratios.bin.lit$trait.binary[which(ratios.bin.lit[,column.index]<0.05)]<-1
keep.these.tips<-which(cooper.fixed$tip.label %in% ratios.bin.lit$tip.name)
drop.these.tips<-cooper.fixed$tip.label[-keep.these.tips]
cooper.trim<-drop.tip(cooper.fixed,drop.these.tips)

source("../code.R")
tree<-cooper.trim
tree$edge.length[tree$edge.length==0] <- quantile(tree$edge.length,0.1)*0.1
trait<-ratios.bin.lit$trait.binary
names(trait)<-ratios.bin.lit$tip.name
trait<-trait[cooper.trim$tip.label]

#like D (and K also, turns out), delta can't handle an invariant trait
deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100)

# # simulation to confirm rtrait works -----
# tree.rando<-rtree(100) #simulate a large tree, larger than small empirical tree
# lambda.0<-0.1
# brown.rando<-rtrait(tree.rando, R = matrix(c(0, lambda.0, lambda.0, 0), 2),nstates=2)
# plot(tree.rando,tip.color=c("red","blue")[brown.rando])
# deltaA.rando <- delta(brown.rando,tree.rando,lambda.0,0.0589,10000,10,100)
# random_delta<-list()
# for (i in 1:100){
#   rtrait <- sample(brown.rando)
#   random_delta[i] <- delta(rtrait,tree.rando,lambda.0,0.0589,10000,10,100)
# } #this way of estimating p seems useless if you have an invariant or little-variant trait:
# # shuffling won't change the value, or amt of information because it'll be same ancestral state either way
# random_delta<-unlist(random_delta)
# p_value <- sum(random_delta>deltaA.rando)/length(random_delta) #p value is 1
# #I'm not sure this is working...
# #shouldn't this example simulated with rtrait definitely have phylogenetic signal?
# boxplot(c(random_delta,deltaA.rando))
# abline(h=deltaA.rando,col="red")
# #oh, okay, so 1 is equivalent to a p values of <0.00001 or sommething small. So you could consider it 
# #significant if p < 0.025 or p > 0.0975
# #importantly, this test confirms that rtrait is doing what we expect: simulating a trait
# #with signal that fits the tree structure. Now, apply to the smaller (less powerful) hyracoid tree:

# simulation to look at impact of increasing trait inertia on delta ------
nreps<-1000
sim.delta<-matrix(NA,nrow = (length(tree$tip.label)-1), ncol = nreps )
lambda.0<-0.1
initial.deltas<-matrix(NA,ncol=2, nrow=nreps)
for (i in 1:nreps){
  brown.traits.0<-rtrait(tree, R = matrix(c(0, lambda.0, lambda.0, 0), 2),nstates=2)
  # plot(tree,tip.color=c("red","blue")[brown.traits.0])
  deltaA.0 <- delta(brown.traits.0,tree,lambda.0,0.0589,10000,10,100)
  #write to one vector for p-values
  initial.deltas[i,1]<-deltaA.0
  trait.counter<-length(which(brown.traits.0==2))
  initial.deltas[i,2]<-trait.counter
  #also write to the simulation table for plot
  sim.delta[trait.counter,i]<-deltaA.0
  # simulate the character becoming more and more retained (conserved) from root to tip
  brown.traits.sim<-brown.traits.0
  for(tip in 1:(length(brown.traits.sim)-1)){
    brown.traits.sim[tip]<-2
    trait.counter<-length(which(brown.traits.sim==2))
    if (is.na(sim.delta[trait.counter,i])) 
    {# print(trait.counter)
      deltaA.sim <- delta(brown.traits.sim,tree,lambda.0,0.0589,10000,10,100)
      # print(deltaA.sim)
      #write to sim.delta
      sim.delta[trait.counter,i]<-deltaA.sim
      #if the next step would create invariant character distribution, top and move to next
      if(length(which(brown.traits.sim==1))==1) {break}
      }
    else {next}
  }
}

write.csv(sim.delta,paste("simulated_delta_",trait.choice,".csv",sep=""))
sim.table<-as.data.frame(sim.delta)
sim.table$number.derived.tips<-seq(1:nrow(sim.delta))

sim.long<-sim.table %>% melt(id.vars="number.derived.tips",variable.name = "replicate")
sim.long$starting.val <- FALSE
for(i in 1:nrow(initial.deltas)){
  sim.long$starting.val[which(sim.long$number.derived.tips==initial.deltas[i,2]&sim.long$value==initial.deltas[i,1])]<-TRUE
}

#plot results
ggplot(sim.long, aes(x = number.derived.tips, y = log(1+value), color = starting.val, group=replicate)) + 
  geom_line(alpha = 0.2) +
  geom_jitter(alpha= 0.4,width=0.25) + 
  geom_point(x = length(which(trait==1)), y = log(1+deltaA), color = "black")+
  scale_color_manual(values=c("gray","red")) + theme_minimal() +
  labs(x = "Number of derived tips", y = "log(1+ \u03B4)",parse=TRUE) 

ggsave(paste("delta_",trait.choice,".pdf",sep=""),device = cairo_pdf, width = single.column.width,
       height = single.column.width,units="in",dpi=600)
# p value for comparison? -----
#In this case, the null is a tree with Brownian motion-type phylogenetic signal
p.delta<-length(which(log(1+initial.deltas[,1])>=(log(1+deltaA))))/nreps

initial.deltas[which(log(1+initial.deltas[,1])>=(log(1+deltaA))),]

#and the alternative is phylogenetic retention, which should have a significantly higher delta values
p.tips<-length(which(initial.deltas[,2]>=length(which(trait==1))))/nreps

p.csv<-cbind(p.delta,p.tips)
write.csv(p.csv,paste("p_delta_",trait.choice,".csv",sep=""))
