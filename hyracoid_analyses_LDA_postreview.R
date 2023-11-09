#Hyracoid locus code after revision

#Request for LDA
#Technically, most LDA intended for multivariate, most tests in paper are univariate, so work up

#Start wtih univariate case, well-sampled empirical:
#Procavia and/or T. domorictus
#Use mean and variance of ratios for 2 loci, following p. 208 of Bookstein 2018 textbook
#Calculate, and visualize with similar diagrams, the distribution and discrimination.
#Then for domorictus, use modeled SD and see if it changes univariate answer
#Then run empirical LDA for capensis
#Then use modeled distribution for domorictus LDA, see how it goes
#If OK, then use modelled LDA for everyone for diagnosibility
#Then apply to Meroehyrax: model as if it did not have reference mandible. Plot together, loook for breakpoints at expected ratio difference

pchisq(10.96254,1,lower.tail=FALSE)
pchisq(4,1,lower.tail=FALSE)
