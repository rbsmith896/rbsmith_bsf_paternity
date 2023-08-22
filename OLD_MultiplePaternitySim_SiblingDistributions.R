#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

#This script produces distributions of kinship coefficients
#for full-sibs and half-sibs with different numbers of loci
#and number of alleles per locus.
#Parameters can be set below, or at the end of the script
#under the header "making distributions"

NUMLOCI <- 15
ALLELESPERLOCUS <- 5
REPS <- 500




#Okay this time I'm gonna try to generalize the allele counts so it can get
#as big as you want it to be

#I think for simplicity I'll start by assuming that all the alleles at a 
#locus have equal frequencies, but this could maybe be fixed later

#creating the parents, with each chromosome corresponding to a different
#microsat locus (since they're in free recombination, as far as I can tell)
#Each locus is defined by 32 sites, for lots and lots of possible alleles
#way more than we need, but works nicely with the intToBits function
founderPop = runMacs(nInd=4, nChr=20, segSites=32)
SP = SimParam$new(founderPop)

#this method will take a population and a number of alleles desired at
#each locus, and set all the haplotypes to have that many alleles in
#equal frequencies
setMicrosatHaplos <- function(founderPop, richness){
  #adjusting the genetic map so there's no recombination
  GM <- SP$genMap
  for(i in 1:founderPop@nChr){
    for(j in 1:32){
      GM[[i]][j] <- GM[[i]][1]
    }
  }
  SP$switchGenMap(GM)
  
  allelesMatrix <- sapply(1:richness,function(x){ as.integer(intToBits(x))})
  Haplos <- pullSegSiteHaplo(founderPop)
  
  for(ind in 1:founderPop@nInd){
    for(loc in 1:founderPop@nChr){
      for(hap in 0:1){
        allele <- allelesMatrix[,sample(ncol(allelesMatrix),1)]
        for(i in 31:0){
          Haplos[(2*ind)-1+hap, (32*loc)-i] <- allele[32-i]
        }
      }
    }
  }
  
  setMarkerHaplo(founderPop, Haplos)
  return(founderPop)
}

#test it
popWithHaplos <- setMicrosatHaplos(founderPop, 10, SP)



#now let's make a function that takes this funny binary-encoded microsat system
#and converts it into a simple matrix of genotypes (with the microsat alleles)
#scored as 1 through the richness
scoreMicrosats <- function(Clutch){
  ClutchHaplo <- pullSegSiteHaplo(Clutch)
  output <- data.frame(matrix(nrow = nrow(ClutchHaplo), ncol = ncol(ClutchHaplo)/32))
  rownames(output) <- rownames(ClutchHaplo)
  colnames(output) <- c(1:Clutch@nChr)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      allele <- ClutchHaplo[i,(32*j-31):(32*j)]
      microsat <- 0
      fact <- 1
      for(k in 1:32){
        microsat <- microsat + allele[k]*fact
        fact <- 2*fact
      }
      output[i,j] <- microsat
    }
  }
  return(output)
}

#test
scoreMicrosats(popWithHaplos)

#now let's see if we're able to discern full sibs from half sibs.
#We'll make 2 pairs of clutches, one of full sibs and one of half sibs
founderPop = runMacs(nInd=3, nChr=30, segSites=32)
SP = SimParam$new(founderPop)
basePop = newPop(founderPop)
basePop <- setMicrosatHaplos(basePop, 10)
fullSib1 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 50)
fullSib2 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 50)
fullSibs <- mergePops(list(fullSib1,fullSib2))

halfSib1 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 50)
halfSib2 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 50)
halfSibs <- mergePops(list(halfSib1,halfSib2))

#now we need a function to get all the pairwise relationships between 2 populations
#but importantly, NOT among individuals in the same population

getPairwiseRelationships <- function(Clutch1, Clutch2){
  Haplos1 <- scoreMicrosats(Clutch1)
  Haplos2 <- scoreMicrosats(Clutch2)
  output <- data.frame(matrix(nrow = Clutch1@nInd, ncol = Clutch2@nInd))
  rownames(output) <- c(1:Clutch1@nInd)
  colnames(output) <- c(1:Clutch2@nInd)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      rel <- 0
      for(loc in 1:ncol(Haplos1)){
        #kind of confusing here, but basically just pulling the 4 microsatellite
        #genotypes, two for the first individual (i11 and i12) and 2 for the second
        #individual (i21 and i22)
        i11 <- Haplos1[2*i-1,loc]
        i12 <- Haplos1[2*i,loc]
        i21 <- Haplos2[2*j-1,loc]
        i22 <- Haplos2[2*j,loc]
        #then if there are 2 matches across individuals 
        #(which can happen 1 of 2 ways), score it as a 2
        if((i11==i21)&(i12==i22) | (i11==i22)&(i21==i12)){
          rel <- rel+2
        }
        #otherwise, if there are ANY matches across individuals,
        #score it was a 1
        else if(i11==i21 | i11==i22 | i12==i21 | i12==i22){
          rel <- rel+1
        }
        #the only other case is there are no matches, so scored a 0
      }
      #once you do that for all the loci, divide by the number of loci
      #and store the value in the corrsponding slot of the relationship matrix
      output[i,j] <- rel/Clutch1@nChr
    }
  }
  return(output)
}

fullsibpairs <- getPairwiseRelationships(fullSib1,fullSib2)

halfsibpairs <- getPairwiseRelationships(halfSib1,halfSib2)

hist(as.matrix(fullsibpairs))

hist(as.matrix(halfsibpairs))

hist(c(as.matrix(halfsibpairs),as.matrix(fullsibpairs)),breaks=25,
     main = "Kinship of 500 full-sibs and 500 half-sibs,\n20 loci with 10 alleles each",
     xlab = "2 * pairwise kinship coefficient")

#finally looks okay! With lots of offspring, lots of loci, lots of alleles


#And to use this, we'll make a method that creates these probability distributions
#It takes a value for the number of loci, and the number of alleles per locus, and then
#returns a list of 2 vectors, one the distribution of 500 full sib pairwise kinships (from parents
#that are newly generated every time), and one that's the distribution of 500 half sib pairwise 
#kinships (again from newly generated parents each time)
getKinshipDists <- function(numLoci, numAlleles, reps){
  fullSibDist <- c()
  parents <- runMacs(nInd=2, nChr=numLoci, segSites=32)
  SP <<- SimParam$new(parents)
  for(i in 1:reps){
    parentsNoHaplo <- newPop(parents)
    parentsWithHaplo <- setMicrosatHaplos(parentsNoHaplo, numAlleles)
    kid1 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kid2 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kinship <- getPairwiseRelationships(kid1,kid2)
    fullSibDist <- c(fullSibDist,kinship[1,1])
  }
  
  halfSibDist <- c()
  parents <- runMacs(nInd=3, nChr=numLoci, segSites=32)
  SP <<- SimParam$new(parents)
  for(i in 1:reps){
    parentsNoHaplo <- newPop(parents)
    parentsWithHaplo <- setMicrosatHaplos(parentsNoHaplo, numAlleles)
    kid1 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kid2 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 1)
    kinship <- getPairwiseRelationships(kid1,kid2)
    halfSibDist <- c(halfSibDist,kinship[1,1])
  }
  
  return(data.frame(FS = fullSibDist, HS = halfSibDist))
}


#MAKING DISTRIBUTIONS

dist <- getKinshipDists(NUMLOCI,ALLELESPERLOCUS, REPS)
par(mfrow = c(1,2))
hist(dist$FS, xlab = mean(dist$FS), breaks = 10)
hist(dist$HS, xlab = mean(dist$HS), breaks = 10)
#note: these are consistently centered at higher than we'd expect
#We want half sibs to be centered around .5, and full sibs around 1, but these
#are always greater - the only time they're coming down is when I turn the
#number of alleles per locus WAY WAY up
#I think this is because I'm calculating kinship using IBS, instead of IBD
#With like 100 alleles per locus, IBS always means IBD, but EVERY so often
#the parents will spawn with the same allele (say like they both have a 43)
#Then one kid could inherit a 43 from mom, and the other could inherit a 43
#from dad, and it would get scored as IBD in kinship when it's really just
#a coincidence

#With fewer alleles, this gets more and more likely, so average kinship gets 
#more and more inflated - I'm pretty sure?




#now let's think about getting a probability of half-sib vs full-sib given some coefficient of kinship
#returns a dataframe with 2 values, the probability of getting a given pairwise kinship value if the 
#2 are full sibs, and the probability of that value if the 2 are full sibs, given a full sib and half sib distribution
getSibProb <- function(kinship, Dist){
  return(data.frame(FullSibProb = sum(Dist[,1]==kinship) / length(Dist[,1]), 
                    HalfSibProb = sum(Dist[,2]==kinship) / length(Dist[,2])))
}

getSibProb(2, dist)
#the problem with this is that as the numbers get big, the chance of having some 
#specific discrete number gets lower and lower
#what we want is the probability that the value is in some kind of bin, right?
#can we assume normality, get the mean and standard deviation from the simulated
#data, and then just calculate a z-score?
#or should we organize the simulated data into bins, with each bin having its own
#probability, and then just sort observed data into the same bins and report
#the probability of each?









