#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

#This alternate script expands the original to incorporate 8 possible
#alleles at each microsatellite locus. I think this is too many, but good 
#to be able to do - allele count can be easily reduced by adjusting the 
#allele distribution

#The pairwise relationship values still seem to have too much overlap
#to confidently distribution full-sib from half-sib pairs though :()

#creating the parents, with each chromosome corresponding to a different
#microsat locus (since they're in free recombination, as far as I can tell)
#Each locus is defined by 3 sites, for 8 possible alleles
founderPop = runMacs(nInd=4, nChr=15, segSites=3)
SP = SimParam$new(founderPop)

#pulling the genetic map, and setting the location of both loci on each chromosome
#to the same value, so there's no recombination between them, 
#then storing it as the new genetic map
GM <- SP$genMap
for(i in 1:founderPop@nChr){
  GM[[i]][2] <- GM[[i]][1]
  GM[[i]][3] <- GM[[i]][1]
}
SP$switchGenMap(GM)

#now to set the haplotypes at each locus
#let's start by assuming all the loci have 8 alleles at equal frequencies
alleleDist <- list(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
Haplos <- pullSegSiteHaplo(founderPop)

#this makes sure that all the haplotypes are sampled from the distribution above
#this distribution can be different for different loci, but for now we're assuming
#it's the same for all
for(ind in 1:founderPop@nInd){
  for(loc in 1:founderPop@nChr){
    for(hap in 0:1){
      allele <- sample(alleleDist,1)
      Haplos[(2*ind)-1+hap,(3*loc)-2] <- allele[[1]][1]
      Haplos[(2*ind)-1+hap,(3*loc)-1] <- allele[[1]][2]
      Haplos[(2*ind)-1+hap,(3*loc)] <- allele[[1]][3]
    }
  }
}
setMarkerHaplo(founderPop,Haplos)


#now let's make a function that takes this funny binary-encoded microsat system
#and converts it into a simple matrix of genotypes (with the microsat alleles)
#scored as 1, 2, 3, 4, 5, 6, 7, 8
scoreMicrosats <- function(Clutch){
  ClutchHaplo <- pullSegSiteHaplo(Clutch)
  output <- data.frame(matrix(nrow = nrow(ClutchHaplo), ncol = ncol(ClutchHaplo)/3))
  rownames(output) <- rownames(ClutchHaplo)
  colnames(output) <- c(1:Clutch@nChr)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      x1 <- ClutchHaplo[i,3*j-2]
      x2 <- ClutchHaplo[i,3*j-1]
      x3 <- ClutchHaplo[i,3*j]
      if(x1==0){
        if(x2==0){
          if(x3==0){
            output[i,j] <- 1
          }
          else{
            output[i,j] <- 2
          }
        }
        else{
          if(x3==0){
            output[i,j] <- 3
          }
          else{
            output[i,j] <- 4
          }
        }
      }
      else{
        if(x2==0){
          if(x3==0){
            output[i,j] <- 5
          }
          else{
            output[i,j] <- 6
          }
        }
        else{
          if(x3==0){
            output[i,j] <- 7
          }
          else{
            output[i,j] <- 8
          }
        }
      }
    }
  }
  return(output)
}

#test
scoreMicrosats(founderPop)

#now for the function that makes the Clutches
#Clutch-making function, that takes:
#a) the number of offspring
#b) the number of sires (indicated by the size of the basePop - individual 1 is the mother, the rest are fathers)
#c) the relative contribution of the primary sire
#and produce a matrix of offspring genotypes
makeClutch <- function(basePop, clutchSize, percentSire1){
  numSires <- nInd(basePop) - 1
  sire1Brood <- round(percentSire1 * clutchSize, digits=0)
  altSireBrood <- round(((1 - percentSire1)/(numSires - 1))*clutchSize, digits=0)
  totalClutch <- sire1Brood + altSireBrood*(numSires-1)
  if(totalClutch != clutchSize){
    sire1Brood <- sire1Brood + (clutchSize - totalClutch)
  }
  Clutch <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = sire1Brood)
  for(i in 3:(numSires+1)){
    altcross <- makeCross(pop = basePop, crossPlan = matrix(c(1, i), ncol = 2), nProgeny = altSireBrood)
    Clutch <- mergePops(list(Clutch,altcross))
  }
  return(Clutch)
}

#now let's see if we're able to discern full sibs from half sibs.
#We'll make 2 pairs of clutches, one of full sibs and one of half sibs
basePop = newPop(founderPop)
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

hist(c(as.matrix(halfsibpairs),as.matrix(fullsibpairs)))

#unfortunately, this has the same problem! The distributions are centered
#on different values, as expected, but the combined distribution is
#not clearly bimodal. Conclusively distinguishing between full-sib and 
#half-sib pairs with a limited sample is still not feasible.


