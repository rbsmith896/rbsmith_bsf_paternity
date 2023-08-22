#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

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
popWithHaplos <- setMicrosatHaplos(founderPop, 16)



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
fullSib1 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 500)
fullSib2 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 500)
fullSibs <- mergePops(list(fullSib1,fullSib2))

halfSib1 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 500)
halfSib2 <- makeCross(pop = basePop, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 500)
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

hist(c(as.matrix(halfsibpairs),as.matrix(fullsibpairs)),breaks=57)

#finally looks okay! With lots of offspring, lots of loci, lots of alleles

