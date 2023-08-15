#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggplot2")
library(ggplot2)

rm(list = ls())

#This script creates a clutch of a set size, with a given number of sires and a given
#proportional contribution from the primary sire (with given genetic parameters) and produces a matrix of kinship
#coefficients for the entire clutch.

#It also can take a sample of the clutch (below, under the header SAMPLE) and report
#the coefficients of kinship just within that sample

NUMLOCI <- 15
ALLELESPERLOCUS <- 5

CLUTCHSIZE <- 500
NUMSIRES <- 3
PROPMAINSIRE <- .8

SAMPLESIZE <- 20


#Genetic functions carrying over from sib-distribution script
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


#Setting up the parents
founderPop = runMacs(nInd = NUMSIRES + 1, nChr = NUMLOCI, segSites=32)
SP = SimParam$new(founderPop)
parents <- newPop(founderPop)
parents <- setMicrosatHaplos(parents, ALLELESPERLOCUS)


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


#creating the clutch
Clutch <- makeClutch(parents, CLUTCHSIZE, PROPMAINSIRE)

ClutchKinships <- getPairwiseRelationships(Clutch, Clutch)

#visualizing all the relationships in the clutch
KinshipPlot <- ClutchKinships %>% rownames_to_column() %>% gather(colname, value, -rowname)
ggplot(KinshipPlot, aes(x = rowname, y = colname, fill = value)) + geom_tile()

#so here we can clearly see the kinship blocks of the main sire, and the subsequent sires


#And then we can use this sampling function
#Sampling function, which takes a certain number of random offspring 
#from a population object of offspring
#I couldn't figure out a nice way to select a sub-population without a phenotype,
#but the attrition function can be finagled to do something similar.
#Not the most elegant, but it works for now
takeSample <- function(Clutch, sampNum){
  Samp <- attrition(Clutch, 1 - (sampNum / nInd(Clutch)))
  while(nInd(Samp) != sampNum){
    Samp <- attrition(Clutch, 1 - (sampNum / nInd(Clutch)))
  }
  return(Samp)
}

Sample <- takeSample(Clutch, SAMPLESIZE)
SamplehKinships <- getPairwiseRelationships(Sample, Sample)

#visualizing the relationships in the sample
SampleKinshipPlot <- SamplehKinships %>% rownames_to_column() %>% gather(colname, value, -rowname)
ggplot(SampleKinshipPlot, aes(x = rowname, y = colname, fill = value)) + geom_tile()

#relationship blocks can be seen here, but not completely obvious. How to decide this statistically?
Sample@father


