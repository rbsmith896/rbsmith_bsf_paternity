#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

#This first part creates a base population (with SNPs) and simulates
#some example clutches manually. After, there are some generalized functions
#to do this automatically with set parameters

#Once multi-allelic microsatellites are encoded, relationships can be 
#estimated using the pairwise relatedness coefficient r_xy, where r_xy = 2*theta_xy,
#and theta_xy is the probability that the two individuals have the same allele 
#at a given locus

#With known allele frequencies, we can simulate genotypes of full sibs, half sibs,
#and unrelated pairs to get distributions of r_xy for each pairing, and then 
#calculate the error rates that full-sibs are mis-classified as half-sibs and vice versa

#creating the parents, with each chromosome corresponding to a different
#microsat locus (since they're in free recombination, as far as I can tell)
#Each locus is defined by 2 sites, for 4 possible alleles
founderPop = runMacs(nInd=4, nChr=15, segSites=2)
SP = SimParam$new(founderPop)

#pulling the genetic map, and setting the location of both loci on each chromosome
#to the same value, so there's no recombination between them, 
#then storing it as the new genetic map
GM <- SP$genMap
for(i in 1:founderPop@nChr){
  GM[[i]][2] <- GM[[i]][1]
}
SP$switchGenMap(GM)

#now to set the haplotypes at each locus
#let's start by assuming all the loci have 3 alleles at equal frequencies
alleleDist <- list(c(1,0),c(0,1),c(1,1),c(0,0))
Haplos <- pullSegSiteHaplo(founderPop)

#this makes sure that all the haplotypes are sampled from the distribution above
#this distribution can be different for different loci, but for now we're assuming
#it's the same for all
for(ind in 1:founderPop@nInd){
  for(loc in 1:founderPop@nChr){
    for(hap in 0:1){
      allele <- sample(alleleDist,1)
      Haplos[(2*ind)-1+hap,(2*loc)-1] <- allele[[1]][1]
      Haplos[(2*ind)-1+hap,(2*loc)] <- allele[[1]][2]
    }
  }
}
setMarkerHaplo(founderPop,Haplos)


#now let's make a function that takes this funny binary-encoded microsat system
#and converts it into a simple matrix of genotypes (with the microsat alleles)
#scored as 1, 2, 3 and 4
scoreMicrosats <- function(Clutch){
  ClutchHaplo <- pullSegSiteHaplo(Clutch)
  output <- data.frame(matrix(nrow = nrow(ClutchHaplo), ncol = ncol(ClutchHaplo)/2))
  rownames(output) <- rownames(ClutchHaplo)
  colnames(output) <- c(1:Clutch@nChr)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      x1 <- ClutchHaplo[i,2*j-1]
      x2 <- ClutchHaplo[i,2*j]
      if(x1==1 & x2==0){
        output[i,j] <- 1
      }
      if(x1==0 & x2==1){
        output[i,j] <- 2
      }
      if(x1==1 & x2==1){
        output[i,j] <- 3
      }
      if(x1==0 & x2==0){
        output[i,j] <- 4
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

#this is a problem: the coefficients of relationship aren't different enough
#between full sibs and half sibs to be able to differentiate them

#need to figure out a better way to assess relationships





#OLDER THINGS

#I guess now we need a function to calculate relatedness between all individuals
#in a population?
getRelationships <- function(Clutch){
  Haplos <- scoreMicrosats(Clutch)
  output <- data.frame(matrix(nrow = Clutch@nInd, ncol = Clutch@nInd))
  rownames(output) <- c(1:Clutch@nInd)
  colnames(output) <- c(1:Clutch@nInd)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      rel <- 0
      for(loc in 1:ncol(Haplos)){
        #kind of confusing here, but basically just pulling the 4 microsatellite
        #genotypes, two for the first individual (i11 and i12) and 2 for the second
        #individual (i21 and i22)
        i11 <- Haplos[2*i-1,loc]
        i12 <- Haplos[2*i,loc]
        i21 <- Haplos[2*j-1,loc]
        i22 <- Haplos[2*j,loc]
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
      output[i,j] <- rel/Clutch@nChr
    }
  }
  return(output)
}

#and let's test
rel_fullSibs <- getRelationships(fullSibs)
rel_halfSibs <- getRelationships(halfSibs)



#and let's make a population and a sample clutch of 100
basePop = newPop(founderPop)
sampleClutch <- makeClutch(basePop, clutchSize = 100, percentSire1 <- .8)

clutchGeno <- scoreMicrosats(sampleClutch)

#now we can get measures of relatedness 



#Sampling function, which takes a certain number of random offspring 
#from a population object of offspring
#I couldn't figure out a nice way to select a sub-population without a phenotype,
#but the atttrition function can be finagled to do something similar.
#Not the most elegant, but it works for now
takeSample <- function(Clutch, sampNum){
  Samp <- attrition(Clutch, 1 - (sampNum / nInd(sampleClutch)))
  while(nInd(Samp) != sampNum){
    Samp <- attrition(Clutch, 1 - (sampNum / nInd(sampleClutch)))
  }
  return(Samp)
}

samp <- takeSample(sampleClutch, 15)










#Generalized functions




#example




#example
samp











#OLD STUFF, PROBABLY DELETE LATER


basePopHaplo = pullSegSiteHaplo(basePop)

basePopHaplo[, 1:5]
basePopHaplo

basePopGeno = pullSegSiteGeno(basePop)

basePopGeno[, 1:5]
basePopGeno

#single paternity model 1: female (ind1) mates with one male (ind2)
M1cross12 = makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1000)

#multiple paternity model 2, equal shares: female (ind1) mates with all 3 males equally
M2cross12 = makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 340)
M2cross13 = makeCross(pop = basePop, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 330)
M2cross14 = makeCross(pop = basePop, crossPlan = matrix(c(1, 4), ncol = 2), nProgeny = 330)

#multiple paternity model 3: female (ind1) mates with one male (ind2) with some contributions from the other males (ind3, ind4)
M3cross12 = makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 800)
M3cross13 = makeCross(pop = basePop, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 100)
M3cross14 = makeCross(pop = basePop, crossPlan = matrix(c(1, 4), ncol = 2), nProgeny = 100)

#progeny genotypes for model 1
M1cross12Geno = pullSegSiteGeno(M1cross12)
M1Geno <- M1cross12Geno
dim(M1Geno)

#progeny genotypes for model 2
M2cross12Geno = pullSegSiteGeno(M2cross12)
M2cross13Geno = pullSegSiteGeno(M2cross13)
M2cross14Geno = pullSegSiteGeno(M2cross14)
M2Geno <- rbind(M2cross12Geno,M2cross13Geno,M2cross14Geno)
dim(M2Geno)

#progeny genotypes for model 3
M3cross12Geno = pullSegSiteGeno(M3cross12)
M3cross13Geno = pullSegSiteGeno(M3cross13)
M3cross14Geno = pullSegSiteGeno(M3cross14)
M3Geno <- rbind(M3cross12Geno,M3cross13Geno,M3cross14Geno)
dim(M3Geno)

#getting correlations among model 1 offspring
M1GenoT <- t(M1Geno)
M1indCor = cor(M1GenoT)
corCols = hcl.colors(n = 21, palette = "RdYlBu",rev = TRUE)
image(M1indCor, xlab = "Individual", ylab = "Individual", axes = FALSE, col = corCols, main = "Model 1")

#then for model 2
M2GenoT <- t(M2Geno)
M2indCor = cor(M2GenoT)
corCols = hcl.colors(n = 21, palette = "RdYlBu",rev = TRUE)
image(M2indCor, xlab = "Individual", ylab = "Individual", axes = FALSE, col = corCols, main = "Model 2")

#and for model 3
M3GenoT <- t(M3Geno)
M3indCor = cor(M3GenoT)
corCols = hcl.colors(n = 21, palette = "RdYlBu",rev = TRUE)
image(M3indCor, xlab = "Individual", ylab = "Individual", axes = FALSE, col = corCols, main = "Model 3")

