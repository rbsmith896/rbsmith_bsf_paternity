#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

library(AlphaSimR)

rm(list = ls())

#initial parameters for setting up and testing functions
NUMLOCI <- 15
ALLELESPERLOCUS <- 3 #this can be between 1 and 31 - since the script encodes alleles as 5 bit integers, it can't handle allele counts above 31
NUMSIRES <- 3
PROPMAINSIRE <- .8 #if this is zero, it'll do equal proportions
MAINALLELEPROP <- .8 #if this is zero, it'll do equal proportions
SAMPLESIZE <- 10 #if this is zero, it'll do the whole clutch
REPS <- 5
KINSHIPPROPCUTOFF <- .9 #this is the confidence with which we want to 
#be able to call half sibs from full sibs - the more stringent this 
#is (closer to 1), the lower our Type I error (less likely to accidentally)
#call full-sibs as half-sibs, but also the lower our power (more likely
#to miss true half-sibs)
#if it's 1, it'll make the cutoff the lowest full-sib relatedness

founderPop = runMacs(nInd=NUMSIRES+1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)

setMicrosatHaplos <- function(founderPop, richness, mainalleleprop){
  #adjusting the genetic map so there's no recombination
  GM <- SP$genMap
  for(i in 1:founderPop@nChr){
    for(j in 1:5){
      GM[[i]][j] <- GM[[i]][1]
    }
  }
  SP$switchGenMap(GM)
  
  allelesMatrix <- sapply(1:richness,function(x){ as.integer(intToBits(x)[1:5])})
  Haplos <- pullSegSiteHaplo(founderPop)
  
  if(mainalleleprop == 0){
    allelesDist <- c(1:richness)
  }
  else{
    allelesDist <- rep(1,100*mainalleleprop)
    remain <- round((1-mainalleleprop)*100 / (richness-1))
    for(i in 2:richness){
      allelesDist <- c(allelesDist, rep(i,remain))
    }
  }
  
  for(ind in 1:founderPop@nInd){
    for(loc in 1:founderPop@nChr){
      for(hap in 0:1){
        allele <- allelesMatrix[,sample(allelesDist,1)]
        for(i in 4:0){
          Haplos[(2*ind)-1+hap, (5*loc)-i] <- allele[5-i]
        }
      }
    }
  }
  
  setMarkerHaplo(founderPop, Haplos)
  return(founderPop)
}

scoreMicrosats <- function(Clutch){
  ClutchHaplo <- pullSegSiteHaplo(Clutch)
  output <- data.frame(matrix(nrow = nrow(ClutchHaplo), ncol = ncol(ClutchHaplo)/5))
  rownames(output) <- rownames(ClutchHaplo)
  colnames(output) <- c(1:Clutch@nChr)
  for(i in 1:nrow(output)){
    for(j in 1:ncol(output)){
      allele <- ClutchHaplo[i,(5*j-4):(5*j)]
      microsat <- 0
      fact <- 1
      for(k in 1:5){
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

getKinshipDists <- function(numLoci, numAlleles,mainAlleleProp, reps){
  fullSibDist <- c()
  parents <- runMacs(nInd=2, nChr=numLoci, segSites=5)
  SP <<- SimParam$new(parents)
  for(i in 1:reps){
    parentsNoHaplo <- newPop(parents)
    parentsWithHaplo <- setMicrosatHaplos(parentsNoHaplo, numAlleles, mainAlleleProp)
    kid1 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kid2 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kinship <- getPairwiseRelationships(kid1,kid2)
    fullSibDist <- c(fullSibDist,kinship[1,1])
  }
  
  halfSibDist <- c()
  parents <- runMacs(nInd=3, nChr=numLoci, segSites=5)
  SP <<- SimParam$new(parents)
  for(i in 1:reps){
    parentsNoHaplo <- newPop(parents)
    parentsWithHaplo <- setMicrosatHaplos(parentsNoHaplo, numAlleles, mainAlleleProp)
    kid1 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kid2 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 3), ncol = 2), nProgeny = 1)
    kinship <- getPairwiseRelationships(kid1,kid2)
    halfSibDist <- c(halfSibDist,kinship[1,1])
  }
  
  return(data.frame(FS = fullSibDist, HS = halfSibDist))
}

getFullSibDists <- function(numLoci, numAlleles,mainAlleleProp, reps){
  fullSibDist <- c()
  parents <- runMacs(nInd=2, nChr=numLoci, segSites=5)
  SP <<- SimParam$new(parents)
  for(i in 1:reps){
    parentsNoHaplo <- newPop(parents)
    parentsWithHaplo <- setMicrosatHaplos(parentsNoHaplo, numAlleles, mainAlleleProp)
    kid1 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kid2 <- makeCross(pop = parentsWithHaplo, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = 1)
    kinship <- getPairwiseRelationships(kid1,kid2)
    fullSibDist <- c(fullSibDist,kinship[1,1])
  }
  
  return(data.frame(FS = fullSibDist))
}


#now for making all the clutches, recording in how many MP is detected
makeClutch <- function(basePop, clutchSize, percentSire1){
  numSires <- nInd(basePop) - 1
  if(numSires == 1){
    Clutch <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), 
                        nProgeny = clutchSize)
  }
  else if(percentSire1 == 0){
    Clutch <- makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), 
                        nProgeny = round(clutchSize / numSires))
    for(i in 3:(numSires+1)){
      altcross <- makeCross(pop = basePop, crossPlan = matrix(c(1, i), ncol = 2), 
                            nProgeny = round(clutchSize / numSires))
      Clutch <- mergePops(list(Clutch,altcross))
    }
  }
  else{
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
  }
  return(Clutch)
}

takeSample <- function(Clutch, sampNum){
  if(sampNum == 0){
    return(Clutch)
  }
  else{
    Samp <- attrition(Clutch, 1 - (sampNum / nInd(Clutch)))
    while(nInd(Samp) != sampNum){
     Samp <- attrition(Clutch, 1 - (sampNum / nInd(Clutch)))
   }
   return(Samp)
  }
}

getCutoffKinship <- function(numLoci, allelesPerLocus, mainAlleleProp, kinshipPropCutoff, reps){
  dist <- getFullSibDists(numLoci,allelesPerLocus,mainAlleleProp, reps)
  FullSibs <- sort(dist$FS)
  print(head(FullSibs))
  if(kinshipPropCutoff == 1){
    CUTOFF <- FullSibs[1]
  }else{
    CUTOFF <- FullSibs[reps*(1-kinshipPropCutoff)]
  }
  return(CUTOFF)
}

getMPprop <- function(numSires, numLoci, allelesPerLocus, mainAlleleProp,
                      propMainSire, sampleSize, kinshipCutoff, reps){
  scoreCard <- rep(0,reps)
  for(rep in 1:reps){
    founderPop = runMacs(nInd = numSires + 1, nChr = numLoci, segSites=5)
    SP = SimParam$new(founderPop)
    parents <- newPop(founderPop)
    parents <- setMicrosatHaplos(parents, allelesPerLocus, mainAlleleProp)
    #assuming clutches are between 300 and 1000 eggs
    clutchSize <- sample(c(300:1000),1)
    Clutch <- makeClutch(parents, clutchSize, propMainSire)
    
    Sample <- takeSample(Clutch, sampleSize)
    
    SampleKinships <- getPairwiseRelationships(Sample, Sample)
    MPcall <- 0
    for(i in 1:sampleSize){
      MP <- rep(0,sampleSize)
      for(j in 1:sampleSize){
        if(SampleKinships[i,j] < kinshipCutoff){
          MP[j] <- 1
        }
      }
      if(sum(MP)/length(MP) > .5){
        MPcall <- 1
      }
    }
    scoreCard[rep] <- MPcall
  }
  
  MPprop <- sum(scoreCard)/length(scoreCard)
  return(MPprop)
}

#alternate version that just looks at the actual fathers represented in the sample
#and reports the "true" proportion of multiple paternity
getTrueMPProp <- function(numSires, numLoci, allelesPerLocus, mainAlleleProp,
                          propMainSire, sampleSize, reps){
  scoreCard <- rep(0,reps)
  for(rep in 1:reps){
    founderPop = runMacs(nInd = numSires + 1, nChr = numLoci, segSites=5)
    SP = SimParam$new(founderPop)
    parents <- newPop(founderPop)
    parents <- setMicrosatHaplos(parents, allelesPerLocus, mainAlleleProp)
    #assuming clutches are between 300 and 1000 eggs
    clutchSize <- sample(c(300:1000),1)
    Clutch <- makeClutch(parents, clutchSize, propMainSire)
    
    Sample <- takeSample(Clutch, sampleSize)
    
    dads <- as.integer(Sample@father)
    MPcall <- 0
    for(d in 1:sampleSize){
      if(dads[d] != dads[1]){
        MPcall <- 1
      }
    }
    
    scoreCard[rep] <- MPcall
  }
  
  MPprop <- sum(scoreCard)/length(scoreCard)
  return(MPprop)
}

#---------------------------------------------------------------------------------------------

NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

#plot1a: showing the proportion of true MP samples with equal sire proportions
plot1a <- data.frame(row.names = c(2, 3, 4, 5, 10, 25, 50, 100))
ssArray <- c(2, 3, 4, 5, 10, 25, 50, 100)
plot1a$sires2 <- rep(0,8)
for(i in 1:8){
  plot1a$sires2[i] <- getTrueMPProp(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], reps = REPS)
}
plot1a$sires3 <- rep(0,8)
for(i in 1:8){
  plot1a$sires3[i] <- getTrueMPProp(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], reps = REPS)
}
plot1a$sires4 <- rep(0,8)
for(i in 1:8){
  plot1a$sires4[i] <- getTrueMPProp(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], reps = REPS)
}


#---------------------------------------------------------------------------------------------


#plot1b: showing the proportion of true MP samples with UNequal sire proportions
plot1b <- data.frame(row.names = c(2, 3, 4, 5, 10, 25, 50, 100))
ssArray <- c(2, 3, 4, 5, 10, 25, 50, 100)
plot1b$sires2 <- rep(0,8)
for(i in 1:8){
  plot1b$sires2[i] <- getTrueMPProp(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], reps = REPS)
}
plot1b$sires3 <- rep(0,8)
for(i in 1:8){
  plot1b$sires3[i] <- getTrueMPProp(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], reps = REPS)
}
plot1b$sires4 <- rep(0,8)
for(i in 1:8){
  plot1b$sires4[i] <- getTrueMPProp(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS,
                                    mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], reps = REPS)
}


#---------------------------------------------------------------------------------------------



NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

#plot2a: showing the proportion of clutches called for MP using the method
plot2a <- data.frame(row.names = c(10, 25, 50, 100))
ssArray <- c(10, 25, 50, 100)
CUTOFF <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                           mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 10000)
plot2a$sire1 <- rep(0,4)
for(i in 1:4){
  plot2a$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                               kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a$sires2 <- rep(0,4)
for(i in 1:4){
  plot2a$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a$sires3 <- rep(0,4)
for(i in 1:4){
  plot2a$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a$sires4 <- rep(0,4)
for(i in 1:4){
  plot2a$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}


#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

#plot2a: showing the proportion of clutches called for MP using the method
plot2a_ls <- data.frame(row.names = c(2, 3, 4, 5))
ssArray <- c(2, 3, 4, 5)
plot2a_ls$sire1 <- rep(0,4)
for(i in 1:4){
  plot2a_ls$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                               kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a_ls$sires2 <- rep(0,4)
for(i in 1:4){
  plot2a_ls$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a_ls$sires3 <- rep(0,4)
for(i in 1:4){
  plot2a_ls$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}
plot2a_ls$sires4 <- rep(0,4)
for(i in 1:4){
  plot2a_ls$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = 0, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF, reps = REPS)
}


#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(10, 25, 50, 100)
CUTOFF_equalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                     mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 1000)
#CUTOFF_unequalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
#mainAlleleProp = .8, kinshipPropCutoff = .9, reps = 10000)

#plot3b: showing proportion called for MP with: UNequal sire pop, equal allele freqs
plot3b <- data.frame(row.names = c(10, 25, 50, 100))

plot3b$sire1 <- rep(0,4)
for(i in 1:4){
  plot3b$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                               kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires2 <- rep(0,4)
for(i in 1:4){
  plot3b$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires3 <- rep(0,4)
for(i in 1:4){
  plot3b$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires4 <- rep(0,4)
for(i in 1:4){
  plot3b$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}

plot3b_bigs <- plot3b

#---------------------------------------------------------------------------------------------



NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(2, 3, 4, 5)
CUTOFF_equalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                     mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 1000)
#CUTOFF_unequalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
#mainAlleleProp = .8, kinshipPropCutoff = .9, reps = 10000)

#plot3b: showing proportion called for MP with: UNequal sire pop, equal allele freqs
plot3b <- data.frame(row.names = c(2, 3, 4, 5))

plot3b$sire1 <- rep(0,4)
for(i in 1:4){
  plot3b$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                               kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires2 <- rep(0,4)
for(i in 1:4){
  plot3b$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires3 <- rep(0,4)
for(i in 1:4){
  plot3b$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}
plot3b$sires4 <- rep(0,4)
for(i in 1:4){
  plot3b$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = 0, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_equalFreq, reps = REPS)
}


plot3b_smalls <- plot3b

#---------------------------------------------------------------------------------------------



NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(10, 25, 50, 100)
#CUTOFF_equalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
#mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 10000)
CUTOFF_unequalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                       mainAlleleProp = .8, kinshipPropCutoff = .9, reps = 1000)

#plot3c: showing proportion called for MP with: equal sire pop, UNequal allele freqs
plot3c_bigs <- data.frame(row.names = c(10, 25, 50, 100))

plot3c_bigs$sire1 <- rep(0,4)
for(i in 1:4){
  plot3c_bigs$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                    mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                    kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_bigs$sires2 <- rep(0,4)
for(i in 1:4){
  plot3c_bigs$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                     mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                     kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_bigs$sires3 <- rep(0,4)
for(i in 1:4){
  plot3c_bigs$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                     mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                     kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_bigs$sires4 <- rep(0,4)
for(i in 1:4){
  plot3c_bigs$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                     mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                     kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}

#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(2, 3, 4, 5)
#CUTOFF_equalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
#mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 10000)
CUTOFF_unequalFreq <- 1.46666666666667

#plot3c: showing proportion called for MP with: equal sire pop, UNequal allele freqs
plot3c_smalls <- data.frame(row.names = c(2, 3, 4, 5))

plot3c_smalls$sire1 <- rep(0,4)
for(i in 1:4){
  plot3c_smalls$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                      mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                      kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_smalls$sires2 <- rep(0,4)
for(i in 1:4){
  plot3c_smalls$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                       mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                       kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_smalls$sires3 <- rep(0,4)
for(i in 1:4){
  plot3c_smalls$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                       mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                       kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3c_smalls$sires4 <- rep(0,4)
for(i in 1:4){
  plot3c_smalls$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                       mainAlleleProp = .8, propMainSire = 0, sampleSize = ssArray[i], 
                                       kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}


#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(2, 3, 4, 5, 10, 25, 50, 100)
#CUTOFF_equalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
#mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 10000)
CUTOFF_unequalFreq <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                       mainAlleleProp = .8, kinshipPropCutoff = .9, reps = 10000)

#plot3d: showing proportion called for MP with: UNequal sire pop, UNequal allele freqs
plot3d <- data.frame(row.names = c(2, 3, 4, 5, 10, 25, 50, 100))

plot3d$sire1 <- rep(0,8)
for(i in 1:8){
  plot3d$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = .8, propMainSire = .8, sampleSize = ssArray[i], 
                               kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3d$sires2 <- rep(0,8)
for(i in 1:8){
  plot3d$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = .8, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3d$sires3 <- rep(0,8)
for(i in 1:8){
  plot3d$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = .8, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}
plot3d$sires4 <- rep(0,8)
for(i in 1:8){
  plot3d$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                mainAlleleProp = .8, propMainSire = .8, sampleSize = ssArray[i], 
                                kinshipCutoff = CUTOFF_unequalFreq, reps = REPS)
}

#---------------------------------------------------------------------------------------------



NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(.8, .85, .9, .95, .99)
kcArray <- c(0, 0, 0, 0, 0)
for(i in 1:5){
  kcArray[i] <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                 mainAlleleProp = 0, kinshipPropCutoff = ssArray[i], reps = 1000)
}

#plot4: showing the effect of changing the kinship cutoff
plot4 <- data.frame(row.names = c(.8, .85, .9, .95, .99))

plot4$sire1 <- rep(0,5)

for(i in 1:5){
  plot4$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                              mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                              kinshipCutoff = kcArray[i], reps = REPS)
}
plot4$sires2 <- rep(0,5)
for(i in 1:5){
  plot4$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot4$sires3 <- rep(0,5)
for(i in 1:5){
  plot4$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot4$sires4 <- rep(0,5)
for(i in 1:5){
  plot4$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}

#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 1000

ssArray <- c(.55, .6, .65, .7, .75)
kcArray <- c(0, 0, 0, 0, 0)
for(i in 1:5){
  kcArray[i] <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                 mainAlleleProp = 0, kinshipPropCutoff = ssArray[i], reps = 1000)
}

#plot4: showing the effect of changing the kinship cutoff
plot4_ls <- data.frame(row.names = c(.55, .6, .65, .7, .75))

plot4_ls$sire1 <- rep(0,5)

for(i in 1:5){
  plot4_ls$sire1[i] <- getMPprop(numSires = 1, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                              mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                              kinshipCutoff = kcArray[i], reps = REPS)
}
plot4_ls$sires2 <- rep(0,5)
for(i in 1:5){
  plot4_ls$sires2[i] <- getMPprop(numSires = 2, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot4_ls$sires3 <- rep(0,5)
for(i in 1:5){
  plot4_ls$sires3[i] <- getMPprop(numSires = 3, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot4_ls$sires4 <- rep(0,5)
for(i in 1:5){
  plot4_ls$sires4[i] <- getMPprop(numSires = 4, numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}

#---------------------------------------------------------------------------------------------


ALLELESPERLOCUS <- 4
REPS <- 1000

ssArray <- c(5, 10, 15, 20, 25)
kcArray <- c(0, 0, 0, 0, 0)
for(i in 1:5){
  founderPop = runMacs(nInd=1, nChr= ssArray[i], segSites=5)
  SP = SimParam$new(founderPop)
  kcArray[i] <- getCutoffKinship(numLoci = ssArray[i], allelesPerLocus = ALLELESPERLOCUS, 
                                 mainAlleleProp = 0, kinshipPropCutoff = .9, reps = 1000)
}

#plot6: showing the effect of changing the number of loci
plot6 <- data.frame(row.names = c(5, 10, 15, 20, 25))

plot6$sire1 <- rep(0,5)

for(i in 1:5){
  founderPop = runMacs(nInd=2, nChr= ssArray[i], segSites=5)
  SP = SimParam$new(founderPop)
  plot6$sire1[i] <- getMPprop(numSires = 1, numLoci = ssArray[i], allelesPerLocus = ALLELESPERLOCUS, 
                              mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                              kinshipCutoff = kcArray[i], reps = REPS)
}
plot6$sires2 <- rep(0,5)
for(i in 1:5){
  founderPop = runMacs(nInd=3, nChr= ssArray[i], segSites=5)
  SP = SimParam$new(founderPop)
  plot6$sires2[i] <- getMPprop(numSires = 2, numLoci = ssArray[i], allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot6$sires3 <- rep(0,5)
for(i in 1:5){
  founderPop = runMacs(nInd=4, nChr= ssArray[i], segSites=5)
  SP = SimParam$new(founderPop)
  plot6$sires3[i] <- getMPprop(numSires = 3, numLoci = ssArray[i], allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}
plot6$sires4 <- rep(0,5)
for(i in 1:5){
  founderPop = runMacs(nInd=5, nChr= ssArray[i], segSites=5)
  SP = SimParam$new(founderPop)
  plot6$sires4[i] <- getMPprop(numSires = 4, numLoci = ssArray[i], allelesPerLocus = ALLELESPERLOCUS, 
                               mainAlleleProp = 0, propMainSire = 0, sampleSize = 25, 
                               kinshipCutoff = kcArray[i], reps = REPS)
}


#---------------------------------------------------------------------------------------------


NUMLOCI <- 15
ALLELESPERLOCUS <- 4
founderPop = runMacs(nInd=1, nChr=NUMLOCI, segSites=5)
SP = SimParam$new(founderPop)
REPS <- 10000

ssArray <- c(.8, .85, .9, .95, .99)

#plot4: showing the effect of changing the kinship cutoff
FinalKCs <- data.frame(row.names = c(.8, .85, .9, .95, .99))

FinalKCs$equalAlleleProps <- rep(0,5)

FinalKCs$UNequalAlleleProps <- rep(0,5)

for(i in 1:5){
  FinalKCs$equalAlleleProps[i] <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                                   mainAlleleProp = 0, kinshipPropCutoff = ssArray[i], reps = REPS)
}

for(i in 1:5){
  FinalKCs$UNequalAlleleProps[i] <- getCutoffKinship(numLoci = NUMLOCI, allelesPerLocus = ALLELESPERLOCUS, 
                                                     mainAlleleProp = .8, kinshipPropCutoff = ssArray[i], reps = REPS)
}








