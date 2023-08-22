#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script


library(AlphaSimR)

rm(list = ls())


#first, we need a kinship cutoff value
founderPop = runMacs(nInd=4, nChr=15, segSites=5)
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

#methods for making clutches and taking samples
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

#taking a sample here is a little different - instead of subsetting the population
#by individual, we'll be taking a sample of the "reads"
#we'll assume that every haplotype from every individual will have an equal probability
#of being a sequenced read, so the function gets all the genotypes and then randomly samples
#the given number of reads from that pool
takeSample <- function(Clutch, depthPerInd){
  genos <- scoreMicrosats(Clutch)
  numReads <- round(Clutch@nInd * depthPerInd)
  
  if(numReads == 0){
    return(genos)
  }
  else{
    output <- data.frame(matrix(0, nrow = numReads, ncol = ncol(genos)))
    for(i in 1:numReads){
      for(j in 1:ncol(genos)){
        output[i,j] <- sample(genos[,j],1)
      }
    }
    return(output)
  }
}

#an alternative version that takes just a raw number of reads,
#instead of a sequencing depth
takeSample_numReads <- function(Clutch, numReads){
  genos <- scoreMicrosats(Clutch)

  if(numReads == 0){
    return(genos)
  }
  else{
    output <- data.frame(matrix(0, nrow = numReads, ncol = ncol(genos)))
    for(i in 1:numReads){
      for(j in 1:ncol(genos)){
        output[i,j] <- sample(genos[,j],1)
      }
    }
    return(output)
  }
}

#now this function reports the proportion of each allele in the pool of reads coming
#from the clutch
getAlleleProps <- function(Reads, allelesPerLocus){
  output <- data.frame(matrix(0, nrow = allelesPerLocus, ncol = ncol(Reads)))
  for(allele in 1:nrow(output)){
    for(locus in 1:ncol(output)){
      count <- 0
      for(read in 1:nrow(Reads)){
        if(Reads[read,locus] == allele){
          count <- count + 1
        }
      }
      rate <- count / nrow(Reads)
      if(is.na(rate)){
        output[allele,locus] <- 0
      }
      else{
        output[allele,locus] <- rate
      }
      
    }
  }
  return(output)
}

#now we want distributions for each allele at each locus after simulating it a bunch
#of times
getAlleleDists <- function(numSires, numLoci, allelesPerLocus, mainAlleleProp, percentSire1, seqDepth, reps){
  output <- array(0, dim = c(reps,allelesPerLocus,numLoci))
  for(rep in 1:dim(output)[1]){
    #so each rep, you make a new pair of parents and make a clutch of full-sibs
    founderPop = runMacs(nInd = numSires + 1, nChr = numLoci, segSites=5)
    SP = SimParam$new(founderPop)
    parents <- newPop(founderPop)
    parents <- setMicrosatHaplos(parents, allelesPerLocus, mainAlleleProp)
    clutchSize <- sample(c(300:1000),1)
    Clutch <- makeClutch(parents, clutchSize, percentSire1)
    #then you take reads from that clutch at some sequencing depth, and get the 
    #proportions of each allele
    fullGeno <- takeSample(Clutch,seqDepth)
    genoProps <- getAlleleProps(fullGeno, allelesPerLocus)
    #then copy that whole guy into the output array
    for(allele in 1:dim(output)[2]){
      for(locus in 1:dim(output)[3]){
        output[rep,allele,locus] <- genoProps[allele,locus]
      }
    }
  }
  return(output)
  }


#okay now for creating the dists

ALLELESPERLOCUS <- 4
REPS <- 1000

#dist1: 1 sire, equal props
dist1 <- getAlleleDists(numSires = 1, numLoci = 1, allelesPerLocus = ALLELESPERLOCUS, 
                        mainAlleleProp = 0, percentSire1 = 0, seqDepth = 2.5, reps = REPS)

#dist2: 2 sires, equal props
dist2 <- getAlleleDists(numSires = 2, numLoci = 1, allelesPerLocus = ALLELESPERLOCUS, 
                        mainAlleleProp = 0, percentSire1 = 0, seqDepth = 2.5, reps = REPS)

#dist3: 3 sires, equal props
dist3 <- getAlleleDists(numSires = 3, numLoci = 1, allelesPerLocus = ALLELESPERLOCUS, 
                        mainAlleleProp = 0, percentSire1 = 0, seqDepth = 2.5, reps = REPS)

#dist3_sireprop: 3 sires, unequal sire props
dist3_sireprops <- getAlleleDists(numSires = 3, numLoci = 1, allelesPerLocus = ALLELESPERLOCUS, 
                        mainAlleleProp = 0, percentSire1 = .8, seqDepth = 2.5, reps = REPS)

#dist3_alleleprop: 3 sires, unequal allele props
dist3_alleleprops <- getAlleleDists(numSires = 3, numLoci = 1, allelesPerLocus = ALLELESPERLOCUS, 
                                  mainAlleleProp = .8, percentSire1 = 0, seqDepth = 2.5, reps = REPS)


nf <- layout( matrix(c(1,2,3), ncol=1) )

hist(dist1[,1,1], breaks = 100, main = "Allele 1 Frequency Distribution\n1 sire",
     xlab = "Allele frequency (Allele 1)", xaxt='n',xlim = c(0,1))
axis(1,at=c(0,.25,.5,.75,1))

hist(dist2[,1,1], breaks = 100, main = "Allele 1 Frequency Distribution\n2 sires",
     xlab = "Allele frequency (Allele 1)", xaxt='n',xlim = c(0,1))
axis(1,at=c(0,.25,.5,.75,1))


hist(dist3[,1,1], breaks = 100, main = "Allele 1 Frequency Distribution with 3 sires\nEqual Sire, Equal Allele Proportions",
     xlab = "Allele frequency (Allele 1)", xaxt='n',xlim = c(0,1))
axis(1,at=c(0,.25,.5,.75,1))

hist(dist3_sireprops[,1,1], breaks = 100, main = "Allele 1 Frequency Distribution with 3 sires\nUnequal Sire, Equal Allele Proportions",
     xlab = "Allele frequency (Allele 1)", xaxt='n',xlim = c(0,1))
axis(1,at=c(0,.25,.5,.75,1))

hist(dist3_alleleprops[,1,1], breaks = 100, main = "Allele 1 Frequency Distribution with 3 sires\nEqual Sire, Unequal Allele Proportions",
     xlab = "Allele frequency (Allele 1)", xaxt='n',xlim = c(0,1))
axis(1,at=c(0,.25,.5,.75,1))



