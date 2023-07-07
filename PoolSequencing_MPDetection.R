#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script


install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())


NUMLOCI <- 1
ALLELESPERLOCUS <- 3
NUMSIRES <- 3
PROPMAINSIRE <- 0 #if this is zero, it'll do equal proportions
MAINALLELEPROP <- 0 #if this is zero, it'll do equal proportions
NUMREADS <- 100 #if this is zero, it'll report one read for each haplotype in the clutch
REPS <- 500

founderPop = runMacs(nInd=NUMSIRES+1, nChr=20, segSites=32)
SP = SimParam$new(founderPop)

#methods for setting and reading multiallelic loci
setMicrosatHaplos <- function(founderPop, richness, mainalleleprop){
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
takeSample <- function(Clutch, numReads){
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
getAlleleProps <- function(Reads, allelicRichness){
  output <- data.frame(matrix(0, nrow = allelicRichness, ncol = ncol(Reads)))
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
getNullAlleleDists <- function(numLoci, allelicRichness, mainAlleleProp,seqDepth,reps){
  output <- array(0, dim = c(reps,allelicRichness,numLoci))
  for(rep in 1:dim(output)[1]){
    #so each rep, you make a new pair of parents and make a clutch of full-sibs
    founderPop = runMacs(nInd = 2, nChr = numLoci, segSites=32)
    SP = SimParam$new(founderPop)
    parents <- newPop(founderPop)
    parents <- setMicrosatHaplos(parents, allelicRichness, mainAlleleProp)
    clutchSize <- sample(c(300:1000),1)
    Clutch <- makeCross(pop = parents, crossPlan = matrix(c(1, 2), ncol = 2), 
              nProgeny = clutchSize)
    #then you take reads from that clutch at some sequencing depth, and get the 
    #proportions of each allele
    fullGeno <- takeSample(Clutch,seqDepth)
    genoProps <- getAlleleProps(fullGeno, allelicRichness)
    #then copy that whole guy into the output array
    for(allele in 1:dim(output)[2]){
      for(locus in 1:dim(output)[3]){
        output[rep,allele,locus] <- genoProps[allele,locus]
      }
    }
  }
  return(output)
  }


#testing
founderPop = runMacs(nInd = 2, nChr = NUMLOCI, segSites=32)
SP = SimParam$new(founderPop)
parents <- newPop(founderPop)
parents <- setMicrosatHaplos(parents, ALLELESPERLOCUS, MAINALLELEPROP)
clutchSize <- sample(c(300:1000),1)
Clutch <- makeClutch(parents, clutchSize, PROPMAINSIRE)
scoreMicrosats(Clutch)
Reads <- takeSample(Clutch, NUMREADS)
dists <- getAlleleProps(Reads, ALLELESPERLOCUS)
nullAlleleDists <- getNullAlleleDists(NUMLOCI, ALLELESPERLOCUS, 0,100,100)
#now we've got a 3 dimensional array, showing the relative proportions
#of each allele in a bunch of sample pool sequences of FULL SIBS

#we can get distributions of each allele at each locus now
hist(nullAlleleDists[,1,1], breaks = 100)



