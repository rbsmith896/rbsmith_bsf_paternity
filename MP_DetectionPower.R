#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

#alrighty here we go boys
#so this guy is gonna take:

#number of loci
#alleles per locus
#number of sires
#main sire proportion
#main allele frequency
#sample size

#and it'll make a bunch of clutches,
#take samples from each, and see if it
#sees multiple paternity 
#then it reports the propotion of 
#clutches in which multiple paternity was
#observed

NUMLOCI <- 15
ALLELESPERLOCUS <- 3
NUMSIRES <- 3
PROPMAINSIRE <- .8 #if this is zero, it'll do equal proportions
MAINALLELEPROP <- .8 #if this is zero, it'll do equal proportions
SAMPLESIZE <- 10 #if this is zero, it'll do the whole clutch
REPS <- 500
KINSHIPPROPCUTOFF <- .9 #this is the confidence with which we want to 
#be able to call half sibs from full sibs - the more stringent this 
#is (closer to 1), the lower our Type I error (less likely to accidentally)
#call full-sibs as half-sibs, but also the lower our power (more likely
#to miss true half-sibs)
#if it's 1, it'll make the cutoff the lowest full-sib relatedness

#first, we need a kinship cutoff value
founderPop = runMacs(nInd=4, nChr=20, segSites=32)
SP = SimParam$new(founderPop)

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
  parents <- runMacs(nInd=2, nChr=numLoci, segSites=32)
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
  parents <- runMacs(nInd=3, nChr=numLoci, segSites=32)
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

#making the distribution and getting the kinship cutoff
dist <- getKinshipDists(NUMLOCI,ALLELESPERLOCUS,MAINALLELEPROP, 500)
FullSibs <- sort(dist$FS)
if(KINSHIPPROPCUTOFF == 1){
  CUTOFF <- FullSibs[1]
}else{
  CUTOFF <- FullSibs[500*(1-KINSHIPPROPCUTOFF)]
}

scoreCard <- rep(0,REPS)

for(rep in 1:REPS){
  founderPop = runMacs(nInd = NUMSIRES + 1, nChr = NUMLOCI, segSites=32)
  SP = SimParam$new(founderPop)
  parents <- newPop(founderPop)
  parents <- setMicrosatHaplos(parents, ALLELESPERLOCUS, MAINALLELEPROP)
  #assuming clutches are between 300 and 1000 eggs
  clutchSize <- sample(c(300:1000),1)
  Clutch <- makeClutch(parents, clutchSize, PROPMAINSIRE)
  
  Sample <- takeSample(Clutch, SAMPLESIZE)
  
  SampleKinships <- getPairwiseRelationships(Sample, Sample)
  MP <- 0
  for(i in 1:(SAMPLESIZE-1)){
    for(j in (i+1):SAMPLESIZE){
      if(SampleKinships[i,j] < CUTOFF){
        MP <- 1
      }
    }
   }
  scoreCard[rep] <- MP
  
}

MPprop <- sum(scoreCard)/length(scoreCard)
MPprop

#problem: there's no happy cutoff. Either it's too stringent and can't 
#observe multiple paternity, or it's too loose and calls multiple paternity
#where there is none

#this is fixed a little bit when the sample size is big, and the
#cutoff is extremely stringent - we can set the cutoff so that it NEVER
#accidentally calls full-sibs as half-sibs, and then it's not great at 
#calling true half-sibs, but it can get them sometimes with a large enough
#sample size


#maybe another way to think about it is to call MP when there's an individual
#that's called as a half sib with more than half the remainders...
#that would look a little like this

scoreCard <- rep(0,REPS)
for(rep in 1:REPS){
  founderPop = runMacs(nInd = NUMSIRES + 1, nChr = NUMLOCI, segSites=32)
  SP = SimParam$new(founderPop)
  parents <- newPop(founderPop)
  parents <- setMicrosatHaplos(parents, ALLELESPERLOCUS, MAINALLELEPROP)
  #assuming clutches are between 300 and 1000 eggs
  clutchSize <- sample(c(300:1000),1)
  Clutch <- makeClutch(parents, clutchSize, PROPMAINSIRE)
  
  Sample <- takeSample(Clutch, SAMPLESIZE)
  
  SampleKinships <- getPairwiseRelationships(Sample, Sample)
  MPcall <- 0
  for(i in 1:SAMPLESIZE){
    MP <- rep(0,SAMPLESIZE)
    for(j in 1:SAMPLESIZE){
      if(SampleKinships[i,j] < CUTOFF){
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
MPprop

#alrighty, I think it actually works like this! With a decent sized sample,
#and a cutoff around .9, we see mostly positive MP calls when there's
#MP, and mostly negative MP calls when there's not

#A little taster data
#with 15 loci, 3 alleles per locus
#80% contribution by main sire, 0.8 frequency for main allele
#0.9 kinship cutoff for calling full sibs
#50 reps

#1 sire, samp size 10:   0.14 called for MP
#1 sire, samp size 25:   0.08 called for MP
#1 sire, samp size 50:   0.04 called for MP
#1 sire, samp size 100:  0.04 called for MP

#2 sires, samp size 10:  0.30 called for MP
#2 sires, samp size 25:  0.54 called for MP
#2 sires, samp size 50:  0.74 called for MP
#2 sires, samp size 100: 0.74 called for MP

#3 sires, samp size 10:  0.38 called for MP
#3 sires, samp size 25:  0.58 called for MP
#3 sires, samp size 50:  0.76 called for MP
#3 sires, samp size 100: 0.84 called for MP

#probability that you get multple in fathers in the SAMPLE
#.     - should add in this part
#      - differentiate the accuracy of the method generally, vs a particular setting
#probability that your "call" is correct

#for each parameter, find one level to be the "default", "base"
#vary one at a time, keeping others constant
#MAYBE vary 2 at a time
#think about plots, sketch them out

#fix cutoff at .9 for now, think about how they look
#maybe redo them more and less stringent, see what the outcome is

#changed to equal proportions of sires and alleles...
#(gotta remember to rerun the cutoff any time the genetics change)

#1 sire, samp size 10:   0.04 called for MP
#1 sire, samp size 25:   0.00 called for MP
#1 sire, samp size 50:   0.02 called for MP
#1 sire, samp size 100:  0.00 called for MP

#2 sires, samp size 10:  0.50 called for MP
#2 sires, samp size 25:  0.70 called for MP
#2 sires, samp size 50:  0.68 called for MP
#2 sires, samp size 100: 0.76 called for MP

#3 sires, samp size 10:  0.78 called for MP
#3 sires, samp size 25:  0.90 called for MP
#3 sires, samp size 50:  0.94 called for MP
#3 sires, samp size 100: 0.96 called for MP

SampleKinshipPlot <- SampleKinships %>% rownames_to_column() %>% gather(colname, value, -rowname)
ggplot(SampleKinshipPlot, aes(x = rowname, y = colname, fill = value)) + geom_tile()


knownKinships <- SampleKinships
fathers <- strtoi(Sample@father)
for(i in 1:nrow(SampleKinships)){
  for(j in 1:ncol(SampleKinships)){
    dad1 <- fathers[i]
    dad2 <- fathers[j]
    if(i == j){
      knownKinships[i,j] <- 0
    }
    else if(dad1 == dad2){
      knownKinships[i,j] <- 1
    }
    else{
      knownKinships[i,j] <- .5
    }
  }
}
knownKinships


install.packages("sequoia")
library(sequoia)

save.image("/home/s2426533/MP_Detection_Test_Script.RData")
