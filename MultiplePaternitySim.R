#RB Smith
#BSF MSc Project
#Multiple Paternity Simulation Script

install.packages("AlphaSimR")
library(AlphaSimR)

rm(list = ls())

#Example making a population (with SNPs) and simulating some crosses

founderGenomes = runMacs(nInd = 4, nChr = 7, segSites = 1000)
SP = SimParam$new(founderGenomes)

basePop = newPop(founderGenomes)

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








#Generalized functions


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
  ClutchGeno <- pullSegSiteGeno(makeCross(pop = basePop, crossPlan = matrix(c(1, 2), ncol = 2), nProgeny = sire1Brood))
  for(i in 3:(numSires+1)){
    altcross <- pullSegSiteGeno(makeCross(pop = basePop, crossPlan = matrix(c(1, i), ncol = 2), nProgeny = altSireBrood))
    ClutchGeno <- rbind(ClutchGeno, altcross)
  }
  return(ClutchGeno)
}


#Sampling function, which takes a certain number of random offspring 
#from a dataframe of offspring genotypes
takeSample <- function(ClutchGeno, sampNum){
  return(ClutchGeno[,sample(ncol(ClutchGeno), sampNum)])
}

#CSD locus function taken directly from SIMply Bee

# Edit genome and controlled mating ----
#' @title Edit the csd locus
#'
#' @description Edits the csd locus in an entire population of individuals to
#'   ensure heterozygosity. The user can provide a list of csd alleles for each
#'   individual or, alternatively, the function samples a heterozygous genotype
#'   for each individual from all possible csd alleles. The gv slot is
#'   recalculated to reflect the any changes due to editing, but other slots
#'   remain the same.
#'
#' @param pop \code{\link{Pop-class}}
#' @param alleles \code{NULL} or list;
#'   If \code{NULL}, then the function samples a heterozygous csd genotype for
#'   each virgin queen from all possible csd alleles.
#'   If not \code{NULL}, the user provides a list of length \code{nInd} with each
#'   node holding a matrix or a data.frame, each having two rows and n columns.
#'   Each row must hold one csd haplotype (allele) that will be assigned to a
#'   virgin queen. The n columns span the length of the csd locus as specified
#'   in \code{\link{SimParamBee}}. The two csd alleles must be different to
#'   ensure heterozygosity at the csd locus.
#' @param simParamBee global simulation parameters.
#'
#' @return Returns an object of \code{\link{Pop-class}}
editCsdLocus <- function(pop, alleles = NULL, simParamBee = NULL) {
  if (is.null(simParamBee)) {
    simParamBee <- get(x = "SP", envir = .GlobalEnv)
  }
  csdSites <- simParamBee$csdPosStart:simParamBee$csdPosStop
  if (is.null(alleles)) {
    # Create a matrix of all possible csd alleles
    alleles <- expand.grid(as.data.frame(matrix(rep(0:1, length(csdSites)), nrow = 2, byrow = FALSE)))
    # Sample two different alleles (without replacement) for each individual
    nAlleles <- simParamBee$nCsdAlleles
    alleles <- sapply(seq_len(pop@nInd), FUN = function(x) list(alleles[sample(nAlleles, size = 2, replace = F), ]))
  }
  
  if (pop@nInd != length(alleles)) {
    stop("The length of the allele list must match the number of individuals in the pop argument")
  }
  if (any(sapply(alleles, FUN = function(x) all(x[1, ] == x[2, ])))) {
    stop("You must provide two different alleles for each individual!")
  }
  
  # Pull out csd haplotype matrix
  csdH = pullMarkerHaplo(pop, markers = paste(simParamBee$csdChr, csdSites, sep="_"))
  # Prepare the haplotype matrix
  alleles <- as.matrix(do.call(rbind, alleles))
  rownames(alleles) <- rownames(csdH)
  colnames(alleles) <- colnames(csdH)
  
  pop <- setMarkerHaplo(pop, haplo=alleles)
  return(pop)
}


