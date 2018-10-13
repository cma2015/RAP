## extract physical chemical properties
#' @export
extractPCP <- function(ProteinSequence, AAindex){
  Len <- nchar(ProteinSequence)
  PhyChemicalFea <- matrix(NA, nrow = 1, ncol = nrow(AAindex))
  for(i in 1:nrow(AAindex)){
    kk <- table(strsplit(ProteinSequence, ""))
    
    tmpMat <- matrix(NA, nrow = length(kk), 2)
    rownames(tmpMat) <- names(kk)
    tmpMat[,1] <- as.numeric(kk)
    tmpMat[,2] <- as.numeric(AAindex[i, names(kk)])
    
    Pmax <- max(tmpMat[,2])
    Pmin <- min(tmpMat[,2])
    
    tmpMat[,2] <- (tmpMat[,2] - Pmin)/(Pmax - Pmin)
    
    Ptot <- sum(tmpMat[,1]*tmpMat[,2])/Len
    PhyChemicalFea[1,i] <- Ptot
  }
  colnames(PhyChemicalFea) <- rownames(AAindex)
  PhyChemicalFea
}




#' @import protr
#' @export
FeatureExtract <- function(ProteinSeq, feature = c("AAC", "DAAC", "PAAC", "APAAC", "PCP"), lambda = 5, w = 0.05){
  #ProteinSeq <- readFASTA(SeqDir)
  featureMat <- list()
  
  ## extract amino acid composition
  if( is.element("AAC", feature) ){
    cat("start computing sequence-based feature AAC", "\n")
    featureAAC <- t( sapply(ProteinSeq, extractAAC) )
    rownames(featureAAC) <- substr(x = rownames(featureAAC), start = 1, stop = 9)
    featureMat[["AAC"]] <- featureAAC
  }
  
  
  ## extract dipeptide composition
  if( is.element("DAAC", feature)){
    cat("start computing sequence-based feature DC", "\n")
    featureDC <- t( sapply(ProteinSeq, extractDC) )
    rownames(featureDC) <- substr(x = rownames(featureDC), start = 1, stop = 9)
    featureMat[["DAAC"]] <- featureDC 
  }
  
  
  ## extract amphiphilic pseudo amino acid composition
  if( is.element("APAAC", feature) ){
    cat("start computing sequence-based feature APAAC", "\n")
    featureAPAAC <- t( sapply(ProteinSeq, extractAPAAC, lambda = lambda, w = w) )
    rownames(featureAPAAC) <- substr(x = rownames(featureAPAAC), start = 1, stop = 9)  
    featureMat[["APAAC"]] <- featureAPAAC
  }
  
  
  ## extract pseudo amino acid composition
  if( is.element("PAAC", feature) ){
    cat("start computing sequence-based feature PAAC", "\n")
    featurePAAC <- t( sapply(ProteinSeq, extractPAAC, lambda = lambda, w = w) )
    rownames(featurePAAC) <- substr(x = rownames(featurePAAC), start = 1, stop = 9)
    featureMat[["PAAC"]] <- featurePAAC
  }
  
  ## extract physical chemical properties
  if( is.element("PCP", feature) ){
    cat("start computing sequence-based feature PCP", "\n")
    system.file('data/AAindex.rda', package = 'RAP')
    AAindex <- AAindex
    featurePCP <- t( sapply(ProteinSeq, extractPCP, AAindex) )
    rownames(featurePCP) <- names(ProteinSeq)
    colnames(featurePCP) <- rownames(AAindex)
    featureMat[["PCP"]] <- featurePCP
  }
  
  
  ##merge feature
  featureAll <- matrix(NA, nrow = nrow(featureMat[[1]]), ncol = 1)
  for(i in 1:length(featureMat)){
    featureAll <- cbind( featureAll, featureMat[[i]] )
  }
  featureAll <- featureAll[,-1]
  featureAll
}
