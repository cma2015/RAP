
## select negative samples which are not connected to any positive samples in protein-protein-interaction(PPI) network
#' @export
selectNegSamples <- function(positives, PPIMat = NULL, balanced = TRUE, ratio = 1, GenomeGeneID){
  
  rownames(PPIMat) <- PPIMat[,1]
  GeneID <- substr(x = rownames(PPIMat), start = 6, stop = 14)
  interGeneID <- intersect(GeneID, positives) 
  
  if( !is.null(PPIMat) ){
    SubLinkMat <- matrix(NA, nrow = 1, ncol = ncol(PPIMat))
    for(i in 1:length(interGeneID)){
      curGeneID <- interGeneID[i]
      tmpMat <- PPIMat[grep(curGeneID, rownames(PPIMat)),]
      SubLinkMat <- rbind(SubLinkMat, tmpMat)
    }
    SubLinkMat <- SubLinkMat[-1,]
    
    nonNegatives <- substr(SubLinkMat[,2], start = 6, stop = 14)
    nonNegatives <- unique(nonNegatives)
    
    negatives <- setdiff(GenomeGeneID, nonNegatives)
    ##excluding positive samples which are not included in network
    negatives <- setdiff(negatives, positives)
    
    
  }else{
    cat("negative samples are selected randomly because PPIMat are not assigned!")
    negatives <- setdiff(GenomeGeneID, positives)
    posLen <- length(positives)
    negLen <- length(negatives)
  }
  
  
  if( balanced == TRUE ){
    if( negLen >= posLen*ratio){
      negatives <- sample(negatives, size = posLen*ratio)
    }
  }
  
  negatives
} 
