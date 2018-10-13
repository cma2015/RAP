### select significant features according to the student's t-test
#' @export
sigFeatureSelection <- function(featureMatrix, positives, negatives, binary = FALSE, level = 0.05){
  positives <- intersect( positives, rownames(featureMatrix) )
  if( length( positives ) < 10 ){
    cat("final positives:", length(positives), "\n")
    stop("too few positives!")
  }
  
  negatives <- intersect( negatives, rownames(featureMatrix) )
  if( length( negatives ) < 10 ){
    cat("final negatives:", length(negatives), "\n")
    stop("too few negatives!")
  }
  
  numFeature <- ncol(featureMatrix)
  FeatureSigMat <- matrix(NA, nrow = numFeature, ncol = 2)
  colnames(FeatureSigMat) <- c("feature", "p-value")
  
  if(binary){
    for(i in 1:numFeature){
      curFeature <- colnames(featureMatrix)[i]
      posFea <- as.numeric( featureMatrix[c(positives), i] )
      negFea <- as.numeric( featureMatrix[c(negatives), i] )
      chiMat <- matrix(NA, nrow = 2, ncol = 2)
      rownames(chiMat) <- c("positives", "negatives")
      colnames(chiMat) <- c("1", "0")
      pos1Len <- length(which(posFea == 1))
      pos0Len <- length(posFea) - pos1Len
      neg1Len <- length(which(negFea == 1))
      neg0Len <- length(negFea) - neg1Len
      chiMat["positives", ] <- c(pos1Len, pos0Len)
      chiMat["negatives", ] <- c(neg1Len, neg0Len)
      class(chiMat) <- "numeric"
      chisq.stat <- chisq.test(x = chiMat)
      FeatureSigMat[i, ] <- c(curFeature, chisq.stat$p.value)
    }
  }else{
    for(i in 1:numFeature){
      curFeature <- colnames(featureMatrix)[i]
      posFea <- as.numeric( featureMatrix[c(positives), i] )
      negFea <- as.numeric( featureMatrix[c(negatives), i] )
      t.stat <- t.test(x = posFea, y = negFea)
      FeatureSigMat[i, ] <- c(curFeature, t.stat$p.value)
    }
  }

  
  sigFeature <- FeatureSigMat[which( as.numeric( FeatureSigMat[,2] ) <= level ), 1]
  
  sigFeatureMat <- featureMatrix[, c(sigFeature)]
  
  sigFeatureMat
}



