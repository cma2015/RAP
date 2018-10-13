#' @import RobustRankAggreg
#' @import missForest
#' @importFrom stats predict t.test
#' @importFrom stats predict chisq.test
#' @importFrom utils read.table
#' @export
RAP <- function(netPredResult , positives, negatives = NULL, featureSel = TRUE, featureMat = NULL,
                ProteinSeq,  PPIMat, GenomeGeneID, ntree = 500){
  
  AraNetPred <- as.matrix(read.table(file = netPredResult , sep = "\t", header = F, quote = ""))
  rownames(AraNetPred) <- AraNetPred[,2]
  if( is.null(negatives) ){
    cat("start selecting negative samples...", "\n")
    negatives <- selectNegSamples(positives = positives, PPIMat = PPIMat, GenomeGeneID = GenomeGeneID)
  }
  
  if( is.null(featureMat) ){
    cat("start extracting sequenced based features...", "\n")
    featureMat <- FeatureExtract(ProteinSeq = ProteinSeq)  
  }
  
  if( featureSel ){
    featureMat <- sigFeatureSelection(featureMatrix = featureMat, positives = positives, negatives = negatives, level = 0.05)  
  }
  
  if( length( intersect( positives, rownames(featureMat) ) ) != length( positives ) ){
    stop("positives must be the subset of the row name of featureMatrix !")
  }
  
  if( length( intersect( negatives, rownames(featureMat) ) ) != length( negatives ) ){
    stop("negatives must be the subset of the row name of featureMatrix !")
  }
  #write.table(featureMat, file = "/home/malab8/zjj/FTgenePrediction/FeatureMatrix/test.txt", sep = "\t", quote = F)
  miss.featureMat <- missForest(xmis = featureMat)
  featureMat <- miss.featureMat$ximp
  
  posLen <- length(positives)
  negLen <- length(negatives)
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positives, negatives), ] )
  cat("start modeling using random forest...")
  obj <- randomForest( x = fmat, y = factor(label), ntree = ntree)
  predSample <- setdiff( rownames(featureMat), positives )
  predict.score <- predict(obj, featureMat[c( predSample ), ], type = "vote")[,"1"]
  predict.score <- predict.score[order( as.numeric(predict.score), decreasing = T )]
  
  GeneRF <- names(predict.score)
  GeneAraNet <- intersect( AraNetPred[,2], GeneRF )
  RRAMatrix <- aggregateRanks(glist = list(GeneAraNet, GeneRF))
  RRAMatrix <- as.matrix(RRAMatrix)
  rank_RRA <- 1:nrow(RRAMatrix)
  
  
  load(system.file('data/description.rda', package = 'RAP'))
  rownames(description) <- description[,1]
  resMat <- matrix(NA, nrow = length(predict.score), ncol = 6)
  rownames(resMat) <- RRAMatrix[,1]
  colnames(resMat) <- c("sequence-based rank", "network-based rank", "RAP-based rank", colnames(description)[2:4])
  resMat[names(predict.score), 1] <- 1:length(predict.score) 
  resMat[intersect( rownames(resMat), AraNetPred[,2] ), 2] <- AraNetPred[intersect( rownames(resMat), AraNetPred[,2] ),1] 
  resMat[,3] <- rank_RRA
  resMat[,4:6] <- description[rownames(resMat), 2:4]
  
  resMat
}
