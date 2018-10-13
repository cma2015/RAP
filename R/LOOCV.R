
#' @import randomForest
#' @import snowfall

.leave_one_out <- function(Len, featureMat, positives, negatives, predictSample){
  i <- Len
  curGeneID <- positives[i]
  train_positives <- positives[-i]
  
  posLen <- length(train_positives)
  negLen <- length(negatives)
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(train_positives, negatives), ] )
  obj <- randomForest(x = fmat, y = factor(label))
  res.predict.score <- predict(obj, featureMat[c(curGeneID, predictSample), ], type = "vote")[,"1"]
}




#' @export
LOOCV <- function(featureMat = NULL, positives, negatives, cpus = 1, predictSample = NULL){
  
  
    positives <- intersect(positives, rownames(featureMat))
    negatives <- intersect(negatives, rownames(featureMat))
  
    if( length( positives ) < 10 ){
      cat("final positives:", length(positives), "\n")
      stop("too few positives!")
    }
    
    if( length( negatives ) < 10 ){
      cat("final negatives:", length(negatives), "\n")
      stop("too few negatives!")
    }
    
  if( is.null(predictSample) ){
    predictSample <- setdiff(rownames(featureMat), positives)
  }

  
  Len <- length(positives)
  loocvRes <- list()
  
  if(cpus > 1){
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".leave_one_out")
    sfLibrary( "randomForest", character.only = TRUE )
    loocvRes <- sfApply(matrix(1:Len, ncol = 1), 1, .leave_one_out, 
                        positives = positives, negatives = negatives, 
                        featureMat = featureMat, predictSample = predictSample)
    rownames(loocvRes) <- 1:nrow(loocvRes)
    sfStop()
  }else{
    for(i in 1:Len){
      loocvRes[[i]] <- .leave_one_out(Len = i, featureMat = featureMat, positives = positives, negatives = negatives, predictSample = predictSample)
    }
  }
  loocvRes
}
