.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}

#' @import pROC
#' 
.one_cross_validation <- function( cv, featureMat, positives, negatives, posSample_cv, negSample_cv) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]

  
  posLen <- length(train_genes_p)
  negLen <- length(train_genes_n)
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(train_genes_p, train_genes_n), ] )
  obj <- randomForest(x = fmat, y = factor(label))
  
  positives.train.score <- predict(obj, featureMat[train_genes_p,], type = "vote")[,"1"]
  negatives.train.score <- predict(obj, featureMat[train_genes_n,], type = "vote")[,"1"]
  positives.test.score <- predict(obj, featureMat[test_genes_p,], type = "vote")[,"1"]
  negatives.test.score <- predict(obj, featureMat[test_genes_n,], type = "vote")[,"1"]
  
  
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n,
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}






#' @export
CrossValidation <- function( seed = 1, featureMat, positives, negatives, cross = 10, cpus = 1){
  
  call <- match.call()
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport(".one_cross_validation", namespace = "RAP")
    sfLibrary( "pROC", character.only=TRUE)
    sfLibrary( "randomForest", character.only=TRUE )
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv)
    }
  }
  cvRes
}

