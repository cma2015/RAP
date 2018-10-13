\name{CrossValidation}
\alias{CrossValidation}



\title{
  Cross validation
}


\description{
  Performance evaluation of the integrative random forest-based gene prioritization algorithm RafSee in distinguishing positives and negatives.
}


\usage{
  CrossValidation(seed = 1, featureMat, positives, 
                  negatives, cross = 10, cpus = 1)
}

                  

\arguments{
  \item{seed}{An integer number specifying a random seed for randomly partitioning dataset.}
  \item{featureMat}{A numeric feature matrix.}
  \item{positives}{A character vector representing positive samples.}
  \item{negatives}{A character vector representing negative samples.}
  \item{cross}{Number of fold for cross validation.}
  \item{cpus}{An integer number specifying the number of cpus to be used for parallel computing.}
}



\value{
  A list containing results from each fold cross validation including:
  \item{positives.train }{positive samples used to train prediction model.}
  \item{negatives.train }{negative samples used to train prediction model.}
  \item{positives.test }{positive samples used to test prediction model.}
  \item{negatives.test }{negative samples used to test prediction model.}
  \item{positives.train.score }{scores of positive samples in training dataset predicted by random foreste.}
  \item{negatives.train.score }{scores of negative samples in training dataset predicted by random forest.}
  \item{positives.test.score }{scores of positive samples in testing dataset predicted by random forest.}
  \item{negatives.test.score }{scores of negative samples in testing dataset predicted by random forest.}
  \item{train.AUC }{AUC value of random forest on training dataset.}
  \item{test.AUC }{AUC value of random forest on testing dataset.s}
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    positives <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610")
    cvRes <- CrossValidation(featureMat = featureMat, positives = positives, 
                      negatives  = negatives, cpus = 1)
    ## featureMat can be calculated by function FeatureExtract
    ## negatives can be calculated by function selectNegSamples
  }
}



\keyword{cross validation}

