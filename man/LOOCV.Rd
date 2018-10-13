\name{LOOCV}
\alias{LOOCV}



\title{
  leave-one-out cross-validation
}


\description{
  Leave-one-out cross-validation algorithm is performed to train and test the integrative random-forest gene prioritization algorithm RafSee.
}


\usage{
  LOOCV(featureMat, positives, negatives, cpus = 1, predictSample = NULL)
}

\arguments{
  \item{featureMat}{A numeric matrix of features where rows represent genes, cols represent features}
  \item{positives}{A character vector of positive samples}
  \item{negatives}{A character vector of negative samples}
  \item{cpus}{an integer number specifying the number of cpus to be used for parallel computing, the default is 1}
  \item{predictSample}{A vector of testing samples, if it is NULL, all genes excluding positive smaples were used}
}



\value{
  Predictive score for each leave-one-out cross-validation
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    positives <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610")
    loocvRes <- LOOCV(featureMat = featureMat, positives = positives, 
                      negatives  = negatives, cpus = 1)
    ## featureMat can be calculated by function FeatureExtract
    ## negatives can be calculated by function selectNegSamples
  }
}



\keyword{ machine learning }

