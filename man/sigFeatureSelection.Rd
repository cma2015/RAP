\name{sigFeatureSelection}
\alias{sigFeatureSelection}



\title{
  Selecting informative features
}


\description{
  This function extracting informative features with feature selection algorithms including the student's t-test and chi-square test feature selection algorithms.
}


\usage{
  sigFeatureSelection(featureMatrix, positives, negatives, binary = FALSE, level = 0.05)
}

\arguments{
  \item{featureMatrix}{A numeric matrix of features where rows represent genes, cols represent features}
  \item{positives}{A character vector of positive samples}
  \item{negatives}{A character vector of negative samples}
  \item{binary}{A logical value, where TRUE represents the features are binary with 0 and 1, the default is FALSE}
  \item{level}{A numeric value recording the significant level, the default is 0.05}
}



\value{
  A numeric feature matrix with only significant features were contained
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    positives <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610")
    sifFeatureMat <- sigFeatureSelection(featureMatrix = featureMat, 
                                   positives = positives, 
                                   negatives = negatives)
    ## featureMat can be calculated by function FeatureExtract
    ## negatives can be calculated by function selectNegSamples
  }
}



\keyword{ selecting significant features }

