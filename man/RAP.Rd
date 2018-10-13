\name{RAP}
\alias{RAP}



\title{
  Rank Aggregation-based Data Fusion for Gene Prioritization
}


\description{
  This function prioritize genes using sequence-based and network-based model
}


\usage{
  RAP(netPredResult, positives, negatives = NULL, featureSel = TRUE,
      featureMat = NULL, ProteinSeq,  PPIMat, GenomeGeneID, 
      ntree = 500)
}

\arguments{
  \item{netPredResult}{Full path of gene prioritization results from the network-based gene prioritization algorithms (e.g., AraNet v2).}
  \item{featureMat}{A numeric matrix of features where rows represent genes, cols represent features.}
  \item{positives}{A character vector of positive samples.}
  \item{negatives}{A character vector of negative samples.}
  \item{featureSel}{A logical value, where TRUE selecting significant features between positive and negative samples.}
  \item{ProteinSeq}{A list of protein sequence, the parameter are not required unless patameter featureMat is NULL.}
  \item{PPIMat}{A matrix of PPI, which contains 3 coloums represent protein1, protein2, score respectively, the patameter is not required unless paramter negatives are NULL.}
  \item{GenomeGeneID}{A vector of genome ID, the parameter, the parameter is not required unless parameter negatives is NULL.}
  \item{ntree}{Number of trees to grow when using random forest, the default is 500.}
}



\value{
  A matrix with the rank of genes and descriptive information
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    positives <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610")
    res <- RAP(netPredResult = "/home/malab/AraNetPred.txt", 
               positives = positives, negatives = negatives, 
               featureSel = TRUE, featureMat = featureMat)
    ## featureMat can be calculated by function FeatureExtract
    ## negatives can be calculated by function selectNegSamples
    ## The sample results of AraNet v2 (AraNetPred.txt) can be
       downloaded from http://bioinfo.nwafu.edu.cn/software
  }
}



\keyword{ gene prioritization }

