\name{FeatureExtract}
\alias{FeatureExtract}



\title{
  Extracting sequence-based features.
}


\description{
  This function generates sequence-based features from protein sequences using five scoring schemes.
}


\usage{
  FeatureExtract(ProteinSeq, feature = c("AAC", "DAAC", "PAAC", "APAAC", "PCP"), 
                 lambda = 5, w = 0.05)
}

\arguments{
  \item{ProteinSeq}{A list of protein sequences.}
  \item{feature}{A vector of encoding schemes.}
  \item{lambda}{The lambda parameter for the PAAC and APAAC-related features, default is 5.}
  \item{w}{The weighting parameter for the PAAC and APAAC-related features, default is 0.05.}
}



\value{
  A feature matrix with genes in rows, features in columns
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    
    ##generate a list of protein sequence
    exampleSeq1 <- "MVQYNFKRITVVPNGKEFVDIILSRTQRQTPTVVHKGYKINRLRQFYMRKVKYTQTNFHAKLSAIIDEFP"
    exampleSeq2 <- "MDSESESKLISFISQLVSRNNTDSENISCMIQTISLVSSMDLKSQPKPESKLMSLVTQTISLFNSM"
    featureMat <- FeatureExtract(ProteinSeq = list(exampleSeq1, exampleSeq2), feature = "AAC")
  }
}



\keyword{ feature extract}

