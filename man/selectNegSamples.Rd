\name{selectNegSamples}
\alias{selectNegSamples}



\title{
  selecting negative samples for training RafSee
}


\description{
  Negative samples are selected based on their connectivity with positive samples in protein-protein interaction network. Of note, the user can also have a try to randomly selected negative samples with a given number.
}


\usage{
  selectNegSamples(positives, PPIMat = NULL, balanced = TRUE,
                   ratio = 1, GenomeGeneID)
}

\arguments{
  \item{positives}{A vector of positive samples.}
  \item{PPIMat}{A matrix of PPI, which contains 3 coloums represent protein1, protein2, score respectively;Of note, negative samples will be selected randomly if this parameter is not assigned.}
  \item{balanced}{A logical value, where TRUE represents balance the positive and negative samples according to the ratio.}
  \item{ratio}{A numuric value of the the ratio between negative and positive samples.}
  \item{GenomeGeneID}{A vector of genome ID.}
}



\value{
  A vector of selected negative samples
}


\author{
  Jingjing Zhai, Chuang Ma
}


\examples{
  \dontrun{
    positives <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610")
    GenomeID <- c("AT1G01060", "AT1G09530", "AT1G09570", "AT1G12610", "AT1G77300", "AT1G79730")
    negatives <- selectNegSamples(positives = positives, PPIMat = PPIMat, GenomeID = GenomeID)
  }
}



\keyword{ selecting negative samples }

