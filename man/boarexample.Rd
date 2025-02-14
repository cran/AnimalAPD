% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{boarexample}
\alias{boarexample}
\title{Camera Trap Observations of wild boar}
\format{
Dataframe with 2 columns and 35 rows
Radians Time of observations, in radians from 0 to 2pi
SamplingPeriod Variable identifying camera trap sampling period
}
\source{
\ Campbell L.A.D. 2017
}
\usage{
boarexample
}
\description{
Example dataset with boar observation times in radians and the sampling period during which each observation was recorded
}
\examples{
\dontrun{ APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians, RE1=wolfexample$SamplingPeriod,
    weibullGLMM=FALSE,frechetGLMM=FALSE,gammaGLMM=FALSE,
    lognormalGLMM=FALSE,invgaussianGLMM=FALSE,
    mean=FALSE, HDI=FALSE, rawmean=TRUE)}

}
\keyword{datasets}
