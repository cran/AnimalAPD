#' Camera Trap Observations of African golden wolves
#' 
#' @name wolfexample
#' @docType data
#' @description Example dataset with wolf observation times in radians and the sampling period during which each observation was recorded
#' 
#' @format Dataframe with 2 columns and 30 rows
#'   Radians Time of observations, in radians from 0 to 2pi
#'   SamplingPeriod variable identifying camera trap sampling period
#'
#' 
#' @source \ Campbell L.A.D. 2017
#' 
#' @keywords datasets
#' 
"wolfexample"


#' Camera Trap Observations of wild boar
#' 
#' @name boarexample
#' @docType data
#' @description Example dataset with boar observation times in radians and the sampling period during which each observation was recorded
#' 
#' @format Dataframe with 2 columns and 35 rows
#'   Radians Time of observations, in radians from 0 to 2pi
#'   SamplingPeriod Variable identifying camera trap sampling period
#' 
#' 
#' @source \ Campbell L.A.D. 2017
#' 
#' @examples
#' \dontrun{ APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians, RE1=wolfexample$SamplingPeriod,
#'     weibullGLMM=FALSE,frechetGLMM=FALSE,gammaGLMM=FALSE,
#'     lognormalGLMM=FALSE,invgaussianGLMM=FALSE,
#'     mean=FALSE, HDI=FALSE, rawmean=TRUE)}
#' 
"boarexample"