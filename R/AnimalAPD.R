#' Activity Probability Density Controlling for Random Effects
#' 
#' @description Calculation of animal activity probability density controlling for nested data with random intercepts using Bayesian GLMMs with 'STAN' and \code{\link[brms]{brm}}. 
#'     The function can automatically select the statistical distribution that is most appropriate for the dataset (weibull, frechet, gamma, lognormal, inverse gaussian) 
#'     using \code{\link[loo]{loo}} and automatically ensures that MCMC chains converge and that a specified minimum effective sample size from the posterior distribution 
#'     is achieved. An APD activity curve plot is provided. 
#' 
#' Package: AnimalAPD
#' Version: 1.0.0
#' Date: 2020-11-08
#' 
#' @author Liz AD Campbell
#' @keywords cameratrap; activity; temporal; Bayesian
#' 
#' @import brms
#' @import circular
#' @import overlap
#' @import activityGCMM
#' 
#' 
#' @param focal Vector of observations in radians of one species/group/individual/etc. for which predictions on another will be made
#' @param contingent Vector of observations in radians or output from generalized circular mixture model of activity curves from \code{link[activityGCMM]{GCMM}} of a species/group/individual/etc. from which predictions will be made
#' @param RE1 Vector identifying a random intercept for observations of the focal to control for hierarchical data (e.g. camera trap IDs)
#' @param RE2 Optional vector identifying levels of a second random effect, for data with additional hierarchical levels (e.g. study sites, sampling periods, data collection seasons); default is NULL
#' @param weibullGLMM Specifies whether to run a weibull GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param frechetGLMM Specifies whether to run a frechet GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param gammaGLMM Specifies whether to run a Gamma GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param lognormalGLMM Specifies whether to run a lognormal GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param invgaussianGLMM Specifies whether to run a inverse.gaussian GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param cores Number of cores to use when running MCMC chains in parallel; default=1
#' @param iter Number of MCMC iteractions per chain; default=5000
#' @param burnin Number of MCMC iteractions discarded as burnin; default=iter/2
#' @param center Value to use as center of graph; default=pi
#' @param adapt_delta Value to use for adapt_delta with brms; default=0.95
#' @param minESS Desired minimum effective sample size; default=1000
#' @param Reloo Whether to use reloo when running leave-one-out cross-validation of models (loo)
#' @param col Specifies colour of points for the focal in the graph
#' @param histcol Specifies colour of the histogram plot of the posterior distribution
#' @param linecol Specifies colour of HDI line in histogram plot of the posterior distribution
#' @param xlimCV A vector of two values indicating the x axis limits for the histCV graph
#' @param min Whether to include minimum APD on the graph; default=TRUE; default=TRUE
#' @param max Whether to include maximum APD on the graph; default=TRUE; default=TRUE
#' @param points Whether to include datapoints for observations of the focal on the graph; default=TRUE
#' @param mean Whether to include the estimated mean APD from the GLMM on the graph; default=TRUE
#' @param HDI Whether to include the estimated 95% highest density interval of mean APD from the GLMM on the graph; default=TRUE
#' @param adjust Smoothing of predicted line; recommended to use default value for observed values and higher value for estimations from circular models 
#' @param rawmean Whether to include the raw mean, not correcting for random effects, on the graph; default=FALSE
#' @param ... Additional parameters
#' 
#' @return Prints graph of activity curve and APD estimates from best-fitting GLMM and prints summary of analysis. Returns object of class \code{APD} is returned, containing a list of analysis results and details: 
#' @return \code{data} List of data used in analysis
#' @return \code{output} Matrix with summary output from selected model
#' @return \code{distribution} Name of distribution of selected model
#' @return \code{model} An object of class \code{brmsfit} containing output from the selected model, including the posterior samples and other information. See \code{\link[brms]{brm}}
#' @return \code{CVposterior} Numeric vector of posterior samples for the calculated family-specific population coefficient of variation (CV)
#' @return \code{allmodels} List of objects of class \code{brmsfit} containing output from all models from the analysis.
#' @return \code{rawvalues} Numeric vector of the raw, uncorrected APD values
#' @return \code{rawsummary} List of summary stats of raw APD values
#' 
#' @seealso \code{\link[activityGCMM]{GCMM}} \code{\link[brms]{brm}} \code{\link[loo]{loo}}
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' \donttest{ APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians, 
#'     RE1=wolfexample$SamplingPeriod, weibullGLMM=TRUE, frechetGLMM=FALSE,
#'     gammaGLMM=FALSE, lognormalGLMM=FALSE, invgaussianGLMM=FALSE,
#'     min=TRUE, max=TRUE, points=TRUE, mean=TRUE, HDI=TRUE, rawmean=FALSE) }
#' 
#' @export
	APDRE<- function(focal, contingent, RE1, RE2=NULL,
		weibullGLMM=TRUE,frechetGLMM=TRUE,gammaGLMM=TRUE,lognormalGLMM=FALSE,invgaussianGLMM=TRUE,
		cores=1, iter=5000, burnin=iter/2, center="pi", Reloo=TRUE, adapt_delta=0.95, adjust=1, minESS=1000,
		col="deeppink3", histcol="deeppink", linecol="black", xlimCV=NULL, min=TRUE, max=TRUE, points=TRUE, mean=TRUE, HDI=TRUE,rawmean=FALSE,...) {

requireNamespace("overlap",quietly=TRUE)
requireNamespace("brms",quietly=FALSE)
requireNamespace("circular",quietly=TRUE)
requireNamespace("stats",quietly=TRUE)
requireNamespace("loo",quietly=TRUE)
	requireNamespace("activityGCMM",quietly=TRUE)

APDout<-NULL

	focal<-as.numeric(focal)
	if(class(contingent)=="GCMM"){ contingent<-as.numeric(contingent$GCMMmixture)
		adjust<-5 } else { contingent<-as.numeric(contingent) }
	RE1<-as.numeric(as.factor(RE1))
	adaptdelta<-adapt_delta
	if (center=="pi") { xcenter<-c("noon"); lim1<-0; lim2<-(2*pi) } else { xcenter<-c("midnight");lim1<-(-1*pi); lim2<-(pi) }
	if (length(RE2)>0) { RE2<-as.numeric(as.factor(RE2))
		datatemp<-data.frame(focal,RE1,RE2); colnames(datatemp)<-c("focal","RE1","RE2") } else {
		datatemp<-data.frame(focal,RE1); colnames(datatemp)<-c("focal","RE1") }
	
	APD<-stats::approxfun(densityPlot(as.numeric(contingent),adjust=adjust,xscale=NA,xcenter=xcenter,extend=NA,lwd=3,main="",ylab="Activity Probability Density",xaxs="i",xlim=c(lim1,lim2)))
		graphics::abline(v=c(-pi/2,-3*pi/2,-pi,0,pi,pi/2,3*pi/2),lty=2,col="gray60")
		if (rawmean==TRUE) { graphics::abline(h=mean(APD(focal)),lty=3,lwd=1,col=col)	}  
		if (min==TRUE) { graphics::abline(h=round(min(APD(focal)),2),lty=3,lwd=1,col=col)	}  
		if (max==TRUE) { graphics::abline(h=max(APD(focal)),lty=3,lwd=1,col=col)	} 
		if (points==TRUE) { points(x=focal,y=APD(focal),pch=16,cex=1.3); points(x=focal,y=APD(focal),pch=16,cex=.7,col=col)  } 
		graphics::box(lwd=2)

	Q1<-as.numeric(stats::quantile(APD(focal))[2])
	Q3<-as.numeric(stats::quantile(APD(focal))[4])
	qcv<-round( (Q3-Q1)/(Q1+Q3) ,2)
	IQR<-round((Q3-Q1),2)
	Range<-round( max(APD(focal)),2)-round(min(APD(focal)),2)

loolist<-c(-100000); loonames<-c("X"); loolistfull<-list()
## Weibull:
if (weibullGLMM==TRUE) {
	message("---------------------")
	message("Running Weibull model")
	message("  ")
	if (length(RE2)>0) { APDw<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=weibull,cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
		APDw<-brms::brm(APD(focal)~1+(1|RE1),data=datatemp,family=weibull,cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) }
		## if rhat > 1.05, will double the number of iterations up to 5 times
		i<-1	# counter for loop & iter
		if ( max(as.numeric(brms::rhat(APDw))) > 1.05) { repeat {
			i<-i+1 
			message(""); message("rhat > 1.05, increasing iterations..."); message("")	
			if ( max(as.numeric(brms::rhat(APDw))) < 1.05) { break }   
	  	      if ( i==5 ) { print("Maximum attempts reached for weibull model. Model stopped before all rhat < 1.05. Try increasing iter"); break }   
		if (length(RE2)>0) { APDw<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,warmup=burnin,family=weibull,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else {
			APDw<-brms::brm(APD(focal)~1+(1|RE1),family=weibull,cores=cores,data=datatemp,chains=3,warmup=burnin,iter=(iter*i),control=list(adapt_delta=adaptdelta)) }
			} }
		message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
		if ( max(as.numeric(rhat(APDw))) < 1.05) { 
			## if minimum ESS < selected minESS, will double the number of iterations up to 5 times
			i2<-1	# counter for loop & iter
			  wESS<-min(as.numeric(neff_ratio(APDw)))*(iter*(i+i2)/2*3)  
			if ( min(as.numeric(neff_ratio(APDw)))*(iter*i/2*3)<minESS ) { repeat {
				i2<-i2+1
				if (length(RE2)>0) { APDw<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,warmup=burnin,family=weibull,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else {
				  APDw<-brms::brm(APD(focal)~1+(1|RE1),cores=cores,family=weibull,data=datatemp,chains=3,warmup=burnin,iter=(iter*(i+i2)),control=list(adapt_delta=adaptdelta)) }
 			  wESS<-min(as.numeric(neff_ratio(APDw)))*(iter*(i+i2)/2*3)
					if ( wESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for weibull model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(wESS,0))); break }
				} }
		if ( min(as.numeric(neff_ratio(APDw)))*(iter*(i+i2)/2*3)>minESS ) {
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(wESS,0))); message("") } 
		}

	if ( max(as.numeric(rhat(APDw))) < 1.05) { if ( wESS>minESS ) {
	summ<-summary(APDw); RE1sdw<-round(c(summ$random$RE1[1,1],summ$random$RE1[1,3:4]),3);RE1sdw<-as.numeric(RE1sdw)
			if (length(RE2)>0) { RE2sdw<-round(c(summ$random$RE2[1,1],summ$random$RE2[1,3:4]),3);RE2sdw<-as.numeric(RE2sdw) }
	looAPDw<-brms::add_criterion(APDw, "loo", reloo=Reloo)
		loow<-loo::loo(looAPDw)$estimates[1,1]
		loolist<-c(loolist,loow); loonames<-c(loonames,"weibull")
		loolistfull<-c(loolistfull,looAPDw)
			message("Weibulll model loo:")
			print(round(loow,2))
			message(" ")
		} }
	}

## Frechet:	
if (frechetGLMM==TRUE) {
	message("---------------------")
	message("Running Frechet model")
	message("  ")
	if (length(RE2)>0) { APDf<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=frechet,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
	   APDf<-brms::brm(APD(focal)~1+(1|RE1),family=frechet,data=datatemp,cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) }
		## if rhat > 1.05, will double the number of iterations up to 5 times
		i<-1	# counter for loop & iter
		if ( max(as.numeric(brms::rhat(APDf))) > 1.05) { repeat {
			i<-i+1
			message(""); message("rhat > 1.05, increasing iterations..."); message("")	
			if (length(RE2)>0) { APDf<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=frechet,chains=3,iter=iter*i,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
			   APDf<-brms::brm(APD(focal)~1+(1|RE1),family=frechet,data=datatemp,cores=cores,chains=3,iter=(iter*i),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
				if ( max(as.numeric(brms::rhat(APDf))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for frechet model. Model stopped before all rhat < 1.05. Try increasing iter."); break }
			} }
		message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
		if ( max(as.numeric(brms::rhat(APDf))) < 1.05) {
			## if minimum ESS < minESS, will double the number of iterations up to 5 times
			i2<-1	# counter for loop & iter
				fESS<-min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)
			if ( min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)<minESS ) { repeat {
				i2<-i2+1
				if (length(RE2)>0) { APDf<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=frechet,chains=3,iter=iter*(i+i2),warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
				  APDf<-brms::brm(APD(focal)~1+(1|RE1),family=frechet,cores=cores,data=datatemp,chains=3,iter=(iter*(i+i2)),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				  fESS<-min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)
					if ( min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for frechet model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(fESS,0))); break }
				} }
		if ( min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)>minESS ) {
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(fESS,0))); message("") }  
		}

	if ( max(as.numeric(rhat(APDf))) < 1.05) { if ( fESS>minESS ) {
	summ<-summary(APDf); RE1sdf<-round(c(summ$random$RE1[1,1],summ$random$RE1[1,3:4]),3);RE1sdf<-as.numeric(RE1sdf)
			if (length(RE2)>0) { RE2sdf<-round(c(summ$random$RE2[1,1],summ$random$RE2[1,3:4]),3);RE2sdf<-as.numeric(RE2sdf) }
	looAPDf<-brms::add_criterion(APDf, "loo", reloo=Reloo)
		loof<-loo::loo(looAPDf)$estimates[1,1]
		loolist<-c(loolist,loof); loonames<-c(loonames,"frechet")
		loolistfull<-c(loolistfull,looAPDf)
			message("Frechet model loo:") 
			print(round(loof,2))  
			message(" ") 
		} }
	}
	 
## Gamma:
if (gammaGLMM==TRUE) {
	message("---------------------")
	message("Running Gamma model")
	message("  ")
 	if (length(RE2)>0) { APDg<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=stats::Gamma(link="log"),cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 		APDg<-brms::brm(APD(focal)~1+(1|RE1),data=datatemp,family=stats::Gamma(link="log"),cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 		## if rhat > 1.05, will double the number of iterations up to 5 times
 		i<-1	# counter for loop & iter
 		if ( max(as.numeric(brms::rhat(APDg))) > 1.05) { repeat {
 			i<-i+1 
			message(""); message("rhat > 1.05, increasing iterations..."); message("")	
	 		if (length(RE2)>0) { APDg<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=stats::Gamma(link="log"),cores=cores,chains=3,iter=iter*i,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 			APDg<-brms::brm(APD(focal)~1+(1|RE1),family=stats::Gamma(link="log"),cores=cores,data=datatemp,chains=3,iter=(iter*i),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				if ( max(as.numeric(brms::rhat(APDg))) < 1.05) { break }
 				if ( i==5 ) { print("Maximum attempts reached for Gamma model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
 			} }
		message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
 		if ( max(as.numeric(rhat(APDg))) < 1.05) { 
 			## if minimum ESS < minESS, will double the number of iterations up to 5 times
 			i2<-1	# counter for loop & iter
 			 gESS<-min(as.numeric(neff_ratio(APDg)))*(iter*(i+i2)/2*3)
 			if ( min(as.numeric(neff_ratio(APDg)))*(iter*i/2*3)<minESS ) { repeat {
 				i2<-i2+1
 				if (length(RE2)>0) { APDg<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=stats::Gamma(link="log"),chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else {
  				  APDg<-brms::brm(APD(focal)~1+(1|RE1),family=stats::Gamma(link="log"),data=datatemp,cores=cores,chains=3,iter=(iter*(i+i2)),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				  gESS<-min(as.numeric(neff_ratio(APDg)))*(iter*(i+i2)/2*3)
 					if ( gESS>minESS ) { break }
 					if ( i2==5 ) { message(paste("Maximum attempts reached for Gamma model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(gESS,0))); break }
 				} }
		if ( min(as.numeric(neff_ratio(APDg)))*(iter*(i+i2)/2*3)>minESS ) {
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(gESS,0))); message("") } 
 		}

	if ( max(as.numeric(rhat(APDg))) < 1.05) { if ( gESS>minESS ) {
 	summ<-summary(APDg)
	RE1sdg<-round(as.numeric(c(summ$random$RE1[1,1],summ$random$RE1[1,3:4]),3))
			if (length(RE2)>0) { RE2sdg<-round(c(summ$random$RE2[1,1],summ$random$RE2[1,3:4]),3);RE2sdg<-as.numeric(RE2sdg) }
	looAPDg<-brms::add_criterion(APDg, "loo", reloo=Reloo)
		loog<-loo::loo(looAPDg)$estimates[1,1]
		loolist<-c(loolist,loog); loonames<-c(loonames,"gamma")
		loolistfull<-c(loolistfull,looAPDg)
			message("Gamma model loo:")
			print(loog);message("") 
		} }
	}
	
#lognormal:
if (lognormalGLMM==TRUE) {
	message("---------------------")
	message("Running lognormal model")
	message("  ")
 	if (length(RE2)>0) { APDln<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=lognormal,cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 		APDln<-brms::brm(APD(focal)~1+(1|RE1),data=datatemp,family=lognormal,cores=cores,chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 		## if rhat > 1.05, will double the number of iterations up to 5 times
 		i<-1	# counter for loop & iter
 		if ( max(as.numeric(brms::rhat(APDln))) > 1.05) { repeat {
 			i<-i+1 
			message(""); message("rhat > 1.05, increasing iterations..."); message("")	
	 		if (length(RE2)>0) { APDln<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=lognormal,cores=cores,chains=3,iter=iter*i,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 			APDln<-brm(APD(focal)~1+(1|RE1),family=lognormal,cores=cores,data=datatemp,chains=3,iter=(iter*i),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				if ( max(as.numeric(brms::rhat(APDln))) < 1.05) { break }
 				if ( i==5 ) { print("Maximum attempts reached for lognormal model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
 			} }
		message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
 		if ( max(as.numeric(rhat(APDln))) < 1.05) {   
 			## if minimum ESS < minESS, will double the number of iterations up to 5 times
 			i2<-1	# counter for loop & iter
			  lnESS<-min(as.numeric(neff_ratio(APDln)))*(iter*(i+i2)/2*3)   
 			if ( min(as.numeric(neff_ratio(APDln)))*(iter*i/2*3)<minESS ) { repeat {
 				i2<-i2+1
 				if (length(RE2)>0) { APDln<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,family=lognormal,cores=cores,chains=3,iter=iter*(i+i2),warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
  				  APDln<-brms::brm(APD(focal)~1+(1|RE1),family=lognormal,cores=cores,data=datatemp,chains=3,iter=(iter*(i+i2)),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				   lnESS<-min(as.numeric(neff_ratio(APDln)))*(iter*(i+i2)/2*3)
 					if ( lnESS>minESS ) { break }
 					if ( i2==5 ) { message(paste("Maximum attempts reached for lognormal model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(lnESS,0))); break }
 				} }
		if ( min(as.numeric(neff_ratio(APDln)))*(iter*(i+i2)/2*3)>minESS ) {
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(lnESS,0))); message("") } 
 		}

	if ( max(as.numeric(rhat(APDln))) < 1.05) { if ( lnESS>minESS ) {
 	summ<-summary(APDln); RE1sdln<-round(c(summ$random$RE1[1,1],summ$random$RE1[1,3:4]),3);RE1sdln<-as.numeric(RE1sdln)
			if (length(RE2)>0) { RE2sdln<-round(c(summ$random$RE2[1,1],summ$random$RE2[1,3:4]),3);RE2sdln<-as.numeric(RE2sdln) }
	looAPDln<-brms::add_criterion(APDln, "loo", reloo=Reloo)
		looln<-loo(looAPDln)$estimates[1,1]
		loolist<-c(loolist,looln); loonames<-c(loonames,"lognormal")
		loolistfull<-c(loolistfull,looAPDln)
			message("Lognormal model loo:")
			print(looln); message("")
		} }
	}
	
#inverse gaussian: 
if (invgaussianGLMM==TRUE) {
	message("---------------------")
	message("Running inverse gaussian model")
	message("  ")
 	if (length(RE2)>0) { APDig<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=stats::inverse.gaussian(link="log"),chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 		APDig<-brms::brm(APD(focal)~1+(1|RE1),data=datatemp,cores=cores,family=stats::inverse.gaussian(link="log"),chains=3,iter=iter,warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 		## if rhat > 1.05, will double the number of iterations up to 5 times
 		i<-1	# counter for loop & iter
 		if ( max(as.numeric(brms::rhat(APDig))) > 1.05) { repeat {
 			i<-i+1 
			message(""); message("rhat > 1.05, increasing iterations..."); message("")	
	 		if (length(RE2)>0) { APDig<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=stats::inverse.gaussian(link="log"),chains=3,iter=iter*i,warmup=burnin,control=list(adapt_delta=adaptdelta)) } else {
 			APDig<-brms::brm(APD(focal)~1+(1|RE1),family=stats::inverse.gaussian(link="log"),data=datatemp,cores=cores,chains=3,iter=(iter*i),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 				if ( max(as.numeric(brms::rhat(APDig))) < 1.05) { break }
 				if ( i==5 ) { print("Maximum attempts reached for inverse Gaussian model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
 			} }
		message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
 		if ( max(as.numeric(rhat(APDig))) < 1.05) {  
 			## if minimum ESS < minESS, will double the number of iterations up to 5 times
 			i2<-1	# counter for loop & iter
			  igESS<-min(as.numeric(neff_ratio(APDig)))*(iter*(i+i2)/2*3)
 			if ( min(as.numeric(neff_ratio(APDig)))*(iter*(i+i2)/2*3)<minESS ) { repeat { 
 				i2<-i2+1
 				if (length(RE2)>0) { APDig<-brms::brm(APD(focal)~1+(1|RE1)+(1|RE2),data=datatemp,cores=cores,family=stats::inverse.gaussian(link="log"),chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else {
  				  APDig<-brms::brm(APD(focal)~1+(1|RE1),family=stats::inverse.gaussian(link="log"),data=datatemp,cores=cores,chains=3,iter=(iter*(i+i2)),warmup=burnin,control=list(adapt_delta=adaptdelta)) }
 	 			  igESS<-min(as.numeric(neff_ratio(APDig)))*(iter*(i+i2)/2*3)
 					if ( igESS>minESS ) { break }
 					if ( i2==5 ) { message(paste("Maximum attempts reached for inverse Gaussian model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(igESS,0))); break }
 				} }
		if ( igESS>minESS ) {
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(igESS,0))); message("") } 
 		}
	if ( max(as.numeric(rhat(APDig))) < 1.05) { if ( igESS>minESS ) {
	summ<-summary(APDig); RE1sdig<-round(c(summ$random$RE1[1,1],summ$random$RE1[1,3:4]),3); RE1sdig<-as.numeric(RE1sdig)
			if (length(RE2)>0) { RE2sdig<-round(c(summ$random$RE2[1,1],summ$random$RE2[1,3:4]),3);RE2sdig<-as.numeric(RE2sdig) }
	looAPDig<-brms::add_criterion(APDig, "loo", reloo=Reloo)
		looig<-loo(looAPDig)$estimates[1,1]
		loolist<-c(loolist,looig); loonames<-c(loonames,"invgaussian")
		loolistfull<-c(loolistfull,looAPDig)
			message("Inverse Gaussian model loo:")
			print(looig);message("")
		} }
	}

	loolow<-max(loolist)

if (weibullGLMM==TRUE) { 
	if ( max(as.numeric(rhat(APDw)))<1.05) { 
		if ( wESS>minESS ) {
			if (loolow==loow) { 
				
		best<-looAPDw; dist="weibull"
		InterceptMean<-fixef(looAPDw)[1]
			InterceptL<-fixef(looAPDw)[3]
			InterceptH<-fixef(looAPDw)[4]
		ShapeMean<-mean(as.data.frame(looAPDw)$shape)
			ShapeH<-stats::quantile((as.data.frame(looAPDw)$shape),.975)
			ShapeL<-stats::quantile((as.data.frame(looAPDw)$shape),.025)
		APDMean<-exp(InterceptMean)
			APDMeanL<-exp(InterceptL)
			APDMeanH<-exp(InterceptH)
			AvgAPD<-c(APDMean,APDMeanL,APDMeanH)

	wdraws<-brms::extract_draws(looAPDw)					
	wMean<-exp(as.numeric(wdraws$dpars$mu$fe$b))			
	wShape<-as.numeric(wdraws$dpars$shape)				
	wScale<-wMean/base::gamma(1+1/wShape)				
	wVar<-(wScale*( base::gamma(1+2/wShape)-(base::gamma(1+1/wShape))^2))/wMean	
	wCVpd<-sqrt(wVar)/wMean										
			CVpd<-wCVpd
	CV<-round(as.numeric(stats::quantile(wCVpd, c(.50, .025, .975))),2)		
	Var<-round(as.numeric(stats::quantile(wVar, c(.50, .025, .975))),2) 
		SD<-round(as.numeric(stats::quantile(sqrt(wVar), c(.50, .025, .975))),2)
	ESS<-wESS

	if (length(RE2)>0) { out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdw),as.numeric(RE2sdw)),3)
				colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd","RE2sd") } else {
		out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdw)),3)
				colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd") }
		message("");message("");message("Output:");message("");message("family: weibull")
		print(out,digits=3)
	} } } 
}
if (frechetGLMM==TRUE) {
	if ( max(as.numeric(rhat(APDf)))<1.05) { 
		if ( fESS>minESS ) {
 
	  if (loolow==loof) { best<-looAPDf; dist="frechet"
	InterceptMean<-fixef(looAPDf)[1]
		InterceptL<-fixef(looAPDf)[3]
		InterceptH<-fixef(looAPDf)[4]
	APDMean<-exp(InterceptMean)
		APDMeanL<-exp(InterceptL)
		APDMeanH<-exp(InterceptH)
		AvgAPD<-c(APDMean,APDMeanL,APDMeanH)
	fdraws<-brms::extract_draws(looAPDf)
		fMean<-exp(as.numeric(fdraws$dpars$mu$fe$b))
		fNu<-as.numeric(fdraws$dpars$nu)
		fVar<-base::gamma(2/fNu+1)-base::gamma(1/fNu+1)*base::gamma(1/fNu+1)
		fCVpd<-sqrt(fVar)/fMean
			CVpd<-fCVpd
		CV<-round(as.numeric(stats::quantile(fCVpd, c(.50, .025, .975))),2)
		Var<-round(as.numeric(stats::quantile(fVar, c(.50, .025, .975))),2)
		SD<-round(as.numeric(stats::quantile(sqrt(fVar), c(.50, .025, .975))),2)
	ESS<-fESS

	if (length(RE2)>0) { out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdf),as.numeric(RE2sdf)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","SD","RE1sd","RE2sd") } else {
		out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdf)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd") }
	message("");message("");message("Output:");message("");message("family: frechet")
	print(out,digits=3)
	} } } 
}
if (gammaGLMM==TRUE) {
	if (max(as.numeric(rhat(APDg)))<1.05) { 
		if (gESS>minESS) {

   if (loolow==loog) { best<-looAPDg; dist="gamma"
	InterceptMean<-fixef(looAPDg)[1]
		InterceptL<-fixef(looAPDg)[3]
		InterceptH<-fixef(looAPDg)[4]
	APDMean<-exp(InterceptMean)
		APDMeanL<-exp(InterceptL)
		APDMeanH<-exp(InterceptH)
		AvgAPD<-c(APDMean,APDMeanL,APDMeanH)
	gdraws<-brms::extract_draws(looAPDg)
		gShape<-as.numeric(gdraws$dpars$shape)
		gMean<-exp(as.numeric(gdraws$dpars$mu$fe$b))
		gSD<-sqrt(gMean/gShape)*gShape	
		gCVpd<-gSD/gMean
			CVpd<-gCVpd
		CV<-round(as.numeric(stats::quantile(gCVpd, c(.50, .025, .975))),2)
		gVar<-(gMean/gShape)*gShape^2		 		
		Var<-round(as.numeric(stats::quantile(gVar, c(.50, .025, .975))),2)
		SD<-round(as.numeric(stats::quantile(gSD, c(.50, .025, .975))),2)
	ESS<-gESS

	if (length(RE2)>0) { out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdg),as.numeric(RE2sdg)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd","RE2sd") } else {
		out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdg)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd") }
	message("");message("");message("Output:");message("");message("family: Gamma")
	print(out,digits=3)
	} } } 
}
if (lognormalGLMM==TRUE) {
	if ( max(as.numeric(rhat(APDln)))<1.05) { 
		if ( lnESS>minESS ) {

   if (loolow==looln) { best<-looAPDln; dist="lognormal"
	InterceptMean<-fixef(looAPDln)[1]
		InterceptL<-fixef(looAPDln)[3]
		InterceptH<-fixef(looAPDln)[4]
	APDMean<-exp(InterceptMean)
		APDMeanL<-exp(InterceptL)
		APDMeanH<-exp(InterceptH)
		AvgAPD<-c(APDMean,APDMeanL,APDMeanH)

	lndraws<-brms::extract_draws(looAPDln)						
		lnMean<-exp(as.numeric(lndraws$dpars$mu$fe$b))			
		lnSigma<-as.numeric(lndraws$dpars$sigma)				
		lnVar<-( exp(lnSigma^2)-1 )*exp(2*lnMean+lnSigma^2)		
		lnCVpd<-sqrt(lnVar)/lnMean						
			CVpd<-lnCVpd
		CV<-round(as.numeric(stats::quantile(lnCVpd, c(.50, .025, .975))),2)
		Var<-round(as.numeric(stats::quantile(lnVar, c(.50, .025, .975))),2)
		SD<-round(as.numeric(stats::quantile(sqrt(lnVar), c(.50, .025, .975))),2)
	ESS<-lnESS

	if (length(RE2)>0) { out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdln),as.numeric(RE2sdln)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd","RE2sd") } else {
		out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdln)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd") }
	message("");message("");message("Output");message("");message("family: lognormal")
	print(out,digits=3)
	} } } 
}
if (invgaussianGLMM==TRUE) {
	if ( max(as.numeric(rhat(APDig)))<1.05) { 
		if ( igESS>minESS ) {
   if (loolow==looig) { best<-looAPDig; dist="invgaussian"
	InterceptMean<-fixef(looAPDig)[1]
		InterceptL<-fixef(looAPDig)[3]
		InterceptH<-fixef(looAPDig)[4]
	APDMean<-exp(InterceptMean)
		APDMeanL<-exp(InterceptL)
		APDMeanH<-exp(InterceptH)
		AvgAPD<-c(APDMean,APDMeanL,APDMeanH)
	igdraws<-brms::extract_draws(looAPDig)
		igMean<-as.numeric(exp(igdraws$dpars$mu$fe$b))
		igShape<-igdraws$dpars$shape
			igVar<- (igMean^3)/igShape
			igCV<- sqrt(igVar)/igMean
			CVpd<-igCV
		CV<-round(as.numeric(stats::quantile(igCV, c(.50, .025, .975))),2)
		Var<-round(as.numeric(stats::quantile(igVar, c(.50, .025, .975))),2)
		SD<-round(as.numeric(stats::quantile(sqrt(igVar), c(.50, .025, .975))),2)
	ESS<-igESS

	if (length(RE2)>0) { out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdig),as.numeric(RE2sdig)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd","RE2sd") } else {
		out<-round(rbind(AvgAPD,CV,SD,as.numeric(RE1sdig)),3); colnames(out)<-c("Mean","HDIlow","HDIhigh"); rownames(out)<-c("AvgAPD","CV","pSD","RE1sd") }
	message("");message("");message("Output");message("");message("family: inverse gaussian")
	print(out,digits=3)
	} } } 
}

	rhatMax<-max(as.numeric(brms::rhat(best)))		
	message(""); message(paste("All MCMC chains reached convergence (rhat<1.05):",rhatMax<1.05,"   Max rhat:",round(rhatMax,4)))
	message(paste("All MCMC chains achieved the minimum effective sample size:",ESS>minESS,"   Minimum ESS:",round(ESS,0)))

  ### Outputs:
	APDraw<-APD(focal)					
	names(loolist)<-loonames[2:length(loolist)]
	loolist<-loolist[2:length(loolist)]
		weibullM<-NULL; frechetM<-NULL; gammaM<-NULL; lognormalM<-NULL; invgaussianM<-NULL
			if (weibullGLMM==TRUE) { if ( max(as.numeric(rhat(APDw)))<1.05) { if (wESS>minESS) { weibullM<-looAPDw }}}
		if (frechetGLMM==TRUE) { if ( max(as.numeric(rhat(APDf)))<1.05) { if (fESS>minESS) { frechetM<-looAPDf }}} 
		if (gammaGLMM==TRUE) { 	if (max(as.numeric(rhat(APDg)))<1.05) { if (wESS>minESS) { gammaM<-looAPDg }}}
		if (lognormalGLMM==TRUE) { if ( max(as.numeric(rhat(APDln)))<1.05) { if (lnESS>minESS) { lognormalM<-looAPDln }}}
		if (invgaussianGLMM==TRUE) { if ( max(as.numeric(rhat(APDig)))<1.05) { if (igESS>minESS) { invgaussianM<-looAPDig }}} 
			allmodels<-list(weibull=weibullM, frechet=frechetM, gamma=gammaM, lognormal=lognormalM, invgaussian=invgaussianM, loolist=loolist)
	raws<-list(rawmean=round(mean(APDraw),2), min=round(min(APDraw),2), max=round(max(APDraw),2),
		range=round(Range,2),IQR=round(IQR,2),qcv=round(qcv,2))
			names(raws)<-c("RawMean","Min","Max","Range","IQR","qcv")
	data<-list(focal=focal, contingent=contingent, RE1=RE1, Re2=RE2)

	raw<-c(mean(APDraw),min(APDraw),max(APDraw),Range,IQR,qcv)
			raw<-round(raw,2)
			names(raw)<-c("Raw Mean","Min","Max","Range","IQR","qcv")

		APDREoutput<-list(data=data, output=out, distribution=dist, 
				model=best, CVposterior=CVpd,
				allmodels=allmodels, rawvalues=APDraw, rawsummary=raws)  
		class(APDREoutput)<-"APD"

	APDout<<-new.env()

	minAPD<-round(min(APD(focal)),2)
		message("");print(raw); message("");
		message(paste("Focal active during contingent inactive time:",minAPD<.01))
		message(paste("Focal active during contingent peak time:",(max(APD(contingent))-max(APD(focal))<0.001)))
		message("");message("")
	if (mean==TRUE) { graphics::abline(h=APDMean,lty=5,lwd=3,col=col)	} 
	if (HDI==TRUE) { graphics::abline(h=c(APDMeanL,APDMeanH),lty=5,lwd=2,col=col)	} 

## Plot:
	..count..<-NULL
	draws<-brms::extract_draws(best)	
	x<-exp(as.numeric(draws$dpars$mu$fe$b))
	hdi<-as.numeric(stats::quantile(x, c(.025,.975)))
	dfAPD<-data.frame(x)
	APDavg<-ggplot2::ggplot(dfAPD, ggplot2::aes(x=x))+ 
		ggplot2::geom_histogram(bins=200,fill=histcol,colour=histcol, ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
		ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="APDavg Estimate",y="Density") +
		ggplot2::geom_segment(ggplot2::aes(x=hdi[1],y=0,xend=hdi[2],yend=0),colour=linecol,size=3.5,lineend="round")

		APDout$hist<<-APDavg

	x<-CVpd
	dfCV<-data.frame(x)
	histCV<-ggplot2::ggplot(dfCV, ggplot2::aes(x=x))+ 
		ggplot2::geom_histogram(bins=200,fill=histcol,colour=histcol, ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
		ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="APD CV Estimate",y="Density") +
		ggplot2::geom_segment(ggplot2::aes(x=CV[2],y=0,xend=CV[3],yend=0),colour=linecol,size=3.5,lineend="round")
	if(length(xlimCV)>0) { histCV<- histCV + xlim(xlimCV) }

		APDout$histCV<<-histCV

	message(""); message("----------------------------------------------------")
			 message("Analysis of Activity Probability Density is complete"); message("")
	
	return(APDREoutput)

}


#' Executable Example of APDRE Function
#' 
#' @description Example of APDRE function using data included in the package
#' @return Provides message with example code for using the APDRE function with data included in the package
#' @export
	exampleAPDRE<-function() {
		message("Try running an example using the data included in the package. Run the following commands:");message("")
		message("WolfBoarAPD<-APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians,")
		message("    RE1=wolfexample$SamplingPeriod, weibullGLMM=TRUE, frechetGLMM=TRUE,")
		message("    invgaussianGLMM=FALSE, gammaGLMM=FALSE)")   }




#' Calculate Raw APD Values
#' 
#' @description Calculates raw APD values, uncorrected for hierarchical data structure
#' @param focal Vector of observations in radians of one species/group/individual/etc. for which predictions on another will be made
#' @param contingent Vector of observations in radians of a species/group/individual/etc. from which predictions will be made
#' @param adjust Smoothing of predicted line; recommended to use default value for observed values and higher value for estimations from circular models; default=1 
#' @return Numeric vector of raw APD values, without correction for nested data structure
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' APDraw(focal=wolfexample$Radians, contingent=boarexample$Radians)
#' 
#' @export
      APDraw<-function(focal, contingent, adjust=1) {
	   xx <- seq(0, 2 * pi, length = 128)
	   bw<-overlap::getBandWidth(contingent, kmax=3)/adjust
	   dens<-overlap::densityFit(contingent, xx, bw)
	   toPlot<-cbind(x=xx, y=dens)
		APD<-stats::approxfun(toPlot); APDvals<-APD(focal)
	}




#' Plot APDavg Posterior Samples
#' @description Plot histogram of samples from posterior distribution for estimated APDavg from \code{\link{APDRE}} function
#' @param model Output from \code{\link{APDRE}} function, an object of class \code{APD}
#' @param histcol Colour for histogram
#' @param linecol Colour for 95% HDI line
#' @return Histogram plot of samples from posterior distribution for estimated APDavg
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' \donttest{ WolfBoarAPD<-APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians,
#'     RE1=wolfexample$SamplingPeriod)
#'     plotAPDavg(WolfBoarAPD) }
#' 
#' @export
	plotAPDavg<-function(model, histcol="deeppink", linecol="black") {
		..count..<-NULL
		draws<-brms::extract_draws(model$model)	
		x<-exp(as.numeric(draws$dpars$mu$fe$b))
		hdi<-as.numeric(stats::quantile(x, c(.025,.975)))
		dfAPD<-data.frame(x)
		ggplot2::ggplot(dfAPD, ggplot2::aes(x=x))+ 
			ggplot2::geom_histogram(bins=200,fill=histcol,colour=histcol, ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="APDavg Estimate",y="Density") +
			ggplot2::geom_segment(ggplot2::aes(x=hdi[1],y=0,xend=hdi[2],yend=0),colour=linecol,size=3.5,lineend="round")
		}




#' Plot APDcv Posterior Samples
#' @description Plot histogram of samples from posterior distribution for estimated APD family-specific population coefficient of variation (CV) from \code{\link{APDRE}} function
#' @param model Output from \code{APDRE} function, an object of class \code{APD}
#' @param histcol Colour for histogram
#' @param linecol Colour for 95% HDI line
#' @return Histogram plot of samples from posterior distribution for estimated APDavg
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' \donttest{ WolfBoarAPD<-APDRE(focal=wolfexample$Radians, contingent=boarexample$Radians,
#'     RE1=wolfexample$SamplingPeriod)
#'     plotAPDcv(WolfBoarAPD) }
#' 
#' @export
	plotAPDcv<-function(model, histcol="deeppink", linecol="black") {
		..count..<-NULL	
		x<-model$CVposterior
		hdi<-as.numeric(stats::quantile(x, c(.025,.975)))
		dfCV<-data.frame(x)
		ggplot2::ggplot(dfCV, ggplot2::aes(x=x))+ 
			ggplot2::geom_histogram(bins=200,fill=histcol,colour=histcol, ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="APD CV Estimate",y="Density") +
			ggplot2::geom_segment(ggplot2::aes(x=hdi[1],y=0,xend=hdi[2],yend=0),colour=linecol,size=3.5,lineend="round")
		}
