#' APDREcorr Correlations between Activity Patterns using APD GLMMs
#' 
#' @description Calculates estimated relationships between activity probability density (APD) of the focal and contingent(s) using Bayesian GLMMs with 'STAN'
#'     using \code{\link[brms]{brm}}, with the option to automatically select the statistical distribution that best fits the dataset (weibull, 
#'     frechet, gamma, lognormal, inverse gaussian) by \code{\link[loo]{loo}}. The function automatically ensures that MCMC chains reach 
#'     convergence and that the specified minimum effective sample size from the posterior distribution is achieved. 
#'     
#' Package: AnimalAPD
#' Version: 1.0.0
#' Date: 2020-11-10
#' @author Liz AD Campbell
#' 
#' @import brms
#' @import circular
#' @import overlap
#' @import ggplot2		
#' @import gridExtra		
#' @import activityGCMM
#' 
#' @keywords cameratrap; activity; temporal; Bayesian
#' 
#' @param focal Vector of observations in radians of one species/group/individual/etc. for which predictions on another will be made.
#' @param cont1 Vector of observations in radians, or output from generalized circular mixture model of activity curves from \code{\link[activityGCMM]{GCMM}}, of a species/group/individual/etc. from which predictions will be made
#' @param cont2 Optional vector of observations in radians, or output from generalized circular mixture model of activity curves from \code{\link[activityGCMM]{GCMM}}, of additional species/group/individual/etc. from which predictions will be made
#' @param cont3 Optional vector of observations in radians, or output from generalized circular mixture model of activity curves from \code{\link[activityGCMM]{GCMM}}, of additional species/group/individual/etc. from which predictions will be made
#' @param cont4 Optional vector of observations in radians, or output from generalized circular mixture model of activity curves from \code{\link[activityGCMM]{GCMM}}, of additional species/group/individual/etc. from which predictions will be made
#' @param RE1 Vector identifying a random intercept for observations of the focal to control for hierarchical data (e.g. camera trap IDs)
#' @param RE2 Optional vector identifying levels of a second random effect, for data with additional hierarchical levels (e.g. study sites, sampling periods, data collection seasons); default is NULL
#' @param weibullGLMM Specifies whether to run a weibull GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param frechetGLMM Specifies whether to run a frechet GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param gammaGLMM Specifies whether to run a Gamma GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param lognormalGLMM Specifies whether to run a lognormal GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param invgaussianGLMM Specifies whether to run a inverse.gaussian GLMM, using the brms package; default is TRUE for all and results from the best-fitting model are returned
#' @param cores Number of cores to use when running MCMC chains in parallel; default=1
#' @param iter Number of MCMC iteractions per chain; burnin is iter/2; default=5000
#' @param thin Thinning rate for saving MCMC draws; default=1
#' @param center Value to use as center of graph; default=pi
#' @param adjust Smoothing of predicted line; recommended to use default value for observed values and higher value for estimations from circular models 
#' @param Reloo Whether to use reloo when running leave-one-out cross-validation of models (loo); see also \code{\link{brms}} and \code{\link{loo}} 
#' @param adapt_delta Value to use for adapt_delta with brms; default=0.95; see also \code{\link{brms}}
#' @param burnin Number of MCMC iterations to be discarded as the burn-in; default=iter/2
#' @param minESS Desired minimum effective sample size; default=1000
#' @param plothist Whether to plot histograms of samples from the posterior distribution for the correlation parameters; default=TRUE		
#' @param ploteffects Whether to plot predicted effects; default=TRUE
#' @param histcol Colour for histogram bars							
#' @param linecol Colour for histogram lines for the 95% HDI and 0			
#' @param effectcol Colour for predicted effect plot 95% HDI				
#' 
#' @return Prints results of best-fitting model and posterior samples and/or predicted effects of parameter estimates if \code{plothist=TRUE} 
#'     and \code{ploteffects=TRUE}, and returns object of class \code{APD} with list of analysis results and information. 
#' @return \code{data} List of data used in analysis
#' @return \code{model} Object of class \code{brmsfit} containing results and information for best-fitting model.
#' @return \code{distribution} Character vector of statistical distribution of best-fitting model
#' @return \code{allmodels} List of output for all tested models; object of class \code{brmsfit}
#' 
#' @seealso \code{\link[activityGCMM]{GCMM}} \code{\link[brms]{brm}} \code{\link[loo]{loo}}
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' \donttest{APDREcorr(focal=wolfexample$Radians,cont1=boarexample$Radians,
#'     RE1=wolfexample$SamplingPeriod)}
#' 
#' 
#' @export
	APDREcorr<- function(focal, cont1, cont2=NULL, cont3=NULL, cont4=NULL, RE1, RE2=NULL, 
		weibullGLMM=TRUE,frechetGLMM=TRUE,gammaGLMM=TRUE,lognormalGLMM=FALSE,invgaussianGLMM=TRUE,
		cores=1, iter=5000, minESS=1000, burnin=iter/2, thin=1, adapt_delta=0.95, center="pi", adjust=1, Reloo=TRUE,
		plothist=TRUE, ploteffects=TRUE, histcol="cyan4", effectcol="cyan4", linecol="red") {

	requireNamespace("overlap",quietly=TRUE)
	requireNamespace("brms",quietly=FALSE)
	requireNamespace("circular",quietly=TRUE)
	requireNamespace("stats",quietly=TRUE)
	requireNamespace("loo",quietly=TRUE)
	requireNamespace("ggplot2",quietly=TRUE)			
	requireNamespace("gridExtra",quietly=TRUE)		
		requireNamespace("activityGCMM",quietly=TRUE)

	APDcOut<-NULL


	focal<-as.numeric(focal)
	RE1<-as.numeric(as.factor(RE1))
	if (length(RE2)>0) { RE2<-as.numeric(as.factor(RE2)) }
	adaptdelta<-adapt_delta

	xx <- seq(0, 2 * pi, length = 128)
	# focal:
	bw<-overlap::getBandWidth(focal, kmax=3)/adjust
	   dens<-overlap::densityFit(focal, xx, bw)
	   toPlot<-cbind(x=xx, y=dens)
		APDfocal<-stats::approxfun(toPlot)
		APDfocal<-APDfocal(focal)
	# C1:
	if(class(cont1)=="GCMM"){ cont1<-cont1$GCMMmixture; adjust<-5 }
	bw<-overlap::getBandWidth(cont1, kmax=3)/adjust
	   dens<-overlap::densityFit(cont1, xx, bw)
	   toPlot<-cbind(x=xx, y=dens)
		APD1<-stats::approxfun(toPlot); APD1<-APD1(focal)
	df<-data.frame(APDfocal,APD1)
	# C2:
	if(length(cont2)>0) {
		if(class(cont2)=="GCMM") { cont2<-cont2$GCMMmixture; adjust<-5 }
	bw<-overlap::getBandWidth(cont2, kmax=3)/adjust
	   dens<-overlap::densityFit(cont1, xx, bw)
	   toPlot<-cbind(x=xx, y=dens)
		APD2<-stats::approxfun(toPlot); APD2<-APD2(focal)
		df<-data.frame(df, APD2)	}
	# C3:
	if(length(cont3)>0) { 
		if(class(cont2)=="GCMM") { cont3<-cont3$GCMMmixture; adjust=5 }
		bw<-overlap::getBandWidth(cont3, kmax=3)/adjust
	   dens<-overlap::densityFit(cont3, xx, bw); toPlot<-cbind(x=xx, y=dens)
		APD3<-stats::approxfun(toPlot); APD3<-APD3(focal); df<-data.frame(df, APD3)	}
	# C4:
	if(length(cont4)>0) { 
		if(class(cont4)=="GCMM") { cont4<-cont4$GCMMmixture; adjust=5 }
		bw<-overlap::getBandWidth(cont4, kmax=3)/adjust
	   dens<-overlap::densityFit(cont4, xx, bw); toPlot<-cbind(x=xx, y=dens)
		APD4<-stats::approxfun(toPlot);APD4<-APD4(focal); df<-data.frame(df, APD4)	}

	if (length(RE2)>0) { df<-data.frame(df,RE1,RE2) } else { df<-data.frame(df,RE1) }


loolist<-c(-100000); loonames<-c("X")
## Weibull:
if (weibullGLMM==TRUE) {
	message("---------------------")
	message("Running Weibull model")
	message("  ")
	if (length(RE2)>0) { 
           if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,warmup=burnin,thin=thin,family=weibull,cores=cores,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,warmup=burnin,thin=thin,family=weibull,cores=cores,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,warmup=burnin,thin=thin,family=weibull,cores=cores,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDw<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,warmup=burnin,thin=thin,family=weibull,cores=cores,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
	if (!length(RE2)>0) { 
           if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDw<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
		i<-1;	if ( max(as.numeric(brms::rhat(APDw))) > 1.05) { repeat {	i<-i+1 
			if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
            		                       APDw<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           				 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,warmup=burnin,family=weibull,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
                                   APDw<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=weibull,warmup=burnin,cores=cores,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           				 } 
				if ( max(as.numeric(brms::rhat(APDw))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for weibull model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
			} }
		if ( max(as.numeric(rhat(APDw))) < 1.05) {
			message(""); message("All MCMC chains reached convergence (rhat<1.05)"); message("")
				i2<-1	
			wESS<-min(as.numeric(neff_ratio(APDw)))*(iter*(i+i2)/2*3)
			if ( wESS<minESS ) { repeat { i2<-i2+1
			  if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=weibull,warmup=burnin,cores=cores,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=weibull,warmup=burnin,cores=cores,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=weibull,warmup=burnin,cores=cores,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
            		                       APDw<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
			  if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=weibull,warmup=burnin,cores=cores,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDw<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=weibull,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
                                   APDw<-brms::brm(APDfocal~APD1+(1|RE1),data=df,family=weibull,thin=thin,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
				    wESS<-min(as.numeric(neff_ratio(APDw)))*(iter*(i+i2)/2*3)
					if ( wESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for weibull model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(wESS,0))); break }
				} }
		if ( wESS>minESS ) { 							 
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(wESS,0))); message("") } 
		}
	if ( max(as.numeric(rhat(APDw))) < 1.05) { if ( wESS>minESS ) {
	looAPDw<-brms::add_criterion(APDw, "loo", reloo=Reloo)
		loow<-loo::loo(looAPDw)$estimates[1,1]
		loolist<-c(loolist,loow); loonames<-c(loonames,"weibull")
			message("Weibull model loo:") 
			print(loow); message("") 
		} }
	} 

## Frechet:	
if (frechetGLMM==TRUE) {
	message("---------------------")
	message("Running Frechet model")
	message("  ")
	if (length(RE2)>0) { 
           if (length(cont4)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDf<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
	if (!length(RE2)>0) { 
           if (length(cont4)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDf<-brms::brm(APDfocal~APD1+(1|RE1),data=df,family=frechet,thin=thin,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
		i<-1;	if ( max(as.numeric(brms::rhat(APDf))) > 1.05) { repeat {	i<-i+1 
			if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
            		                       APDf<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
                                               APDf<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
				if ( max(as.numeric(brms::rhat(APDf))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for frechet model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
			} }
		if ( max(as.numeric(rhat(APDf))) < 1.05) {  i2<-1	
				  fESS<-min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)
			if ( fESS<minESS ) { repeat { i2<-i2+1
			if (length(RE2)>0) { 
      		      if (length(cont4)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDf<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
            		                       APDf<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDf<-brms::brm(APDf~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDf<-brms::brm(APDf~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDf<-brms::brm(APDf~APD1+APD2+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
                                               APDf<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=frechet,cores=cores,chains=3,warmup=burnin,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
				  fESS<-min(as.numeric(neff_ratio(APDf)))*(iter*(i+i2)/2*3)
					if ( fESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for frechet model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(fESS,0))); break }
				} }
		if ( fESS>minESS ) { 							 
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(fESS,0))); message("") } 
		}
	if ( max(as.numeric(rhat(APDf))) < 1.05) { if ( fESS>minESS ) {
	  looAPDf<-brms::add_criterion(APDf, "loo", reloo=Reloo)
		loof<-loo::loo(looAPDf)$estimates[1,1]
		loolist<-c(loolist,loof); loonames<-c(loonames,"frechet")
			message("Frechet model loo:") 
			print(loof); message("") 
		} }
	} 



## Gamma:
if (gammaGLMM==TRUE) {
	message("---------------------")
	message("Running Gamma model")
	message("  ")
	if (length(RE2)>0) { 
           if (length(cont4)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDg<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
	if (!length(RE2)>0) { 
           if (length(cont4)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDg<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
		i<-1;	if ( max(as.numeric(brms::rhat(APDg))) > 1.05) { repeat {	i<-i+1 
			if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
            		                       APDg<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
                                               APDg<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
				if ( max(as.numeric(brms::rhat(APDg))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for gamma model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
			} }
		if ( max(as.numeric(rhat(APDg))) < 1.05) {  i2<-1	
			gESS<-min(as.numeric(neff_ratio(APDg)))*(iter*(i+i2)/2*3)
			if ( gESS<minESS ) { repeat { i2<-i2+1
			if (length(RE2)>0) { 
      		      if (length(cont4)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDg<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
            		                       APDg<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDg<-brms::brm(APDf~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDg<-brms::brm(APDf~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDg<-brms::brm(APDf~APD1+APD2+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
                                               APDg<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=stats::Gamma(link="log"),cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
				  gESS<-min(as.numeric(neff_ratio(APDg)))*(iter*(i+i2)/2*3)
					if ( gESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for gamma model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(gESS,0))); break }
				} }
		if ( gESS>minESS ) { 							 
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(gESS,0))); message("") } 
		}
	if ( max(as.numeric(rhat(APDg))) < 1.05) { if ( gESS>minESS ) {
	  looAPDg<-brms::add_criterion(APDg, "loo", reloo=Reloo)
		loog<-loo::loo(looAPDg)$estimates[1,1]
		loolist<-c(loolist,loog); loonames<-c(loonames,"gamma")
			message("Gamma model loo:") 
			print(loog); message("") 
		} }
	} 


#lognormal:
if (lognormalGLMM==TRUE) {
	message("---------------------")
	message("Running lognormal model")
	message("  ")
	if (length(RE2)>0) { 
           if (length(cont4)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDln<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
	if (!length(RE2)>0) { 
           if (length(cont4)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDln<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,chains=3,warmup=burnin,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
		i<-1;	if ( max(as.numeric(brms::rhat(APDln))) > 1.05) { repeat {	i<-i+1 
			if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
            		                       APDln<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
                                               APDln<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
				if ( max(as.numeric(brms::rhat(APDln))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for lognormal model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
			} }
		if ( max(as.numeric(rhat(APDln))) < 1.05) {  i2<-1	
			lnESS<-min(as.numeric(neff_ratio(APDln)))*(iter*(i+i2)/2*3)
			if ( lnESS<minESS ) { repeat { i2<-i2+1
			if (length(RE2)>0) { 
      		      if (length(cont4)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDln<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
            		                       APDln<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDln<-brms::brm(APDf~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDln<-brms::brm(APDf~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDln<-brms::brm(APDf~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
                                               APDln<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
				  lnESS<-min(as.numeric(neff_ratio(APDln)))*(iter*(i+i2)/2*3)
					if ( lnESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for gamma model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(lnESS,0))); break }
				} }
		if ( lnESS>minESS ) { 							 
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(lnESS,0))); message("") }
		}
	if ( max(as.numeric(rhat(APDln))) < 1.05) { if ( lnESS>minESS ) {
	  looAPDln<-brms::add_criterion(APDln, "loo", reloo=Reloo)
		looln<-loo::loo(looAPDln)$estimates[1,1]
		loolist<-c(loolist,looln); loonames<-c(loonames,"lognormal")
			message("Lognormal model loo:") 
			print(looln); message("") 
		} }
	} 


#inverse gaussian: 
if (invgaussianGLMM==TRUE) {
	message("---------------------")
	message("Running inverse gaussian model")
	message("  ")
	if (length(RE2)>0) { 
           if (length(cont4)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDig<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
	if (!length(RE2)>0) { 
           if (length(cont4)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont3)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
		if (length(cont2)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta)) } else
                                   APDig<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter,control=list(adapt_delta=adaptdelta))
            } 
		i<-1;	if ( max(as.numeric(brms::rhat(APDig))) > 1.05) { repeat {	i<-i+1 
			if (length(RE2)>0) { 
      		     if (length(cont4)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
            		                       APDig<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta)) } else
                                               APDig<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*i,control=list(adapt_delta=adaptdelta))
           		 } 
				if ( max(as.numeric(brms::rhat(APDig))) < 1.05) { break }
				if ( i==5 ) { print("Maximum attempts reached for inverse gaussian model. Model stopped before all rhat < 1.05. Try increasing iter"); break }
			} }
		if ( max(as.numeric(rhat(APDig))) < 1.05) {  i2<-1	
			igESS<-min(as.numeric(neff_ratio(APDig)))*(iter*(i+i2)/2*3)
			if ( igESS<minESS ) { repeat { i2<-i2+1
			if (length(RE2)>0) { 
      		      if (length(cont4)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+APD4+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+APD3+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDig<-brms::brm(APDfocal~APD1+APD2+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
            		                       APDig<-brms::brm(APDfocal~APD1+(1|RE1)+(1|RE2),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
			if (!length(RE2)>0) { 
      		      if (length(cont4)>0) { APDig<-brms::brm(APDf~APD1+APD2+APD3+APD4+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont3)>0) { APDig<-brms::brm(APDf~APD1+APD2+APD3+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
				if (length(cont2)>0) { APDig<-brms::brm(APDf~APD1+APD2+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta)) } else
                                               APDig<-brms::brm(APDfocal~APD1+(1|RE1),data=df,thin=thin,family=lognormal,cores=cores,warmup=burnin,chains=3,iter=iter*(i+i2),control=list(adapt_delta=adaptdelta))
           		 } 
				  igESS<-min(as.numeric(neff_ratio(APDig)))*(iter*(i+i2)/2*3)
					if ( igESS>minESS ) { break }
					if ( i2==5 ) { message(paste("Maximum attempts reached for inverse gaussian model. Model stopped before all ESS > minESS. Try increasing iter. Minimum ESS =",round(igESS,0))); break }
				} }
		if ( igESS>minESS ) { 							 
			message(""); message(paste("All MCMC chains achieved the minimum effective sample size. Minimum ESS =",round(igESS,0))); message("") } 
		}
	if ( max(as.numeric(rhat(APDig))) < 1.05) { if ( igESS>minESS ) {
	  looAPDig<-brms::add_criterion(APDig, "loo", reloo=Reloo)
		looig<-loo::loo(looAPDig)$estimates[1,1]
		loolist<-c(loolist,looig); loonames<-c(loonames,"inv gaussian")
			message("Inverse gaussian model loo:") 
			print(looig); message("") 
		} }
	} 

 
	loolow<-max(loolist) 		
	  if (weibullGLMM==TRUE) { if (loolow==loow) { best<-looAPDw; ESS<-wESS; dist="weibull"  }  }
	  if (frechetGLMM==TRUE) { if (loolow==loof) { best<-looAPDf; ESS<-fESS; dist="frechet"  }  }
	  if (gammaGLMM==TRUE) { if (loolow==loog) { best<-looAPDg; ESS<-gESS; dist="gamma"	}  }
	  if (lognormalGLMM==TRUE) { if (loolow==looln) { best<-looAPDln; ESS<-lnESS; dist="lognormal"	}  }
	  if (invgaussianGLMM==TRUE) { if (loolow==looig) { best<-looAPDig; ESS<-igESS; dist="invgaussian" } } 
	APDcorrmodel<-best
	print(APDcorrmodel)
		rhatMax<-max(as.numeric(brms::rhat(best)))	
		message(""); message(paste("All MCMC chains reached convergence (rhat<1.05):",rhatMax<1.05,"   Max rhat:",round(rhatMax,4)))	
		message(paste("All MCMC chains achieved the minimum effective sample size:",ESS>minESS,"   Minimum ESS:",round(ESS,0)))

	### Outputs:
	names(loolist)<-loonames					
		loolist<<-loolist[2:length(loolist)]			
	weibullM<-NULL; frechetM<-NULL; gammaM<-NULL; lognormalM<-NULL; invgaussianM<-NULL
		if (weibullGLMM==TRUE) { if (max(as.numeric(rhat(APDw)))<1.05) { if (wESS>minESS) { weibullM<-looAPDw }}} 
		if (frechetGLMM==TRUE) { if (max(as.numeric(rhat(APDf)))<1.05) { if (fESS>minESS) { frechetM<-looAPDf }}} 
		if (gammaGLMM==TRUE) { if (max(as.numeric(rhat(APDg)))<1.05) { if (gESS>minESS) { gammaM<-looAPDg }}} 
		if (lognormalGLMM==TRUE) { if (max(as.numeric(rhat(APDln)))<1.05) { if (lnESS>minESS) { lognormalM<-looAPDln }}} 
		if (invgaussianGLMM==TRUE) { if (max(as.numeric(rhat(APDig)))<1.05) { if (igESS>minESS) { invgaussianM<-looAPDig }}} 
	allmodels<-list(weibull=weibullM, frechet=frechetM, gamma=gammaM, lognormal=lognormalM, invgaussian=invgaussianM, loolist=loolist)
	data<-list(focal=focal, cont1=cont1, cont2=cont2, cont3=cont3, cont4=cont4, RE1=RE1, Re2=RE2)

	APDcOutput<-list(data=data,  model=best, distribution=dist, allmodels=allmodels) 
		class(APDcOutput)<-"APD"


#### Posterior distribution plots: 
	..count..<-NULL
	if (plothist==TRUE) {
		draws<-brms::extract_draws(APDcorrmodel)
		draws<-draws$dpars$mu$fe$b
		x<-draws[,2]
		hdi1<-as.numeric(stats::quantile(draws[,2], c(.025, .975))) 
		dfdraws<-data.frame(x) 
			h1<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+ 
				ggplot2::geom_histogram(bins=200,fill="cyan4",colour="cyan4", ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
				ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="C1 Parameter Estimate",y="Density") +
				ggplot2::geom_vline(xintercept=0, color="red", size=1.5) +  
				ggplot2::geom_segment(ggplot2::aes(x=hdi1[1],y=0,xend=hdi1[2],yend=0),colour = "red",size=3.5,lineend="round")

	if(length(cont2)>0) {
		x<-draws[,3]
		dfdraws<-data.frame(x) 
		hdi2<-as.numeric(stats::quantile(draws[,3], c(.025, .975)))
		h2<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+
			ggplot2::geom_histogram(bins=200,fill="cyan4",colour="cyan4", ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="C2 Parameter Estimate",y="Density") +
			ggplot2::geom_vline(xintercept=0, color="red", size=1.5) + #linetype, 
			ggplot2::geom_segment(ggplot2::aes(x=hdi2[1],y=0,xend=hdi2[2],yend=0),colour = "red",size=3.5,lineend="round")
		}
	if(length(cont3)>0) {
		x<-draws[,4]
		dfdraws<-data.frame(x)
		hdi3<-as.numeric(stats::quantile(draws[,4], c(.025, .975)))
		h3<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+
			ggplot2::geom_histogram(bins=200,fill="cyan4",colour="cyan4", ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="C3 Parameter Estimate",y="Density") +
			ggplot2::geom_vline(xintercept=0, color="red", size=1.5) +  
			ggplot2::geom_segment(ggplot2::aes(x=hdi3[1],y=0,xend=hdi3[2],yend=0),colour = "red",size=3.5,lineend="round")
		}
	if(length(cont4)>0) {
		x<-draws[,5]
		dfdraws<-data.frame(x)#,hdi)
		hdi4<-as.numeric(stats::quantile(draws[,5], c(.025, .975)))
		h4<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+
			ggplot2::geom_histogram(bins=200,fill="cyan4",colour="cyan4", ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x="C4 Parameter Estimate",y="Density") +
			ggplot2::geom_vline(xintercept=0, color="red", size=1.5) + 
			ggplot2::geom_segment(ggplot2::aes(x=hdi4[1],y=0,xend=hdi4[2],yend=0),colour = "red",size=3.5,lineend="round")
		}
  if(ploteffects==FALSE) {
	if(length(cont4)>0) { gridExtra::grid.arrange(h1,h2,h3,h4, nrow=1) } else
	if(length(cont3)>0) { gridExtra::grid.arrange(h1,h2,h3, nrow=1) } else
	if(length(cont2)>0) { gridExtra::grid.arrange(h1,h2, nrow=1) } else
				    h1
	}
 }

	
#### Effect Plots:
if(ploteffects==TRUE) {
	plt<-brms::marginal_effects(APDcorrmodel)
	x<-plt$APD1$APD1; y<-plt$APD1$estimate
	c1 <- data.frame(x=x, y=y)
		eb1 <- ggplot2::aes(ymax=plt$APD1$upper, ymin=plt$APD1$lower)
		p1<-ggplot2::ggplot(data=c1, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size = 10) +
			ggplot2::geom_ribbon(eb1,fill=effectcol,colour="white",alpha=.5)+
			ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x="C1 APD",y="Focal APD")
if(length(cont2)>0) {
	x<-plt$APD2$APD2; y<-plt$APD2$estimate
	c2 <- data.frame(x=x, y=y)
		eb2 <- ggplot2::aes(ymax=plt$APD2$upper, ymin=plt$APD2$lower)
		p2<-ggplot2::ggplot(data=c2, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size = 10) +
			ggplot2::geom_ribbon(eb2,fill=effectcol,colour="white",alpha=.5)+
			ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x="C2 APD",y="Focal APD")
	}
if(length(cont3)>0) {
	x<-plt$APD3$APD3; y<-plt$APD3$estimate
	c3 <- data.frame(x=x, y=y)
		eb3 <- ggplot2::aes(ymax=plt$APD3$upper, ymin=plt$APD3$lower)
		p3<-ggplot2::ggplot(data=c3, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size = 10) +
			ggplot2::geom_ribbon(eb3,fill=effectcol,colour="white",alpha=.5)+
			ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x="C3 APD",y="Focal APD")
	}
if(length(cont4)>0) {
	x<-plt$APD4$APD4; y<-plt$APD4$estimate
	c4 <- data.frame(x=x, y=y)
		eb4 <- ggplot2::aes(ymax=plt$APD4$upper, ymin=plt$APD4$lower)
		p4<-ggplot2::ggplot(data=c4, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size = 10) +
			ggplot2::geom_ribbon(eb4,fill=effectcol,colour="white",alpha=.5)+
			ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x="C4 APD",y="Focal APD")
	}
  if (plothist==FALSE) {
	if(length(cont4)>0) { gridExtra::grid.arrange(p1,p2,p3,p4, nrow=1) } else
	if(length(cont3)>0) { gridExtra::grid.arrange(p1,p2,p3, nrow=1) } else
	if(length(cont2)>0) { gridExtra::grid.arrange(p1,p2, nrow=1) } else
				    p1
	} else
	if(length(cont4)>0) { gridExtra::grid.arrange(h1,h2,h3,h4,p1,p2,p3,p4, nrow=2) } else
	if(length(cont3)>0) { gridExtra::grid.arrange(h1,h2,h3,p1,p2,p3, nrow=2) } else
	if(length(cont2)>0) { gridExtra::grid.arrange(h1,h2,p1,p2, nrow=2) } else
				    gridExtra::grid.arrange(h1,p1, nrow=2)
}


	message("");message("Analysis of Activity Probability Density correlation is complete")

	return(APDcOutput)

}


#' Executable Example of APDREcorr Function 
#' @description Example of APDREcorr function using data included in the package
#' @return Prints message with example code for using the APDRE function using data included in the package
#' @export
	exampleAPDREcorr<-function() {
		message("Try running an example using the data included in the package. Run the following commands:"); message("")
		message("WolfBoarAPDcorr<-APDREcorr(focal=wolfexample$Radians, cont1=boarexample$Radians,")
		message("    RE1=wolfexample$SamplingPeriod, iter=2000, weibullGLMM=TRUE, frechetGLMM=FALSE,")
		message("    invgaussianGLMM=FALSE, gammaGLMM=TRUE)")   }





#' Plot APDcorr Posterior Samples and Predicted Effects
#'
#' @description Plot histogram of posterior distribution samples and/or predicted effects of parameter estimated from \code{\link{APDREcorr}} function
#' @param model Object of class \code{APD} of model output from \code{\link{APDREcorr}} function
#' @param hist Logical argument for whether to plot histogram of posterior distribution samples; default=TRUE
#' @param effects Logical argument for whether to plot predicted effects; default=TRUE
#' @param yname Character vector of name of focal, to be used in the y axis labels
#' @param xname Character vector of name(s) of contingent(s), to be used in the x axis labels
#' @param histcol Colour for posterior distribution histogram
#' @param linecol Colour for 95% HDI line on posterior distribution histogram
#' @param effectcol Colour for 95% HDI ribbon in line plot of predicted effects
#' 
#' @return Histogram plot(s) of samples from posterior distribution for estimated relationship between focal and contingent(s), if plothist=TRUE
#' @return Plot of predicted effects (mean and 95% highest density interval) for relationship between focal APD and contingent(s) APD, if ploteffects=TRUE
#' 
#' @examples
#' data(wolfexample)
#' data(boarexample)
#' \donttest{ WolfBoarAPDc<-APDREcorr(focal=wolfexample$Radians, cont1=boarexample$Radians,
#'     RE1=wolfexample$SamplingPeriod)
#'     plotAPDcorr(WolfBoarAPDc, yname="Wolf", xname="Boar") }
#' 
#' @export
	plotAPDcorr <- function (model, hist=TRUE, effects=TRUE, histcol=c("cyan4","cyan4","cyan4","cyan4"),
		linecol="red", effectcol=c("cyan4","cyan4","cyan4","cyan4"), yname="Focal", xname=c("C1","C2","C3","C4") ) {

	..count..<-NULL

## Posterior distribution plots:
if (hist==TRUE) {
	draws<-brms::extract_draws(model$model)
	draws<-draws$dpars$mu$fe$b

	x<-draws[,2]
	hdi1<-as.numeric(stats::quantile(draws[,2], c(.025, .975)))
	dfdraws<-data.frame(x)
		h1<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+ 
			ggplot2::geom_histogram(bins=200,fill=histcol[1],colour=histcol[1], ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x=(paste(xname[1],"Parameter Estimate")),y="Density") +
			ggplot2::geom_vline(xintercept=0, color=linecol, size=1.5) +  
			ggplot2::geom_segment(ggplot2::aes(x=hdi1[1],y=0,xend=hdi1[2],yend=0),colour=linecol,size=3.5,lineend="round")

	if(ncol(draws)>2) {
		x<-draws[,3]
		dfdraws<-data.frame(x)
		hdi2<-as.numeric(stats::quantile(draws[,3], c(.025, .975)))
		h2<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x))+
			ggplot2::geom_histogram(bins=200,fill=histcol[2],colour=histcol[2], ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x=(paste(xname[2],"Parameter Estimate")),y="Density") +
			ggplot2::geom_vline(xintercept=0, color=linecol, size=1.5) +  
			ggplot2::geom_segment(ggplot2::aes(x=hdi2[1],y=0,xend=hdi2[2],yend=0),colour=linecol,size=3.5,lineend="round")
		}
	if(ncol(draws)>3) {
		x<-draws[,4]
		dfdraws<-data.frame(x)
		hdi3<-as.numeric(stats::quantile(draws[,4], c(.025, .975)))
		h3<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x)) + 
			ggplot2::geom_histogram(bins=200,fill=histcol[3],colour=histcol[3], ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size = 10)+ ggplot2::labs(x=(paste(xname[3],"Parameter Estimate")),y="Density") +
			ggplot2::geom_vline(xintercept=0, color=linecol, size=1.5) + 
			ggplot2::geom_segment(ggplot2::aes(x=hdi3[1],y=0,xend=hdi3[2],yend=0),colour=linecol,size=3.5,lineend="round")
		}
	if(ncol(draws)>4) {
		x<-draws[,5]
		dfdraws<-data.frame(x)
		hdi4<-as.numeric(stats::quantile(draws[,5], c(.025, .975)))
		h4<- ggplot2::ggplot(dfdraws, ggplot2::aes(x=x)) + 
			ggplot2::geom_histogram(bins=200,fill=histcol[4],colour=histcol[4], ggplot2::aes(y = ggplot2::stat(..count..) / sum(..count..))) +
			ggplot2::theme_classic(base_size=10)+ ggplot2::labs(x=(paste(xname[4],"Parameter Estimate")),y="Density") +
			ggplot2::geom_vline(xintercept=0, color=linecol, size=1.5) +  
			ggplot2::geom_segment(ggplot2::aes(x=hdi4[1],y=0,xend=hdi4[2],yend=0),colour=linecol,size=3.5,lineend="round")
		}
  if(effects==FALSE) {
	if(ncol(draws)>4) { gridExtra::grid.arrange(h1,h2,h3,h4, nrow=1) } else
	if(ncol(draws)>3) { gridExtra::grid.arrange(h1,h2,h3, nrow=1) } else
	if(ncol(draws)>2) { gridExtra::grid.arrange(h1,h2, nrow=1) } else
				    h1
	}
   }

	
## Effect Plots:
  if(effects==TRUE) {
		plt<-brms::marginal_effects(model$model)

		x<-plt$APD1$APD1; y<-plt$APD1$estimate
		c1 <- data.frame(x=x, y=y)
		eb1 <- ggplot2::aes(ymax=plt$APD1$upper, ymin=plt$APD1$lower)
		p1<-ggplot2::ggplot(data=c1, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size=10) +
			ggplot2::geom_ribbon(eb1,fill=effectcol[1],colour="white",alpha=.5)+
			ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x=(paste(xname[1],"APD")),y=(paste(yname,"APD")))

		if(length(plt$APD2)>0) {
			x<-plt$APD2$APD2; y<-plt$APD2$estimate
			c2 <- data.frame(x=x, y=y)
			eb2 <- ggplot2::aes(ymax=plt$APD2$upper, ymin=plt$APD2$lower)
			p2<-ggplot2::ggplot(data=c2, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size=10) +
				ggplot2::geom_ribbon(eb2,fill=effectcol[2],colour="white",alpha=.5)+
				ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x=(paste(xname[1],"APD")),y=(paste(yname,"APD")))
			}

		if(length(plt$APD3)>0) {
			x<-plt$APD3$APD3; y<-plt$APD3$estimate
			c3 <- data.frame(x=x, y=y)
			eb3 <- ggplot2::aes(ymax=plt$APD3$upper, ymin=plt$APD3$lower)
			p3<-ggplot2::ggplot(data=c3, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size=10) +
				ggplot2::geom_ribbon(eb3,fill=effectcol[3],colour="white",alpha=.5)+
				ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x=(paste(xname[1],"APD")),y=(paste(yname,"APD")))
			}

		if(length(plt$APD4)>0) {
			x<-plt$APD4$APD4; y<-plt$APD4$estimate
			c4 <- data.frame(x=x, y=y)
			eb4 <- ggplot2::aes(ymax=plt$APD4$upper, ymin=plt$APD4$lower)
			p4<-ggplot2::ggplot(data=c4, ggplot2::aes(x=x, y=y)) + ggplot2::theme_classic(base_size=10) +
				ggplot2::geom_ribbon(eb4,fill=effectcol[4],colour="white",alpha=.5)+
				ggplot2::geom_line(colour="black",size=1.2)+ ggplot2::labs(x=(paste(xname[1],"APD")),y=(paste(yname,"APD")))
			}

	if (hist==FALSE) {
		if(length(plt$APD4)>0){ gridExtra::grid.arrange(p1,p2,p3,p4, nrow=1) } else
		if(length(plt$APD3)>0){ gridExtra::grid.arrange(p1,p2,p3, nrow=1) } else
		if(length(plt$APD2)>0){ gridExtra::grid.arrange(p1,p2, nrow=1) } else
						p1
			} else
		if(ncol(draws)>4){ gridExtra::grid.arrange(h1,h2,h3,h4,p1,p2,p3,p4, nrow=2) } else
		if(ncol(draws)>3){ gridExtra::grid.arrange(h1,h2,h3,p1,p2,p3, nrow=2) } else
		if(ncol(draws)>2){ gridExtra::grid.arrange(h1,h2,p1,p2, nrow=2) } else
					 gridExtra::grid.arrange(h1,p1, nrow=2)
		}
}


