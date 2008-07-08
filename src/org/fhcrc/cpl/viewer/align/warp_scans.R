##Adapted from align_pair.R.  This does less than that did.  Given pairs of features that
##are thought to be associated, create a warping function that maps scans in one featureset
##to scans in another.

map_file_scans_onto_file=function(featurePairsFilename,outFilename,maxScan=10000,bigIntensity=100,topN=NA,df=20) {
    ## Create a mapping of scans, based on feature pairs in featurePairsFilename,
    ## writing the mapping out to outFilename
    feature_pairs=read.table(featurePairsFilename, header=T, fill=T, sep="\t")
	warp=align_pair(feature_pairs,df)
	source_scan = (c(1:maxScan))
	dest_scan = warp(source_scan)
#	dest_scan = round(dest_scan)
	scan_map = cbind(source_scan, dest_scan)
	write.table(scan_map,outFilename,row.names=F,quote=F,na="",sep="\t")
    }

readf=function(fn) {
	f=read.table(fn, header=T, fill=T, sep="\t")
	}

linear_align=function(s1,s2,sigma) {
	### Make a coarse linear alignment
	ds=s2-s1
	med.l=median((ds)[s1<median(s1)])
	med.r=median((ds)[s1>median(s1)])
	q1=quantile(s1,.25)
	q3=quantile(s1,.75)
	s1t=(s1-q1)*(med.r-med.l)/(q3-q1)+med.l+s1

	### Repeatedly try to refine the linear alignment by shifting the
        ### constant and then the slope to minimize the sum of min( (residuals/sigma)^2 ,1)
	for (reps in 1:3) {
	cv=seq(-1000,1000,length=200)
	cscore=0*cv
	for (i in 1:length(cv)) {
		ds=s2-(s1t+cv[i])
		cscore[i]=sum(pmin((ds/sigma)^2,1))
		}
	c=cv[which.min(cscore)]
	s1t=s1t+c
	bv=exp(seq(-2,2,length=200))
	bscore=0*bv
	m1=median(s1t)
	for (i in 1:length(bv)) {
		ds=s2-(bv[i]*(s1t-m1)+m1)
		bscore[i]=sum(pmin((ds/sigma)^2,1))
		}
	b=bv[which.min(bscore)]
	s1t=b*(s1t-m1)+m1
	}
	##plot(s1t,s2-s1t,ylim=c(-300,300))

	##### Now s1t is a suitable linear transform of s1; recover the coefficients

	lin.coef=lm(s1t~s1)$coef
	}

## Robust smoother that uses weighted smooth splines to simulate
## a smoother with either a Huber-type or a bounded loss function
## control the smoothness with the df (degrees of freedom) parameter
## for tau<=0, tau is ignored and the loss is of Huber type:
##       the loss is quadratic between -sigma and sigma and extends linearly
## for tau>0, the central region is still quadratic but the formerly linear part
##       now gradually flattens out to a constant (tau is a half-life of the slope)
msmooth=function(x,y,df0=20,sigma=0.1,tau=10,iter=10,ww=rep(1,length(x))) {

  resid.y=y

  if (tau <= 0 ) {
	  for (j in 1:iter) {
	    w=pmin(1,sigma/abs(resid.y))*ww
	    s1=smooth.spline(x,y,w,df=df0)
	    resid.y=y-predict(s1,x)$y
	  }
  }
  if (tau>0) {
	  for (j in 1:iter) {
	    w=pmin(1,sigma*2^(-(abs(resid.y)-sigma)/tau)/abs(resid.y))*ww
	    s1=smooth.spline(x,y,w,df=df0)
	    resid.y=y-predict(s1,x)$y
	  }
	}
  s1
  }

## A subset of features is used to build the warp function. The subset is
## either formed by intensity threshold (which can be problematic) or by
## picking the most intense N from each run.
align_pair=function(P,df=20,sigma2=20,sigma=100) {
	s1=P$source;
	s2=P$dest;

	coef=linear_align(s1,s2,sigma)
	s1t=coef[1]+coef[2]*s1
	sm1=msmooth(s1t,s2-s1t,df,sigma2,sigma2)
	sm=function(x) {
		xt=coef[1]+coef[2]*x
		xtt=predict(sm1,xt)$y+xt
		}
	if (0) {
		plot(s1t,s2-s1t)
		lines(sm1,col='red')
		o=order(s1t)
		os1=s1t[o]
		os2=s2[o]
		lines(os1,runmed(os2-os1,51))
		lines(msmooth(s1t,s2-s1t,15,20,0),col='blue')
		}
	sm
	}


