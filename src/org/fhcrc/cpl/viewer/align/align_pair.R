map_file_scans_onto_file=function(fname1,fname2,fname3,bigIntensity=100,topN=NA) {
    ## dhmay creating 9-7-2006
    ## Create a mapping of file 1's scans to match file 2, producing file 3
	f1=readf(fname1)
	f2=readf(fname2)
	warp=align_pair(f1, f2, bigIntensity=bigIntensity, topN=topN)
	source_scan = (c(1:max(f1$scan)))
	dest_scan = warp(source_scan)
	dest_scan = round(dest_scan)
	scan_map = cbind(source_scan, dest_scan)
	write.table(scan_map,fname3,row.names=F,quote=F,na="",sep="\t")
    }

align_file_onto_file=function(fname1,fname2,fname3,bigIntensity=100,topN=NA) {
	## Align file 1's scans to match file 2 producing file 3
	f1=readf(fname1)
	f2=readf(fname2)
	warp=align_pair(f1, f2, bigIntensity=bigIntensity, topN=topN)
	f3=read.table(fname1, header=T, fill=T,sep="\t")
	f3[,'scan']=warp(f3[,'scan'])
	f3[,'scan']=round(f3[,'scan'])
	write.table(f3,fname3,row.names=F,quote=F,na="",sep="\t")
	}

readf=function(fn) {
	f=read.table(fn, header=T, fill=T, sep="\t")
	f=f[order(f[,'mz']),]
	}

get_pairs=function(f1,f2,mzcut=0.1) {
	mz1=f1[,'mz']
	mz2=f2[,'mz']

	mlist=numeric(0)
	i2top=1;I2=length(mz2)
	upseq=function(i,j) {if (j>=i) {i:j} else {c()}}
	for (i1 in 1:length(mz1)) {
		i2=i2top
		while (i2<=I2) {
			if (mz2[i2]<mz1[i1]-mzcut) {i2=i2top=i2+1;next}
			if (mz2[i2]>mz1[i1]+mzcut) break
			mlist=rbind(mlist,c(i1,i2))
			i2=i2+1
			}
		}
	mlist
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

# Return a vector of booleans that flag at *least* the top N entries in a list.
# This will return more than N entries when more than one entry shares the
# value of the boundry element; we do this to avoid making random exclusions
# of "equivalent" values. Note that in the degenerate case, this will return
# all TRUEs if all members of the input list have the same value.
topn <- function(x, n=300)
  {
    if (length(x) <= n)
      return(array(TRUE, c(length(x))))

    thresh = sort(x)[length(x) - n + 1]
    x >= thresh
  }

## A subset of features is used to build the warp function. The subset is
## either formed by intensity threshold (which can be problematic) or by
## picking the most intense N from each run.
align_pair=function(f1,f2,df=20,sigma2=20,sigma=100,bigIntensity=100,topN=NA) {
  
        if (is.na(topN))
          {
            g1=f1[ f1[,'intensity']>bigIntensity, ]
            g2=f2[ f2[,'intensity']>bigIntensity, ]
          }
        else
          {
            g1=f1[ topn(f1$intensity, topN), ]
            g2=f2[ topn(f2$intensity, length(g1$intensity)), ] # to keep the lengths comparable; may not be necessary
          }

	P=get_pairs(g1,g2)
	s1=g1[P[,1],'scan']
	s2=g2[P[,2],'scan']
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

if (0) { ## Testing code
	names=c('10e','11h')
	pattern='X:/EDI/ISB/ovarian/good_data/ovarin%s.mzXML.batPeptides.txt'
	filenames=sapply(names,function(x){sprintf(pattern,x)})

	align_file_onto_file(filenames[1],filenames[2],'test12.txt')

	fl=lapply(filenames,readf)
	f1=fl[[1]]
	f2=fl[[2]]

	f2[,'scan']=f2[,'scan']+0.2*pmax(f2[,'scan']-1500,0)
	sm=align_pair(f1,f2)
	plot(sm(s1),s2-sm(s1),ylim=c(-100,100))
}


