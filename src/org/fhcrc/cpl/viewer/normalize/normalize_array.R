#################### function to normalize peptide array

###############################################
####### Run normalization on an array of intensities (columns 1:7 must have been stripped out!)
###############################################
driveNormalization<-function(pepArray, filterCut1=NULL, filterCut2=NULL, PlugMissing=TRUE, label=NULL, iter.min=NULL)
{
  n<-nrow(pepArray)
  m<-ncol(pepArray)

  pepArray[is.na(pepArray)]<-0
  f.c<-apply(pepArray!=0, 1, sum)

################## get filter threshold for intensity
  if(is.null(filterCut1)) filterCut1<-max(floor((m+1)/4),1) # features aligned in 1/4 runs considered good
  curArray<-pepArray[f.c>filterCut1,]
  curpool<-curArray[curArray!=0]
  filter<-quantile(curpool, 0.01)

################### filter
  PA.original<-pepArray # preserve the old data

  pepArray[pepArray<filter]=0
  f.c.old<-apply(pepArray!=0, 1, sum)
  PA.f<-pepArray[f.c.old>0,] # use filtered result to estimate normalization parameters
  
################### do normalization
  temp<-normalize.array.new(PA.f, cut.fc1=filterCut1, cut.fc2=filterCut2)
  PA.norm<-temp$PeptideArray
  PA.scale<-temp$scale

#################### if not doing missing value plug-in, just scale and return normalized array
  if(!PlugMissing)
     {
     scal.m<-matrix(PA.scale, nrow=nrow(pepArray), ncol=ncol(pepArray), byrow=T)
     result<-PA.original*scal.m  # apply normlization to original array
     return(result)
     }

##################### otherwise, perform missing value plug in
  f.v<-rep(filter, m)* exp(PA.scale) 
  feature.pick<-apply(PA.norm>0, 1, sum)>filterCut1 ### using good feature to look for KNN
  if(is.null(label)) ### unsupervised 
     {
       
       cur.data<-PA.norm
       # feature.pick<-apply(PA.norm>0, 1, sum)>filterCut1 ### using good feature to look for KNN
       
       if(is.null(iter.min)) 
             iter.min<-min(m, floor((m/2)+1)) ### the smallest number of nearest neighborhood number     
        
       for(kk in m:iter.min)
         {
           cur.dist<-dist(t(normlines(log(1+cur.data[feature.pick,]))))### work at log scale, and normalize each line
           d0<-as.matrix(cur.dist)
           cur.data<-KNNplugin(PA.norm, d0, f.v=f.v, lumda=1, k=kk) ### plug in missing with the nearest kk samples
          }
        PA.nm<-cur.data

    } else{ ### supervised 
         label.set<-unique(label)
         PA.nm<-PA.norm        
         for(i in label.set)
           {
              cur.pick<-label==i
              if(sum(cur.pick)>1) ### plug in missing only when a group has more than one sample
               {
                 cur.data<-PA.norm[, cur.pick]
                 cur.dist<-dist(t(normlines(log(1+cur.data[feature.pick,])))) 
                 d0<-as.matrix(cur.dist)
                 PA.nm[, cur.pick]<-KNNplugin(cur.data, d0, f.v=f.v, lumda=1, k=ncol(cur.data))
               }   
            } 
    }
 
  ######## apply computed scale factors to the original pepArray
   scal.m<-matrix(PA.scale, nrow=nrow(pepArray), ncol=ncol(pepArray), byrow=T)
   PA.result<-PA.original*scal.m  # apply normlization to original array
   PA.result[f.c.old>0, ]<-PA.nm  # plug in missing only for filtered features.

   return(PA.result)  
}

#######################


normlines<-function(xm)
{
   scale<-apply(xm, 1, nonzero.median)
   scale.m<-matrix(scale, nrow=nrow(xm), ncol=ncol(xm), byrow=F)
   xm/scale.m
}

nonzero.median<-function(xv)
{
  if(sum(xv!=0, na.rm=T)>0)
    return(median(xv[xv!=0]))
  return(NA)
}



###############################################
####### normalization function
###############################################
normalize.array.new<-function(PeptideArray, cut.fc1=NULL, cut.fc2=NULL, cut.ave.med=NULL,first=NULL,downquan=0.75, upquan=0.95, filt=0)
{
  n<-nrow(PeptideArray)
  m<-ncol(PeptideArray)
 
  ave.med<-apply(PeptideArray, 1, my.nonzero.median)
  f.c<-apply(PeptideArray!=0, 1, sum)
  
####################### for random missing
 
  if(is.null(cut.fc1)) cut.fc1=max(floor((m+1)/4),1) # ????? was cut.fc1=floor((m+1)/2)
  if(is.null(cut.ave.med)) cut.ave.med<-quantile(ave.med, 0.925)
 
  pick<-f.c>cut.fc1 & ave.med>cut.ave.med # ???? default was floor((m+1)/2)
  
  RM.alpha<-apply(PeptideArray[pick,]!=0, 2, sum)/sum(pick)
 
####################### for intensity dependent missing
 
  if(is.null(cut.fc2)) cut.fc2=max(floor((m+1)/5),1) ### added new parameter.
  PA.order<-apply(PeptideArray[f.c>cut.fc2, ], 2, sort)

  miss.count<-apply(PeptideArray>0.001, 2, sum)
  min.n<-min(miss.count)
 
  mean.logarray<-1:ncol(PA.order)
  n.sub<-nrow(PA.order)
  for(i in 1:ncol(PA.order))
    {
      if(!is.null(first))
        {
          cur.first<-floor(first*RM.alpha[i])
          pick<-(n.sub-first+1):n
        } else {
          cur.n<-min.n*RM.alpha[i]    
          pick<-floor(n.sub-cur.n*(1-downquan)):floor(n.sub-cur.n*(1-upquan))
        }
      mean.logarray[i]<-median(log(PA.order[pick,i]))############ take the largest peaks for normalization
    }
 
##################################
  mean.end<-mean(mean.logarray)
  alpha<-exp(mean.end-mean.logarray)
  alpha.matrix<-matrix(alpha, nrow=n, ncol=m, byrow=T)
  result<-PeptideArray * alpha.matrix
  result[result<filt]<-0
  return(list(PeptideArray=result, scale=alpha, RM.alpha=RM.alpha))
}

################################
my.nonzero.median.old<-function(x)
{
  y<-x[x>0]
  mean(y)
}

################################
my.nonzero.median<-function(x)
{
  return(nonzero.mean(x)) # Just call nonzero.mean ....
}

################################
nonzero.mean<-function(x)
{
  temp<-x[x!=0]
  if(length(temp)==0)
    return(0)
  return(mean(temp))
}

#################################
Missing.plugin<-function(norm.array, f.v)
{
  xdist<-as.matrix(dist(t(norm.array)))
  K=ncol(norm.array)
  result<-KNNplugin(norm.array, xdist, f.v, k=K)
  result
}

####################################
KNNplugin<-function(x, xdist, f.v, lumda=1, k)
{
  n<-nrow(x)
  m<-ncol(x)
  d<-xdist
  f.m<-matrix(f.v, nrow=n, ncol=m, byrow=TRUE)

  gm<-x
  obs<-x
  
  for(i in 1:m)
    {
      neighbor<-min.k(d[i,], k)
      temp<-x[, neighbor]
      temp1<-apply(temp,1, nonzero.mean)  ### estimate peak value
      gm[,i]<-temp1
    }

  tot<-x
  for(i in 1:m)
    {
      neighbor<-min.k(d[i,], k)
      temp<-x[, neighbor]
      gm.1<-gm[, neighbor]
      obs[,i]<-apply(temp!=0, 1, sum)     
      tot[,i]<-apply(gm.1>f.m[, neighbor] | temp!=0, 1, sum)
    }
  
  pp<-(obs+lumda)/(tot+lumda)  ### prob to have a peak
  pp[is.na(pp)]<-0
  result<-x
  result[gm<f.m]<- (pp*gm)[gm<f.m]
  result
}

min.k<-function(x,k)
{
  n<-length(x)
  cut<-sort(x)[k]
  (1:n)[x<=cut]
}

# ======================================================================

test1<-function()
  {
    # Pairwise alignment
    tmp <- read.delim('c:/data/ms2quant/RFP/plasma/2f.pepArray.tsv')
    inputIntensities <- tmp[,-(1:7)]
    normalizedIntensities <- driveNormalization(inputIntensities)
    write.table(normalizedIntensities, '2f.peiArray.tsv', quot=F, sep='\t', row.names=F)
  }

test2<-function()
  {
    # 10-way alignment
    tmp <- read.delim('c:/data/ms2quant/RFP/plasma/10f.pepArray.tsv')
    inputIntensities <- tmp[,-(1:7)]
    normalizedIntensities <- driveNormalization(inputIntensities)
    write.table(normalizedIntensities, '10f.peiArray.tsv', quot=F, sep='\t', row.names=F)
  }

test3<-function()
  {
    # Francisella
    tmp <- read.delim('y:/Projects/msInspectSupplement/alignment/FN.pepArray.tsv')
    inputIntensities <- tmp[,-(1:7)]
    normalizedIntensities <- driveNormalization(inputIntensities)
    write.table(normalizedIntensities, 'FN.peiArray.tsv', quot=F, sep='\t', row.names=F)
  }

FN.noplugin <-function()
  {
    # Francisella without missing value plugin
    tmp <- read.delim('FN-over1000.pepArray.tsv')
    inputIntensities <- tmp[,-(1:7)]
    normalizedIntensities <- driveNormalization(inputIntensities, PlugMissing=FALSE)
    write.table(normalizedIntensities, 'FN.peiArray.noplugin.tsv', quot=F, sep='\t', row.names=F)
  }

FN.pluginGroups <-function()
  {
    # Francisella without missing value plugin
    tmp <- read.delim('FN-over1000.pepArray.tsv')
    inputIntensities <- tmp[,-(1:7)]
    # First for runs are mutant, next are wild-type
    normalizedIntensities <- driveNormalization(inputIntensities, label=c(1,1,1,1,2,2,2,2))
    write.table(normalizedIntensities, 'FN.peiArray.groups.tsv', quot=F, sep='\t', row.names=F)
  }

inputIntensities <- read.delim('ArrayIntensities.tsv')
normalizedIntensities <- driveNormalization(inputIntensities, PlugMissing=FALSE)
write.table(normalizedIntensities, 'NormalizedIntensities.tsv', quot=F, sep='\t', row.names=F)
