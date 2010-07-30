#This script reads a GSE (geneset enrichment) file of the format:
#peptide	ratio	tscore	intmeancase	intmeancontrol	protein1	protein2	...
#AMIAYWTNFAK     0.64924484      -1.0779228      97130.67        149605.61       0       1	...
#Using the peptide t-scores associated with each protein, it performs a Wilcox test, comparing that 
#protein's peptide t-scores with all other peptide t-scores in the file.  It calculates p-values and
#q-values for each protein, and writes the results to outfile

#requires q-value library to be loaded

#required variables:
#file: input file
#outfile: output file

gsea<-read.table(file,head=TRUE);
minpresent<-2;

library('qvalue');

ttestcol<-function(gsea, colnum, minpresent) 
{ 
  res<-NA;
  if (length(gsea$tscore[gsea[,colnum]==1]) >= minpresent) 
  {
    res<-t.test(gsea$tscore[gsea[,colnum]==1],gsea$tscore[gsea[,colnum]==0])$p.value;
  } 
  as.double(res);
}

wilcoxcol<-function(gsea,colnum, minpresent) 
{ 
  res<-NA;
  if (length(gsea$tscore[gsea[,colnum]==1]) >= minpresent) 
  {
    res<-wilcox.test(gsea$tscore[gsea[,colnum]==1],gsea$tscore[gsea[,colnum]==0])$p.value;
  } 
  as.double(res);
}

numinfocols<-5

numgenesets<-length(colnames(gsea))-numinfocols;
geneset_probs<-array(dim=c(numgenesets,6));
geneset_probs[,1]<-colnames(gsea)[(numinfocols+1):length(colnames(gsea))];

meanmeanints<-apply(gsea[,4:5],1,mean);
for (colnum in 1:numgenesets)
{
       geneset<-colnames(gsea)[colnum+numinfocols];
#        geneset_probs[colnum,2]<-geneset_genecounts[geneset_genecounts$geneset==geneset,]$numgenes;
        geneset_probs[colnum,2]<-length(unique(gsea[gsea[,colnum+numinfocols]==1,1]));
        geneset_probs[colnum,3]<-mean(gsea$tscore[gsea[,colnum+numinfocols]==1]);
	geneset_probs[colnum,4]<-ttestcol(gsea,colnum+numinfocols,minpresent);
	geneset_probs[colnum,5]<-wilcoxcol(gsea,colnum+numinfocols,minpresent);
       geneset_probs[colnum,6]<-mean(meanmeanints[gsea[,colnum+numinfocols]==1]);
}

#colnames(geneset_probs)<-c('geneset','numgenes','numobserved','mean_tscore','tp','wilcoxp','meanmeanint');
colnames(geneset_probs)<-c('geneset','numobserved','mean_tscore','tp','wilcoxp','meanmeanint');
write.table(geneset_probs, "temp.tsv",na="NA",sep="\t",row.names=FALSE,quote=FALSE);
geneset_probs<-read.table("temp.tsv",head=TRUE);

  qvalues_nona<-qvalue(geneset_probs$wilcoxp[!is.na(geneset_probs$wilcoxp)])$qvalues;
  qindex=1;
  qvalues<-array(data=NA,dim=length(geneset_probs$wilcoxp));
  for (i in 1:length(geneset_probs$wilcoxp))
  {
    if (!is.na(geneset_probs$wilcoxp[i]))
    {
     qvalues[i]<-qvalues_nona[qindex];
     qindex = qindex + 1;
    }
  }
geneset_probs[,7]=qvalues;
#colnames(geneset_probs)<-c('protein','numpeptides','numobserved','mean_tscore','tp','wilcoxp','meanmeanint','wilcoxq');
colnames(geneset_probs)<-c('protein','numobserved','mean_tscore','tp','wilcoxp','meanmeanint','wilcoxq');
write.table(geneset_probs, outfile,na="",sep="\t",row.names=FALSE,quote=FALSE);
#hist(geneset_probs$wilcoxq,breaks=30);
