#file<-'consensus_intersect.peparray.tsv';
#minPerGroup<-2;
#runs_np<-c(6,10,11,1,9);
#runs_pneg<-c(2,12,5,8,13);
#runs_phgd<-c(4,7,14,3,15);
#
#cases<-runs_phgd;
#cases<-runs_pneg;
#controls<-runs_np;

library('qvalue');

ttest_array<-function(file, minPerGroup, cases, controls)
{
  arr = read.delim(file,comment.char='#')
  arr_intonly = arr[seq(from=9, to=ncol(arr), by=3)];

  #arr_intcase = arr_intonly[cases];
  #arr_intcontrol = arr_intonly[controls];

  hasEnoughVals <- function(x) {x_cases<-x[cases]; x_controls<-x[controls]; length(x_cases[!is.na(x_cases)]) >= minPerGroup && length(x_controls[!is.na(x_controls)]) >= minPerGroup;}
  nonEquivalentVals <- function(x) { length(unique(x[c(cases,controls)])) > 1; }

#arr_hasEnoughVals <- apply(arr_intonly, 1, hasEnoughVals)

#arr = arr[arr_hasEnoughVals,];
#arr_intonly = arr[seq(from=9, to=ncol(arr), by=3)];
  runTTest <- function(x)
  {
      rowcases<-x[cases];
      rowcontrols = x[controls];
      result=c(NA,NA);
      if (hasEnoughVals(x) & nonEquivalentVals(x))
      {
        tresult<-t.test(rowcases, rowcontrols);
        result=c(tresult$p.value, tresult$statistic);
if (result[2]==1) {print(rowcases)}
if (result[2]==1) {print(rowcontrols)}
if (result[2]==1) {print('.')}
      }
      result;
  }
#par(mfrow=c(3,1))
#print('running t-tests...');
  tresults<-apply(arr_intonly, 1, runTTest);
#print('done running t-tests.');
  pvalues<-tresults[1,];
  tstats<-tresults[2,];
#print('calculating q-values...');
  qvalues_nona<-qvalue(pvalues[!is.na(pvalues)])$qvalues;
  qindex=1;
  qvalues<-array(data=NA,dim=length(pvalues));
  for (i in 1:length(pvalues))
  {
    if (!is.na(pvalues[i]))
    {
     qvalues[i]<-qvalues_nona[qindex];
     qindex = qindex + 1;
    }
  }

  result=matrix(nrow=length(pvalues), ncol=3);
  result[,1]<-pvalues;
  result[,2]<-tstats;
  result[,3]<-qvalues;
  result;
}


