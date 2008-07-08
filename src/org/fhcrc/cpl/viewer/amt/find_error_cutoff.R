#proportion_of_max_delta=0.25;
#errordata<-read.delim("/home/dhmay/projects/AMT/data_driven_tolerances/masserror_data.tsv");
#showcharts<-TRUE;
#chart_out_file<-("/home/dhmay/temp/charts.png");

#requires parameters to be set:
#proportion_of_max_delta:  the proportion of the cutoff delta-density value
#  to the maximum delta-density
#errordata: the data to find cutoffs for
#showcharts: show charts?  TRUE or FALSE
#chart_out_file: output chart file, if showing charts
dens<-density(errordata);
delta<-dens$y[-1]-dens$y[-length(dens$y)];
cut<-max(delta)*proportion_of_max_delta;

fullx<-dens$x[-length(dens$x)];

py<-delta[(fullx < fullx[delta==max(delta)])];
px<-fullx[(fullx < fullx[delta==max(delta)])];
p<-rep(T, length(py));
for (i in 1:length(py))
 { p[i]<-all(py[i:length(py)]>(cut)) };
leftx<-px[px==min(px[p])];
rm(py, px, p);

py<- delta[(fullx > fullx[delta==min(delta)])];
px<-fullx[(fullx > fullx[delta==min(delta)])];
p<-rep(F, length(py));
for (i in 1:length(py))
 {p[i]<-all(py[1:i]< -(cut)) };
rightx<-px[px==max(px[p])];
rm(py, px, p);

cutoff<-c(leftx, rightx);

if (showcharts)
{
  png(chart_out_file);
  par(mfrow=c(3,1));
  
  hist(errordata, breaks=50, xlab="X", main="");
  plot(dens$x, dens$y, pch=20, xlab="X", ylab="Density");
  abline(v=cutoff[1], col=2);
  abline(v=cutoff[2], col=2);
  
  delta<-dens$y[-1]-dens$y[-length(dens$y)];
  plot(dens$x[-length(dens$x)], delta, pch=20, type="p", xlab="X", ylab="Delta.Density");
  abline(v=cutoff[1], col=2);
  abline(v=cutoff[2], col=2);
  dev.off();
}

cutoff;
