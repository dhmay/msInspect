data1_1 = read.table(file = all_pep_text_file, header=T)
dataS_all = read.table(file = summary_text_file, header=T)
dataP_all_1 = read.table(file = all_prot_text_file, header=T)

data1 <- data1_1[substr(data1_1$Protein, 1,3) != "rev",]

dataP_all <- dataP_all_1[substr(dataP_all_1$Indistinguishable_Proteins, 1,3) != "rev",]

num_levels = length(levels(as.factor(data1$X.fraction)))

##############################################
#Legends
png(file = out_legend_file, height = 800, width = 600, pointsize=12, bg = "white")
plot(0,1, ann=FALSE, axes=FALSE)
title("legend")
n1 <- c(levels(as.factor(data1$X.fraction)))
n2 <- c(1:num_levels)
n3 <- paste(n2, n1, sep =": ")
legend("top", legend = n3, cex = 0.9, box.lty= 0)
dev.off()


##############################################
#delta mass
data1$FrDeltaMass <- data1$Delta_Mass
data1$FrDeltaMass <- ifelse(data1$Delta_Mass > 0.5 & data1$Delta_Mass <= 1.5, data1$Delta_Mass-1, data1$FrDeltaMass)
data1$FrDeltaMass <- ifelse(data1$Delta_Mass > 1.5 & data1$Delta_Mass <= 2.5, data1$Delta_Mass-2, data1$FrDeltaMass)
data1$FrDeltaMass <- ifelse(data1$Delta_Mass > 2.5 & data1$Delta_Mass <= 3.5, data1$Delta_Mass-3, data1$FrDeltaMass)
data1$FrDeltaMass <- ifelse(data1$Delta_Mass >-1.5 & data1$Delta_Mass <= -0.5, data1$Delta_Mass+1, data1$FrDeltaMass)
data1$FrDeltaMass <- ifelse(data1$Delta_Mass >-2.5 & data1$Delta_Mass <= -1.5, data1$Delta_Mass+2, data1$FrDeltaMass)
data1$FrDeltaMass <- ifelse(data1$Delta_Mass >-3.5 & data1$Delta_Mass <= -2.5, data1$Delta_Mass+3, data1$FrDeltaMass)

data1$FrDeltaMassppm <- data1$FrDeltaMass * 1000000 / data1$Precursor_Mass;

levels(data1$X.fraction) = c(1:num_levels)

plot_medium <- function(dataF, name1, name2)
{
	bins = max(as.numeric(dataF[,name1]))
	y = as.factor(dataF[,name1])
	z<-matrix(data=NA, nrow=bins, ncol=1)
	for (i in 1:bins) {
	z[i,]<-quantile(dataF[y==i,name2], probs=c(0.50))
	}
	z2 = c(NA)
	for (i in 1:bins){
	z2[i] <- median(as.numeric(dataF[y==i,name1]),na.rm=T)
	}
	for (i in 1:1) {
	lines(z2,z[,i], col=i+2)
	}
}

png(file = delta_mass_plots_file, height = 800, width = 1000, pointsize=12, bg = "white")
op = par(mfrow=c(2,1))
plot(as.numeric(data1$X.fraction), data1$Delta_Mass, xaxt="n", xlab = "fractions", ylab = "Delta_Mass (Da)",  main = "Delta Mass Per Fraction",pch = 19,cex=0.3)
legend("topleft",col="green",lty=1,legend="median")
axis(1,xaxp=c(1,100,99))
abline(h=-1, col=2)
abline(h=0, col=2)
abline(h=1, col=2)
abline(h=2, col=2)
abline(h=3, col=2)
data_s1 = data1[data1$Delta_Mass < -0.5 & data1$Delta_Mass >=-1.5,]
if(nrow(data_s1) > 0) {plot_medium(data_s1, "X.fraction", "Delta_Mass")}
data_s2 = data1[data1$Delta_Mass <= 0.5 & data1$Delta_Mass >=-0.5,]
if(nrow(data_s2) > 0) {plot_medium(data_s2, "X.fraction", "Delta_Mass")}
data_s3 = data1[data1$Delta_Mass <= 1.5 & data1$Delta_Mass >=0.5,]
if(nrow(data_s3) > 0) {plot_medium(data_s3, "X.fraction", "Delta_Mass")}
data_s4 = data1[data1$Delta_Mass <= 2.5 & data1$Delta_Mass >=1.5,]
if(nrow(data_s4) > 0) {plot_medium(data_s4, "X.fraction", "Delta_Mass")}
data_s5 = data1[data1$Delta_Mass <= 3.5 & data1$Delta_Mass >=2.5,]
if(nrow(data_s5) > 0) {plot_medium(data_s5, "X.fraction", "Delta_Mass")}

####plot fractional delta mass
plot(as.numeric(data1$X.fraction), data1$FrDeltaMassppm, xaxt="n", xlab = "fractions", ylab = "Fractional_Delta_Mass (ppm)",  main = "Fractional Delta Mass Per Fraction",pch = 19,cex=0.3)
legend("topleft",col="green",lty=1,legend="median")
axis(1,xaxp=c(1,100,99))
abline(h=0, col=2)
plot_medium(data1, "X.fraction", "FrDeltaMassppm")
par(op)
dev.off()

png(file = delta_mass_boxplots_file, height = 800, width = 1000, pointsize=12, bg = "white")
boxplot(data1$FrDeltaMassppm~as.factor(data1$X.fraction), xlab = "fractions", ylab = "Fractional_Delta_Mass (ppm)",  main = "Fractional Delta Mass Per Fraction")
abline(h=0, col = 2)
dev.off()


###############################################
## plots for peptides/proteins/genes
dataS = dataS_all[dataS_all$X.fractions != "total",]
dataS$number <- c(1:length(dataS$X.fractions))
r = nrow(dataS)

png(file = peptides_file, height = 600, width = 800, pointsize=12, bg = "white")
plot(dataS$total_peptide~dataS$number, xaxt="n", xlab = "fractions", ylab = "counts", ylim=c(0, max(dataS$total_peptide)), col="white")
lines(dataS$number, dataS$total_peptide, lty=1, lwd=2, col=1)
lines(dataS$number, dataS$unique_peptide, lty=2, lwd=2, col=2)
lines(dataS$number, dataS$total_quantified_peptide, lty=3, lwd=2, col=3)
lines(dataS$number, dataS$unique_quantified_peptide, lty=4, lwd=2, col=4)
n1 <- c("total peptides", "unique peptides", "total quantified peptides", "unique quantified peptides")
n2 <- c(dataS_all[r+1, 2], dataS_all[r+1,3], dataS_all[r+1,6], dataS_all[r+1,7])
n3 <- paste(n1, n2, sep ="=")
legend("topleft",legend = n3, lty=c(1,2,3,4), bty='n', col=1:4, lwd=rep(2,4), cex=0.65)
axis(1,xaxp=c(1,100,99))
dev.off()


png(file = gene_protein_groups_file, height = 600, width = 800, pointsize=12, bg = "white")
if(max(dataS$protein_groups) > max(dataS$gene_symbols))
{plot(dataS$protein_groups~dataS$number, xaxt="n", xlab = "fractions", ylim=c(0, max(dataS$protein_groups)),ylab = "counts",col="white")}
if(max(dataS$protein_groups) <= max(dataS$gene_symbols))
{plot(dataS$gene_symbol~dataS$number, xaxt="n", xlab = "fractions", ylab = "counts", ylim=c(0, max(dataS$gene_symbol)),col="white")}
lines(dataS$number, dataS$protein_groups, lty=1, lwd=2, col=1)
lines(dataS$number, dataS$gene_symbols, lty=2, lwd=2, col=2)
lines(dataS$number, dataS$quantified_proteins, lty=3, lwd=2, col=3)
lines(dataS$number, dataS$quantified_genes, lty=4, lwd=2, col=4)
n1 <- c("protein groups", "genes", "quantified protein groups", "quantified genes")
n2 <- c(dataS_all[r+1, 4], dataS_all[r+1,5], dataS_all[r+1,8], dataS_all[r+1,9])
legend_titles <- paste(n1, n2, sep ="=")
legend("topleft",legend = legend_titles, lty=c(1,2,3,4), bty='n', col=1:4, lwd=rep(2,4), cex=0.65)
axis(1,xaxp=c(1,100,99))
dev.off()

###############################################
##quantification plots for peptides
dataA1 <- data1[data1$Decimal_Ratio != -666,]
dataA <- dataA1[!(dataA1$Light_Area == 0 & dataA1$Heavy_Area == 0),]
dataA$Ratio = replace(dataA$Decimal_Ratio, ((dataA$Light_Area <= 0 & dataA$Heavy_Area > 0) | dataA$Decimal_Ratio <= 0), 0.001)
dataA$Ratio = replace(dataA$Decimal_Ratio, (dataA$Light_Area > 0 & dataA$Heavy_Area <= 0), 999)
dataA$logL2H =log2(dataA$Ratio)
dataA$logLH = (log10(dataA$Light_Area*dataA$Heavy_Area))/2.0

if(nrow(dataA) > 0){

plot_quantile <- function(dataF, name1, name2, name3)  #name1=x value, name2=y value
{
	if(nrow(dataF) > 0){
	plot(dataF[,name1],dataF[,name2], xlab = name1, ylab = name2, main = name3,pch = 19,cex=0.3)}

	if(nrow(dataF) > 10)
	{
		temp <- rank(dataF[,name1])
		len <- length(temp)
		bins <- 1
		if(len >= 2000){	bins = 20}
		if(len >=1000 & len < 2000) {bins = 15}
		if(len >= 500 & len < 1000){bins = 8}
		if(len >= 100 & len < 500){bins = 4}
		if(len < 100){bins = 2}

		ctemp = cut(temp, breaks = bins)
		levels(ctemp) = c(1:bins)
		y = factor(ctemp)

		z<-matrix(data=NA,nrow=bins+2, ncol=5)

		for (i in 1:bins) {
		z[i+1,]<-quantile(dataF[y==i,name2], probs=c(0.05,0.25,0.50,0.75,0.95))
		}
		z[1,] = quantile(dataF[y==1,name2], probs=c(0.05,0.25,0.50,0.75,0.95))
		z[bins+2,] = quantile(dataF[y==bins, name2], probs=c(0.05,0.25,0.50,0.75,0.95))

		z2 = c(min(dataF[,name1]))
		for (i in 1:bins){
		z2 = c(z2,median(dataF[y==i,name1]))
		}
		z2[bins+2] = c(max(dataF[,name1]))

		lines(z2,z[,1], col=3)
		lines(z2,z[,2], col=4)
		lines(z2,z[,3], col=2)
		lines(z2,z[,4], col=4)
		lines(z2,z[,5], col=3)
		legend("topright",legend = c("median","mean","quater","5 percent"), lty=1,col=c(2,5,4,3), cex=0.65)

		abline(h=0)
		vv = summary(dataF[,name1])
		abline(v=vv[2], col=4)
		abline(v=vv[3], col=2)
		abline(v=vv[4], col=5)
		abline(v=vv[5], col=4)
	}
}

png(file = quantified_peps_file, height = 600, width = 800, pointsize=12, bg = "white")
op = par(mfrow=c(1,2))
plot_quantile(dataA, "logLH", "logL2H", "MA plot")
plot_quantile(dataA, "Precursor_Mass", "logL2H", "ratio vs. mass")
par(op)
dev.off()


#this doesn't currently appear to be working... why?

num_slides = floor(num_levels/4)
if(num_slides>=1){
for (i in 1:num_slides)
{
	file_name = paste(qa_dir,"quantified_peps_", sep="")
	file_name = paste(file_name,i, sep="")
#file_name = paste("quantified_peps_", i, sep="")
	file_name = paste(file_name,".png", sep="")
	png(file = file_name, height = 600, width = 800, pointsize=12, bg = "white")
	op = par(mfrow=c(2,2))
	for (j in 1:4)
	{
		frxn = (i-1)*4+j
		dataF <- dataA[as.integer(dataA$X.fraction) == frxn,]
		plot_name = paste("MA_plot_", frxn, sep="")
		plot_quantile(dataF, "logLH", "logL2H", plot_name)
	}
	par(op)
	dev.off()
}
}

if(num_levels > num_slides*4)
{
file_name = paste(qa_dir,"quantified_peps_", sep="")
file_name = paste(file_name,num_slides+1, sep="")
file_name = paste(file_name,".png", sep="")
png(file = file_name, height = 600, width = 800, pointsize=12, bg = "white")
op = par(mfrow=c(2,2))
for (j in 1:(num_levels-num_slides*4))
{
	frxn = num_slides*4+j
	dataF <- dataA[as.integer(dataA$X.fraction) == frxn,]
	plot_name = paste("MA_plot_", frxn, sep="")
	plot_quantile(dataF, "logLH", "logL2H", plot_name)
}
par(op)
dev.off()
}

##################################################
###########log ratio versus num of labels############
dataA$num_labels <- round((dataA$Heavy_Mass - dataA$Light_Mass)/label_mass_diff)

plot_quantile2 <- function(dataF, name1, name2)
{
	y = factor(dataF[,name1])
	bins = max(dataF[,name1])
	z<-matrix(data=NA,nrow=bins, ncol=5)
	for (i in 1:bins) {
	z[i,]<-quantile(dataF[y==i,name2], probs=c(0.05,0.25,0.50,0.75,0.95))
	}
	z2 = c(NA)
	for (i in 1:bins){
	z2[i] = median(dataF[y==i,name1],na.rm=T)
	}
	plot(dataF[,name1],dataF[,name2], xlab = name1, ylab = name2,pch = 19,cex=0.3)
	lines(z2,z[,1], col=3)
	lines(z2,z[,2], col=4)
	lines(z2,z[,3], col=2)
	lines(z2,z[,4], col=4)
	lines(z2,z[,5], col=3)
	legend("topright",legend = c("median","quater","5 percent"), lty=1,col=c(2,4,3), cex=0.65)
	abline(h=0)
}
png(file = quantified_peps_numlabels_file, height = 600, width = 800, pointsize=12, bg = "white")
plot_quantile2(dataA, "num_labels", "logL2H")
dev.off()




####################Generate spread sheets for quantile summary
write_quantile <- function(dataF, name1, name2, File)  #name1=x value, name2=y value
{
	temp <- rank(dataF[,name1])
	bins <- 3
	ctemp = cut(temp, breaks = bins)
	levels(ctemp) = c(1:bins)
	yy = factor(ctemp)
	zz<-matrix(data=NA,nrow=bins+1, ncol=9)
	for (i in 1:bins) {
	zz[i,]<-quantile(dataF[yy==i,name2], probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
	}
	zz[bins+1,] <- quantile(dataF[,name2], probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
	zz2 <- c(NA)
	for (i in 1:bins){
	zz2[i] <- median(dataF[yy==i,name1],na.rm=T)
	}
	zz2[bins+1] <- median(dataF[,name1],na.rm=T)
	n1 <- c("low", "mid", "high", "all")
	catagory <- paste(n1, name1, sep = "_")
	data_points <- as.data.frame(table(ctemp))
	frequency <- data_points$Freq
	frequency[bins+1] <- nrow(dataF)
	mid_value <- round(zz2,2)
	percent90 <- paste(round(2^zz[,1],2), round(2^zz[,9],2), sep = " - ")
	percent75 <- paste(round(2^zz[,2],2), round(2^zz[,8],2), sep = " - ")
	percent50 <- paste(round(2^zz[,3],2), round(2^zz[,7],2), sep = " - ")
	percent25 <- paste(round(2^zz[,4],2), round(2^zz[,6],2), sep = " - ")
	medium_ratio <- round(2^zz[,5],2)
	mid_value_name <- paste("median_", name1,sep = "_")
	dataW1 <-data.frame(catagory, frequency, mid_value, medium_ratio, percent25, percent50, percent75, percent90)
	colnames(dataW1)[3] <- mid_value_name
	write.csv(dataW1, File, row.names = FALSE)
}

write_quantile2 <- function(dataF, name1, name2, File)  #name1=x value, name2=y value
{
	y = factor(dataF[,name1])
	bins = max(dataF[,name1])
	zz<-matrix(data=NA,nrow=bins+1, ncol=9)
	for (i in 1:bins) {
	zz[i,]<-quantile(dataF[y==i,name2], probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
	}
	zz[bins+1,] <- quantile(dataF[,name2], probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
	zz2 = c(NA)
	for (i in 1:bins){
	zz2[i] = median(dataF[y==i,name1],na.rm=T)
	}
	zz2[bins+1] <- median(dataF[,name1],na.rm=T)
	n1 <- c(1:bins)
	catagory <- c(paste(n1, name1, sep = "_"), "all")
	data_points <- as.data.frame(table(dataF[,name1]))
	frequency <- data_points$Freq
	frequency[bins+1] <- nrow(dataF)
	percent90 <- paste(round(2^zz[,1],2), round(2^zz[,9],2), sep = " - ")
	percent75 <- paste(round(2^zz[,2],2), round(2^zz[,8],2), sep = " - ")
	percent50 <- paste(round(2^zz[,3],2), round(2^zz[,7],2), sep = " - ")
	percent25 <- paste(round(2^zz[,4],2), round(2^zz[,6],2), sep = " - ")
	medium_ratio <- round(2^zz[,5],2)
	dataW1 <-data.frame(catagory, frequency, medium_ratio, percent25, percent50, percent75, percent90)
	write.csv(dataW1, File, row.names = FALSE)
}

write_quantile(dataA, "Precursor_Mass", "logL2H", quantified_mass_summary_file)
write_quantile(dataA, "logLH", "logL2H", quantile_intensity_summary_file)
write_quantile2(dataA, "num_labels", "logL2H",quantile_cys_summary_file)


}## end of if nrow > 0


##############################################
#quantification plots for proteins
dataP = dataP_all[dataP_all$L2H_Mean != -666,]
dataP$logG = log2(dataP$Ratio_Peps)
dataP$logR = log2(dataP$L2H_Mean)

if(nrow(dataP) > 0){
ctemp = cut(dataP$Ratio_Peps, breaks = c(0, 1, 2, 3, 4, 6, 10, 20, max(dataP$Ratio_Peps)))
len = length(table(ctemp))
bins = len
levels(ctemp) = c(1:len)
y = factor(ctemp)
z<-matrix(data=NA,nrow=bins+1, ncol=5)
for (i in 1:bins) {
z[i,]<-quantile(dataP[y==i,]$logR, probs=c(0.05,0.25,0.50,0.75,0.95))
}
z[bins+1,] <- z[bins,]
z2 = c(NA)
for (i in 1:bins){
z2[i] = median(dataP[y==i,]$logG,na.rm=T)
}
z2[bins+1] = max(dataP$logG)
png(file = quantified_protein_groups_file, height = 600, width = 800, pointsize=12, bg = "white")
op = par(mfrow=c(1,2))
plot(dataP$logG, dataP$logR, xlab = "log2(#PeP in protein group)", ylab = "log2L/H")
lines(z2,z[,1], col=3)
lines(z2,z[,2], col=4)
lines(z2,z[,3], col=2)
lines(z2,z[,4], col=4)
lines(z2,z[,5], col=3)
legend("topright",legend = c("median","quater","5 percent"), lty=1,col=c(2,4,3), cex=0.65)
abline(h=0)
plot(density(log2(dataP$L2H_Mean)), main="log2(L2H Mean)")
par(op)
dev.off()

### write a summary spread sheet
ctemp1 = cut(dataP$Ratio_Peps, breaks = c(0, 1, 3, 6, 10, max(dataP$Ratio_Peps)))
len1 = length(table(ctemp1))
bins1 = len1
levels(ctemp1) = c(1:len1)
y = factor(ctemp1)
z3 <-matrix(data=NA,nrow=bins1+1, ncol=9)
for (i in 1:bins1) {
z3[i,]<-quantile(dataP[y==i,]$logR, probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
}
z3[bins1+1,] <- quantile(dataP$logR, probs=c(0.05,0.125,0.25,0.375,0.50,0.625,0.75,0.875,0.95))
z4 = c(NA)
for (i in 1:bins1){
z4[i] = median(dataP[y==i,]$Ratio_Peps,na.rm=T)
}
z4[bins1+1] <- median(dataP$Ratio_Peps,na.rm=T)
mid_value <- round(z4,1)
catagory <- c("1 peptides", "2~3 peptides","4~6 peptides","7~10 peptides",">10 peptides" ,"all")
data_points <- as.data.frame(table(ctemp1))
frequency <- data_points$Freq
frequency[bins1+1] <- nrow(dataP)
percent90 <- paste(round(2^z3[,1],2), round(2^z3[,9],2), sep = " - ")
percent75 <- paste(round(2^z3[,2],2), round(2^z3[,8],2), sep = " - ")
percent50 <- paste(round(2^z3[,3],2), round(2^z3[,7],2), sep = " - ")
percent25 <- paste(round(2^z3[,4],2), round(2^z3[,6],2), sep = " - ")
medium_ratio <- round(2^z3[,5],2)
mid_value_name <- "median_num_peps"
dataW1 <-data.frame(catagory, frequency, mid_value, medium_ratio, percent25, percent50, percent75, percent90)
colnames(dataW1)[3] <- mid_value_name
write.csv(dataW1, "quantile_summary_protein_groups.csv", row.names = FALSE)
}



7;
