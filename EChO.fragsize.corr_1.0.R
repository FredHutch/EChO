#!/usr/bin/Rscript

rm(list=ls())
## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 4) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--first=somevalue	- Input line to be compared
			--second=someValue   - Data frame to compare to
			--comb=someValue   - Line number
			--span=someValue	- Output prefix
")
 
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
 
## Arg1 default
#if(is.null(args[1])){
if(is.null(argsL$first) | is.null(argsL$second) | is.null(argsL$int) | is.null(argsL$output)) {
  stop("Argument is missing!
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--first=somevalue	- Input line to be compared
			--second=someValue   - Data frame to compare to
			--comb=someValue   - Line number
			--span=someValue	- Output prefix
")
 
  q(save="no")
}

first<-read.table(argsL$first)
second<-read.table(argsL$second)
int<-as.vector(read.table(argsL$int))
index<-seq(as.numeric(int[2]),as.numeric(int[3]),1)
indexfirst<-seq(as.numeric(int[5]),as.numeric(int[6]),1)
indexsecond<-seq(as.numeric(int[9]),as.numeric(int[10]),1)
lofirst<-loess(first$V2~first$V1, span=first$V3[1])
outfirst<-predict(lofirst, indexfirst, se=TRUE)
outfirstfit<-outfirst$fit[indexfirst>=min(index) & indexfirst<=max(index)]
outfirstse<-outfirst$se.fit[indexfirst>=min(index) & indexfirst<=max(index)]
losecond<-loess(second$V2~second$V1, span=second$V3[1])
outsecond<-predict(losecond, indexsecond, se=TRUE)
outsecondfit<-outsecond$fit[indexsecond>=min(index) & indexsecond<=max(index)]
outsecondse<-outsecond$se.fit[indexsecond>=min(index) & indexsecond<=max(index)]
match<-na.omit(data.frame(first=outfirstfit, firstse=outfirstse, second=outsecondfit, secondse=outsecondse, index=index))
if(nrow(match) > 401){
	outfirst2<-match[c(200:(length(match$first)-200)),c(1:2)]
	outfirstmin<-min(outfirst2$first)[1]
	outfirstminse<-outfirst2$firstse[outfirst2$first==min(outfirst2$first)][1]
	minfirst<-match$index[outfirst2$first==min(outfirst2$first)][1]
	#head(minfirst)
	outsecond2<-match[c(200:(length(match$second)-200)),c(3:4)]
	outsecondmin<-min(outsecond2$second)[1]
	outsecondminse<-outsecond2$secondse[outsecond2$second==min(outsecond2$second)][1]
	minsecond<-match$index[outsecond2$second==min(outsecond2$second)][1]
	#head(minsecond)
	matchfirst<-match[match$index >= minfirst-200 & match$index <= minfirst+200,]
	matchsecond<-match[match$index >= minsecond-200 & match$index <= minsecond+200,]
	corrfirst<-as.vector(cor(matchfirst[,c(1,3)]))[2]
	corrsecond<-as.vector(cor(matchsecond[,c(1,3)]))[2]
	ccfirst<-ccf(matchfirst$first, matchfirst$second, lag.max=ceiling(length(index)/2), plot=FALSE)
	maxcorrfirst<-max(ccfirst$acf)[1]
	maxshiftfirst<-ccfirst$lag[ccfirst$acf==max(ccfirst$acf)][1]
	ccsecond<-ccf(matchsecond$first, matchsecond$second, lag.max=ceiling(length(index)/2), plot=FALSE)
	maxcorrsecond<-max(ccsecond$acf)[1]
	maxshiftsecond<-ccsecond$lag[ccsecond$acf==max(ccsecond$acf)][1]
	all<-rbind(c(as.vector(int[1,1]),minfirst,minfirst+1,outfirstmin,outfirstminse,corrfirst,maxcorrfirst,maxshiftfirst),c(as.vector(int[1,1]),minsecond,minsecond+1,outsecondmin,outsecondminse,corrsecond,maxcorrsecond,maxshiftsecond))
	write.table(all, file=paste(argsL$output, ".crosscor.txt", sep=""), append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	#write.table(all, file="", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
