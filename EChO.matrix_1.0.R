#!usr/bin/Rscript

rm(list=ls())
## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--frame=someValue   - Input offset vs. fragment size data frame
			--center=someValue   - Peak center
			--chr=someValue   - Peak chromosome
			--span=someValue   - span (auto-select with GCV if non-numeric entry)
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
if(is.null(argsL$frame) | is.null(argsL$center)) {
  stop("Argument is missing!
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--frame=someValue   - Input offset vs. fragment size data frame
			--center=someValue   - Peak center
			--chr=someValue   - Peak chromosome
			--line=someValue	- Line number being processed
			--span=someValue   - span (auto-select with GCV if non-numeric entry)
")
 
  q(save="no")
}

if(is.null(argsL$span)) {
  stop("Must supply span! Numeric for direct value or non-numeric for GCV-calculated span
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--frame=someValue   - Input offset vs. fragment size data frame
			--center=someValue   - Peak center
			--chr=someValue   - Peak chromosome
			--line=someValue	- Line number being processed
			--span=someValue   - span (auto-select with GCV if non-numeric entry)
")
 
  q(save="no")
}

frame<-read.table(argsL$frame)
center<-as.vector(argsL$center)
index<-seq(min(frame$V1), max(frame$V1), 1)
if(is.na(as.numeric(argsL$span))==FALSE){
	if(argsL$span <= 1){
		span<-as.numeric(argsL$span)
		lo<-loess(frame$V2 ~ frame$V1, span=span)
	}else{
		span<-as.numeric(argsL$span)/nrow(frame)
		lo<-loess(frame$V2 ~ frame$V1, span=span)
	}
}else{
	values<-seq(1,0.05,-0.01)
	span<-c()
	wss<-data.frame(matrix(ncol=3,nrow=0))
	for(i in values){
		if(i*nrow(frame) > 1){
		lo<-loess(frame$V2 ~ frame$V1, span=i)
		if(is.na(lo$s)==FALSE){
			out2<-predict(lo,index)
			if(sum(out2<0)==0){
				rough<-sd(diff(out2))
				new3<-data.frame(i,lo$s,rough)
				wss<-rbind(wss,new3)
			}
		}
		}
	}
	wss<-wss[!is.infinite(rowSums(wss)),]
	wss$loscale<-scale(wss$lo.s)																								#Distance to line def. by minimum lo.s and rough values
	wss$roughscale<-scale(wss$rough)																							#Distance to line def. by minimum lo.s and rough values
	dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-cbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}								#Distance to line def. by minimum lo.s and rough values
	a2<-c(wss$loscale[wss$loscale==min(wss$loscale)][1],wss$roughscale[wss$loscale==min(wss$loscale)][1])						#Distance to line def. by minimum lo.s and rough values
	b2<-c(wss$loscale[wss$roughscale==min(wss$roughscale)][1],wss$roughscale[wss$roughscale==min(wss$roughscale)][1])			#Distance to line def. by minimum lo.s and rough values
	wss$dist<-apply(wss,1,function(x) dist2d(c(x[4],x[5]),a2,b2))																#Distance to line def. by minimum lo.s and rough values
	span<-wss$i[wss$dist==max(wss$dist)][1]																						#Distance to line def. by minimum lo.s and rough values
	lo<-loess(frame$V2 ~ frame$V1, span=span)
}
out3<-predict(lo, index, se=T)
spanout<-data.frame(argsL$line,span,lo$s,out3$fit[index==0],out3$se[index==0])
write.table(spanout, file=paste(argsL$output, ".spans.txt", sep=""), append=TRUE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
out<-predict(lo, index)
new<-data.frame(index=index, value=out)
new2<-new$value[new$index >= -200 & new$index <= 200]
if(min(new$index) > -200){
	new3<-data.frame(index=seq(-200,min(new$index)-1,1), value=rep(500,times=length(seq(-200,min(new$index)-1,1))))
	new<-rbind(new3,new)
	new2<-new$value[new$index >= -200 & new$index <= 200]
}
if(max(new$index) < 200){
	new3<-data.frame(index=seq(max(new$index)+1,200,1), value=rep(500,times=length(seq(max(new$index)+1,200,1))))
	new<-rbind(new,new3)
	new2<-new$value[new$index >= -200 & new$index <= 200]
}
if(exists("new2") == FALSE){
	new2<-rep(NA, 401)
}
write.table(paste(new2, collapse="\t"), file="", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
