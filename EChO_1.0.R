#!usr/bin/Rscript

rm(list=ls())
## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 3) {
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
			--output=someValue   - Output prefix
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
if(is.null(argsL$frame) | is.null(argsL$center) | is.null(argsL$chr)) {
  stop("Argument is missing!
     Calculate area under the curve threshold for CUT&RUN peaks 
 
      Arguments:
			--frame=someValue   - Input offset vs. fragment size data frame
			--center=someValue   - Peak center
			--chr=someValue   - Peak chromosome
			--output=someValue   - Output prefix
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
			--span=someValue   - span (auto-select with GCV if non-numeric entry)
")
 
  q(save="no")
}

frame<-read.table(argsL$frame)
center<-as.vector(argsL$center)
chr<-as.vector(argsL$chr)
frame<-frame[order(frame$V1),]
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
	wss$loscale<-scale(wss$lo.s)																																													#Distance to line def. by minimum lo.s and rough values
	wss$roughscale<-scale(wss$rough)																																											#Distance to line def. by minimum lo.s and rough values
	dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-cbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}													#Distance to line def. by minimum lo.s and rough values
	a2<-c(wss$loscale[wss$loscale==min(wss$loscale)][1],wss$roughscale[wss$loscale==min(wss$loscale)][1])								#Distance to line def. by minimum lo.s and rough values
	b2<-c(wss$loscale[wss$roughscale==min(wss$roughscale)][1],wss$roughscale[wss$roughscale==min(wss$roughscale)][1])			#Distance to line def. by minimum lo.s and rough values
	wss$dist<-apply(wss,1,function(x) dist2d(c(x[4],x[5]),a2,b2))																													#Distance to line def. by minimum lo.s and rough values
	span<-wss$i[wss$dist==max(wss$dist)][1]																																								#Distance to line def. by minimum lo.s and rough values
	lo<-loess(frame$V2 ~ frame$V1, span=span)
}
out<-predict(lo, index, se=T)
infl<-c(FALSE, diff(diff(out$fit)>0)==1, FALSE)
new<-data.frame(index=index, value=out$fit, disp=out$se/out$fit, TF=infl)
offsets<-new$index[new$TF==TRUE]
lengths<-new$value[new$TF==TRUE]
disps<-new$disp[new$TF==TRUE]
summits<-offsets+as.numeric(center[1])
new2<-data.frame(chrom=rep(chr[1], each=length(summits)), start=summits, stop=summits+1, lengths=lengths, disps=disps, span=rep(span, each=length(summits)))
new2<-new2[new2$lengths > 0,]
if(length(new2) > 0){
	write.table(new2, file="", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
