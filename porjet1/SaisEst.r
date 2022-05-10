SaisEst <- function(x,freq){
	SaisEstfun <- function(i,x,freq){
		n <- length(x)
		ind <- seq(i,n,freq)
		mean(x[ind],na.rm=TRUE)
	}
	motif <- sapply(1:freq,SaisEstfun,x,freq)
	list(motif=motif,serie=ts(rep(motif,length.out=length(x)),start=start(x),frequency=frequency(x)))
}
