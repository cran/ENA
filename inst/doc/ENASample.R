### R code from vignette source 'ENASample.Rnw'

###################################################
### code chunk number 1: simulateNoise
###################################################
simul <- matrix(rnorm(50*20), ncol=50)
colnames(simul) <- paste("s", 1:50, sep="")
rownames(simul) <- paste("g", 1:20, sep="")


###################################################
### code chunk number 2: AddConnections
###################################################
#generate a few obvious "trends" which should establish correlations on a few genes.
simul[5,] <- seq(from=-1, to=1, length.out = ncol(simul))
simul[10,] <- seq(from=1, to=-1, length.out = ncol(simul))
simul[15,] <- seq(from=-1, to=1, length.out = ncol(simul))
simul[20,] <- seq(from=1, to=-1, length.out = ncol(simul))


###################################################
### code chunk number 3: bootstrapSpace
###################################################
library(ENA)
invisible(capture.output(sp <- bootstrap(simul, "buildSpace")))

head(sp)


###################################################
### code chunk number 4: ReconstructAll
###################################################
gn <- buildGenenet(simul)
wg <- buildWgcna(simul)
ar <- buildAracne(simul)


###################################################
### code chunk number 5: symmetricize
###################################################
gn <- symmetricize(abs(gn), "max", adjacencyList=TRUE)
colnames(gn)[3] <- "GeneNet"
wg <- symmetricize(abs(wg), "max", adjacencyList=TRUE)
colnames(wg)[3] <- "WGCNA"
ar <- symmetricize(abs(ar), "max", adjacencyList=TRUE)
colnames(ar)[3] <- "Aracne"


###################################################
### code chunk number 6: mergeMethods
###################################################
add <- getTableAddressing(rownames(simul))

#if all the addressing columns are identical, we can use the faster method.
if (identical(add[,1:2], gn[,1:2]) && identical(add[,1:2], wg[,1:2]) && identical(gn[,1:2], sp[,1:2]) && identical(gn[,1:2], ar[,1:2])){
	agg <- cbind(add, gn[,3], wg[,3], sp[,3], ar[,3])
} else{
	#build up using merge in case the order of the addressing is incorrect.
	agg <- merge(add, gn)	
	agg <- merge(agg, wg)	
	agg <- merge(agg, sp)	
	agg <- merge(agg, ar)	
}

head(agg)



###################################################
### code chunk number 7: ena
###################################################
ena <- ena(agg[,3:6])
agg$ENA <- ena


###################################################
### code chunk number 8: inspect
###################################################
cbind(agg[,1:2], round(agg[,3:7], 2))


