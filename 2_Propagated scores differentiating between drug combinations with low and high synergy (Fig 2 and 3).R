
############################################################
## Data preparation
############################################################

dc.info <- read.csv('filesData/drug_combination_infomation.csv', header=T, row.names=NULL)
# Equivalent to Supplementary Data S1.csv

X <- read.csv('filesData/Network_Popagation_Result_alpha0.5(-log).csv', header=T, row.names=1)
# This is the -log10 of the network propagation result with alpha = 0.5
# Equivalent to -log10 of data in Supplementary Data S3.csv

###############################################################################
## Can the propagated score of each gene effectively distinguish 
## between Quarter 1 and Quarter 4 synergy scores?
###############################################################################

get.df.FC.pvalue <- function(synergy_score)
{
	quantile( dc.info[, synergy_score], probs = c(0.25, 0.50, 0.75)) ->> Q
	Q.score <<- c()
	for( i in 1:nrow(dc.info))
	{
		if( dc.info[i, synergy_score] <= Q[1] )
		{
			Q.score <<- c(Q.score, 1)			
		} else if( dc.info[i, synergy_score] <= Q[2] )
		{
			Q.score <<- c(Q.score, 2)	
		} else if( dc.info[i, synergy_score] <= Q[3] )
		{
			Q.score <<- c(Q.score, 3)
		} else
		{
			Q.score <<- c(Q.score, 4)
		}
	}
	
	c( which( Q.score == 1 ), which( Q.score == 4 ) ) -> index
	Q.score[index] -> Q.score.sub
	X[, index] -> X.sub
	
	df.FC <- data.frame()
	for( i in 1:nrow(X.sub) )
	{
		# For Wilcoxson tests, using propagated scores, log10(propagated scores), -propagated scores, -log10(propagated scores) yield the same results
		wilcox.test( as.numeric(X.sub[i, ]) ~ Q.score.sub, alternative='two.sided' ) -> test2
		aggregate( 10^(-as.numeric(X.sub[i, ])) ~ Q.score.sub, FUN = median) -> df.aggregate
		df.aggregate[2,2] / df.aggregate[1,2] -> fold.change
		df.FC <- rbind( df.FC, c(fold.change, test2$p.value) )
	}
	colnames(df.FC) <- c('fold.change','p.value')
	
	df.FC$p.adjusted.value <- p.adjust( df.FC$p.value, method='fdr' ) 
	df.FC[ order( df.FC$p.adjusted.value ), ] -> df.FC
	return( df.FC )
}

## Boxplot: Propagated scores: Quarter 1-4
pdf(file='Fig2_boxplot.pdf', width=10, height=8)
par(mfrow = c(2, 2))
get.df.FC.pvalue('synergy_hsa') -> df.HSA
as.numeric( rownames(df.HSA)[1] ) -> index
boxplot( 10^(-as.numeric( X[index,])) ~ paste0('Quarter ', Q.score), xlab='HSA', ylab='Propagated score', main=rownames(X)[index] )

get.df.FC.pvalue('synergy_zip') -> df.ZIP
as.numeric( rownames(df.ZIP)[1] ) -> index
boxplot( 10^(-as.numeric( X[index,])) ~ paste0('Quarter ', Q.score), xlab='ZIP', ylab='Propagated score', main=rownames(X)[index] )

get.df.FC.pvalue('synergy_bliss') -> df.BLISS
as.numeric( rownames(df.BLISS)[1] ) -> index
boxplot( 10^(-as.numeric( X[index,])) ~ paste0('Quarter ', Q.score), xlab='BLISS', ylab='Propagated score', main=rownames(X)[index] )

get.df.FC.pvalue('synergy_loewe') -> df.LOEWE
as.numeric( rownames(df.LOEWE)[1] ) -> index
boxplot( 10^(-as.numeric( X[index,])) ~ paste0('Quarter ', Q.score), xlab='LOEWE', ylab='Propagated score', main=rownames(X)[index] )

dev.off()


## Scatter plot: Propagated scores versus Synergy scores
pdf(file='Fig3_scatter.pdf', width=8, height=8)
par(mfrow = c(2, 2))
get.df.FC.pvalue('synergy_hsa') -> df.HSA
as.numeric( rownames(df.HSA)[1] ) -> index
plot( 10^(-as.numeric( X[index,])), dc.info[,'synergy_hsa'], col=Q.score, pch=16, main=rownames(X)[index], xlab='Propagated score', ylab='HSA' )

get.df.FC.pvalue('synergy_zip') -> df.ZIP
as.numeric( rownames(df.ZIP)[1] ) -> index
plot( 10^(-as.numeric( X[index,])), dc.info[,'synergy_zip'], col=Q.score, pch=16, main=rownames(X)[index], xlab='Propagated score', ylab='ZIP' )

get.df.FC.pvalue('synergy_bliss') -> df.BLISS
as.numeric( rownames(df.BLISS)[1] ) -> index
plot( 10^(-as.numeric( X[index,])), dc.info[,'synergy_bliss'], col=Q.score, pch=16, main=rownames(X)[index], xlab='Propagated score', ylab='BLISS' )

get.df.FC.pvalue('synergy_loewe') -> df.LOEWE
as.numeric( rownames(df.LOEWE)[1] ) -> index
plot( 10^(-as.numeric( X[index,])), dc.info[,'synergy_loewe'], col=Q.score, pch=16, main=rownames(X)[index], xlab='Propagated score', ylab='LOEWE' )

dev.off()


cutoff <- 1E-4
length( which( df.HSA$p.adjusted.value <= cutoff ) )
length( which( df.ZIP$p.adjusted.value <= cutoff ) )
length( which( df.BLISS$p.adjusted.value <= cutoff ) )
length( which( df.LOEWE$p.adjusted.value <= cutoff ) )


