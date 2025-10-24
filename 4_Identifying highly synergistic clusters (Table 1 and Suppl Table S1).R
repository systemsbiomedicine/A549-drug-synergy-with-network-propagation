
############################################################
## Data preparation
############################################################

dc.info <- read.csv('filesData/drug_combination_infomation.csv', header=T, row.names=NULL)
# Equivalent to Supplementary Data S1.csv

X <- read.csv('filesData/Network_Popagation_Result_alpha0.5(-log).csv', header=T, row.names=1)
# This is the -log10 of the network propagation result with alpha = 0.5
# Equivalent to -log10 of data in Supplementary Data S3.csv

H <- read.csv('filesData/H.csv', header=T, row.names=1)
colnames(H) <- colnames(X)
rownames(H) <- paste0('Metagene', 1:26)

############################################################
## Hierarchical clustering
############################################################

hc <- hclust( dist(t(H), method = 'euclidean'), method='ward.D2')

############################################################
## As the 'height' of cutree varies, identify clusters whose average synergy scores ≥ Q3
############################################################

identify.clusters <- function(synergy_score)
{
	quantile( dc.info[, synergy_score], probs = c(0.25, 0.50, 0.75)) -> Q
	df <- data.frame()
	my.block.id <- c()
	last.ct <- 0
	for( h in rev(hc$height) )
	{
		cutree( hc, h=h ) -> ct
		dc.info$new.group <- ct
		
		if( prod( ct == last.ct ) == 1 ) 
		{ 	
			# The cut gives the same partitions as the previous cut.
			next
		}
		
		for( i in 1:max(ct) )
		{
			dc.info[ dc.info$new.group == i, ] -> sub.group.info		
			length( sub.group.info[ , synergy_score] ) -> N
			if( N < 10 ) { next }
			mean( sub.group.info[ , synergy_score] ) -> mean1
			if(mean1 > Q[3]) 
			{
				paste( sub.group.info$block_id, collapse=',') -> new.block.id
				if( length(my.block.id) > 0 )
				{
					if( new.block.id %in% my.block.id )
					next;
				}
				
				sd( sub.group.info[ , synergy_score] ) -> sd1
				min( sub.group.info[ , synergy_score] ) -> min1
				max( sub.group.info[ , synergy_score] ) -> max1
				table( c( sub.group.info[ , 'drug_row'], sub.group.info[ , 'drug_col'] ) ) -> tab
				names(sort(tab, decreasing=T))[1] -> main.drug
				as.numeric( sort(tab, decreasing=T)[1] ) -> main.drug.N
				names(sort(tab, decreasing=T))[2] -> main.drug.2
				as.numeric( sort(tab, decreasing=T)[2] ) -> main.drug.2.N
				
				df <- rbind( df,
								c( synergy_score, h, i, N, mean1, sd1, min1, max1, main.drug, main.drug.N, main.drug.2, main.drug.2.N)	
							)
				my.block.id <- c(my.block.id, new.block.id)				
			}
		}
		ct -> last.ct
	}
	colnames( df ) <- c('synergy_score', 'h', 'cluster_id', 'N', 'mean', 'sd', 'min', 'max', 'main.drug', 'main.drug.N', 'main.drug.2', 'main.drug.2.N')
	df$cluster_id <- as.integer(df$cluster_id)
	df$h <- as.numeric(df$h)
	df$N <- as.integer(df$N)
	df$mean <- as.numeric(df$mean)
	df$sd <- as.numeric(df$sd)
	df$min <- as.numeric(df$min)
	df$max <- as.numeric(df$max)
	df$main.drug.N <- as.numeric(df$main.drug.N)
	df$main.drug.2.N <- as.numeric(df$main.drug.2.N)
	
	ceiling( df$h*10^2 )/10^2 -> df$h  ## round up (ceiling) at the digits = 2 to avoid the rounding down problem when cutting the dendrogram
	round( df$mean, digits=2 ) -> df$mean
	round( df$sd, digits=2 ) -> df$sd
	
	return(df)
}

identify.clusters('synergy_hsa') -> df1
identify.clusters('synergy_zip') -> df2
identify.clusters('synergy_bliss') -> df3
identify.clusters('synergy_loewe') -> df4

library(dplyr)
df.list <- list(df1, df2, df3, df4)
key.column <- 'h'
df.combined <- Reduce( function(x, y) full_join(x, y, by = key.column), df.list )
df.combined[ order( -df.combined$h ), ] -> df.combined

df.combined[, c(
			'h',
			'synergy_score.x','synergy_score.y','synergy_score.x.x','synergy_score.y.y',
			'cluster_id.x','cluster_id.y','cluster_id.x.x','cluster_id.y.y',
			'N.x','N.y','N.x.x','N.y.y',
			'mean.x','mean.y','mean.x.x','mean.y.y',
			'sd.x','sd.y','sd.x.x','sd.y.y',
			'main.drug.x','main.drug.y','main.drug.x.x','main.drug.y.y',
			'main.drug.N.x','main.drug.N.y','main.drug.N.x.x','main.drug.N.y.y'
)] -> df.combined

colnames(df.combined) <- c(
			'h',
			'synergy_score.hsa','synergy_score.zip','synergy_score.bliss','synergy_score.loewe',
			'cluster_id.hsa','cluster_id.zip','cluster_id.bliss','cluster_id.loewe',
			'N.hsa','N.zip','N.bliss','N.loewe',
			'mean.hsa','mean.zip','mean.bliss','mean.loewe',
			'sd.hsa','sd.zip','sd.bliss','sd.loewe',
			'main.drug.hsa','main.drug.zip','main.drug.bliss','main.drug.loewe',
			'main.drug.N.hsa','main.drug.N.zip','main.drug.N.bliss','main.drug.N.loewe'
)


write.csv( df.combined, file='highly_synergistic_clusters_results.csv', row.names=F)
