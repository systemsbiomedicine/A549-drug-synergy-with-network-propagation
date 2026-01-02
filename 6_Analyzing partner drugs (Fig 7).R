
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
## Function
############################################################

swap.drug.names <- function(my.drug, index)
{
	for( i in index )
	{
		if( dc.info[i, 'drug_col'] == my.drug )
		{
			dc.info[ i , 'drug_col' ] <<- dc.info[ i, 'drug_row' ]
			dc.info[ i, 'drug_row' ] <<- my.drug		
		}		
	}		
}

############################################################

partner.drugs <- c( 'Lapatinib', 'Vemurafenib', 'Sunitinib', 'Axitinib', 'Vismodegib', 'Sorafenib', 'Lenalidomide', 'Vandetanib', 'Imatinib' )

df.new <- data.frame()
for( i in partner.drugs )
{
	c( which(dc.info$drug_row == i), which(dc.info$drug_col == i) ) -> index
	dc.info[ index, ] -> sub.dc.info
	sub.dc.info$class <- i
	df.new <- rbind( df.new, sub.dc.info )
}

# Check if the synergy scores are signifantly higher than 0
my.p.values <- c()
names <- c()
for( drug in unique(df.new$class) )
{
	for( synergy_score in c('synergy_hsa', 'synergy_zip', 'synergy_bliss', 'synergy_loewe') )
	{
		df.new[ df.new$class == drug, synergy_score] -> dat	
		dev.new(width=3, height=3)
		qqnorm( dat )
		qqline( dat )
		my.p.values <- c(my.p.values, wilcox.test( dat, mu=0, alternative='greater' )$p.value)
		names <- c(names, paste0( drug, ':', synergy_score))
	}
}
p.adjust(my.p.values, method='fdr')
names

# next convert to a long table format
df.long <- data.frame()
for( i in 1:nrow(df.new) )
{
	df.long <- rbind( df.long,
						rbind(
							c( df.new[i, 'class'], 'HSA', df.new[i, 'synergy_hsa'] ),
							c( df.new[i, 'class'], 'ZIP', df.new[i, 'synergy_zip'] ),
							c( df.new[i, 'class'], 'BLISS', df.new[i, 'synergy_bliss'] ),
							c( df.new[i, 'class'], 'LOEWE', df.new[i, 'synergy_loewe'] )
							)
						)
}
colnames(df.long) <- c('drug','synergy_type','synergy_score')
df.long$drug <- as.factor(df.long$drug)
df.long$synergy_type <- as.factor(df.long$synergy_type)
df.long$synergy_score <- as.numeric(df.long$synergy_score)

pdf(file='Fig7a_partner_drugs.pdf', width=8, height=6)
par(mar = c(9, 4, 4, 2) + 0.1)
bp <- boxplot( synergy_score ~ synergy_type * drug, data=df.long, plot=F )
boxplot( synergy_score ~ synergy_type * drug, data=df.long, xlab=NA, ylab='Synergy score',
col = rep(c("lightblue", "lightgreen", "lightcoral", "gold") ), 
xaxt='n',
at=c(1:4, 6:9, 11:14, 16:19, 21:24, 26:29, 31:34, 36:39, 41:44),
ylim=c( min(bp$out)-1, max(bp$out)+5 )
 )

axis( side=1, las=2,
at=c(2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5),
labels=c( 'Axitinib', 'Imatinib', 'Lapatinib', 'Lenalidomide', 'Sorafenib', 'Sunitinib', 'Vandetanib', 'Vemurafenib',  'Vismodegib' ), cex.axis=1 )

legend("bottomright",
       legend = c('BLISS','HSA','LOEWE','ZIP'),
       fill = c("lightblue", "lightgreen", "lightcoral", "gold"),    
       cex = 0.8 # Font size control
)

text( 36, max(bp$out)+4, '*', cex=2 )

# Reset margins to default after plotting (important!)
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()

############################################################
my.drug <- 'Lapatinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names

library(pheatmap)
pdf(file='Fig7b_heatmap_Lapatinib.pdf')
pheatmap( t(H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=dc.names[index], fontsize_col=8, fontsize_row=8)
dev.off()
############################################################
############################################################
my.drug <- 'Imatinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names

library(pheatmap)
pdf(file='Fig7c_heatmap_Imatinib.pdf')
pheatmap( t(H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=dc.names[index], fontsize_col=8, fontsize_row=8)
dev.off()
############################################################

