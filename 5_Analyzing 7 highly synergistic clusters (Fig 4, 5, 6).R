
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
rownames(H) <- paste0('Metagene ', 1:26)

############################################################
## Hierarchical clustering
############################################################

hc <- hclust( dist(t(H), method = 'euclidean'), method='ward.D2')

############################################################
## Functions
############################################################

my.cut <- function( h, group )
{
	cutree( hc, h=h ) -> ct
	dc.info$new.group <- ct
	return( dc.info[dc.info$new.group==group,] )
}

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

# offset the h value by 1e-5 to avoid the rounding problem

drug1 <- 'Crizotinib'
h <- 148.47063 + 1e-5
g <- 26
my.cut(h, g) |> rownames() |> as.numeric() -> index1

drug2 <- 'Dasatinib'
h <- 115.6835 + 1e-5
g <- 28
my.cut(h, g) |> rownames() |> as.numeric() -> index2

drug3 <- 'Erlotinib'
h <- 104.1198 + 1e-5
g <- 27
my.cut(h, g) |> rownames() |> as.numeric() -> index3

drug4 <- 'Gefitinib'
h <- 102.69 + 1e-5
g <- 34
my.cut(h, g) |> rownames() |> as.numeric() -> index4

drug5 <- 'Mitotane'
h <- 66.85398 + 1e-5
g <- 16
my.cut(h, g) |> rownames() |> as.numeric() -> index5

drug6 <- 'Paclitaxel'
h <- 133.2709 + 1e-5
g <- 19
my.cut(h, g) |> rownames() |> as.numeric() -> index6

drug7 <- 'Pazopanib'
h <- 155.4543 + 1e-5
g <- 22
my.cut(h, g) |> rownames() |> as.numeric() -> index7

drug8 <- 'Quinacrine'
h <- 124.1335 + 1e-5
g <- 2
my.cut(h, g) |> rownames() |> as.numeric() -> index8

drug9 <- 'Ruxolitinib'
h <- 78.18202 + 1e-5
g <- 69
my.cut(h, g) |> rownames() |> as.numeric() -> index9

drug10 <- 'Vemurafenib'
h <- 190.4986 + 1e-5
g <- 21
my.cut(h, g) |> rownames() |> as.numeric() -> index10

############################################################
 
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
 
############################################################

t(H)[ 	c(index1, 
		index2,
		index3,
		index4,
		index5,
		index6,
		index7,
		index8,
		index9,
		index10
		), ] -> df.wide

pdf(file='Fig4_boxplot_all.pdf', width=6, height=6)
par(mar = c(10, 4, 4, 2) + 0.1)
boxplot( df.wide, las=2, ylab='Metagene value' )
# Reset margins to default after plotting
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()

## Next, we show that metagene values are NOT normally distributed
for( i in 1:26 )
{
	dev.new(width=4, height=4)
	qqnorm( df.wide[, i], main=i )
	qqline( df.wide[, i] )
}

my.p.values <- c()
index <- 1:26
for( i in index[-2] )
{
	my.p.values <- c( 	my.p.values, 	
						wilcox.test( df.wide[,2], df.wide[,i], paried=TRUE, alternative='greater' )$p.value 
					)
}
p.adjust(my.p.values) < 0.05

############################################################

library(pheatmap)

## Strong Metagene 2 Panel A
swap.drug.names('Dasatinib', index2)
swap.drug.names('Paclitaxel', index6)
swap.drug.names('Quinacrine', index8)
swap.drug.names('Erlotinib', index4)
swap.drug.names('Gefitinib', index4)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
c(index2, index6, index8, index4) -> index

pdf(file='Fig5a_heatmap_metagene2_panelA.pdf')
pheatmap( t(H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index2), length(index2)+length(index6), length(index2)+length(index6)+length(index8)) )
dev.off()

## Strong Metagene 2 Panel B
swap.drug.names('Crizotinib', index1)
swap.drug.names('Pazotinib', index7)
swap.drug.names('Ruxolitinib', index9)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
c(index1, index7, index9) -> index

pdf(file='Fig5b_heatmap_metagene2_panelB.pdf')
pheatmap( t(H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index1), length(index1)+length(index7)) )
dev.off()

## The rest with no Metagene 2 signal
swap.drug.names('Erlotinib', index3)
swap.drug.names('Gefitinib', index3)
swap.drug.names('Mitotane', index5)
swap.drug.names('Vemurafenib', index10)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
c(index3, index5, index10) -> index

pdf(file='heatmap_no_metagene2.pdf')
pheatmap( t(H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c( length(index3), length(index3)+length(index5) ) )
dev.off()

############################################################

############################################################
# All Dasatinib pairs cut into 2 groups using matrix H

my.drug <- 'Dasatinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='less' ) -> test.result
format( test.result$p.value, digits=2 ) -> print.p.value.1

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='less' ) -> test.result
format( test.result$p.value, digits=2 ) -> print.p.value.2

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='less' ) -> test.result
format( test.result$p.value, digits=2 ) -> print.p.value.3

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='less' ) -> test.result
format( test.result$p.value, digits=2 ) -> print.p.value.4

draw_bracket <- function(i, j, y_pos, label, text.offset) {
  # Draw the horizontal line
  segments(i, y_pos, j, y_pos, lwd = 1.5)
  # Draw the vertical ticks
  segments(i, y_pos, i, y_pos - (y_pos * 0.02), lwd = 1.5)
  segments(j, y_pos, j, y_pos - (y_pos * 0.02), lwd = 1.5)
  # Add the label
  text((i + j) / 2, y_pos + (y_pos * text.offset), label, cex = 0.9)
}

pdf(file='Fig6b_boxplot_dasatinib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))

### Sub-plot 1
y_max <- max(sub.dc.info$synergy_hsa, na.rm = TRUE)
y_limit <- y_max * 1.25
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA', ylim = c(min(sub.dc.info$synergy_hsa), y_limit) )
draw_bracket(	1, 2, 
				y_pos = y_max + 1, 
				label = paste0('p-value = ', print.p.value.1 ),
				text.offset = 0.05
			)

### Sub-plot 2
y_max <- max(sub.dc.info$synergy_zip, na.rm = TRUE)
y_limit <- y_max * 1.25
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP', ylim = c(min(sub.dc.info$synergy_zip), y_limit) )
draw_bracket(	1, 2, 
				y_pos = y_max + 1, 
				label = paste0('p-value = ', print.p.value.2 ),
				text.offset = 0.05
			)

### Sub-plot 3
y_max <- max(sub.dc.info$synergy_bliss, na.rm = TRUE)
y_limit <- y_max * 1.25
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS', ylim = c(min(sub.dc.info$synergy_bliss), y_limit) )
draw_bracket(	1, 2, 
				y_pos = y_max + 1, 
				label = paste0('p-value = ', print.p.value.3 ),
				text.offset = 0.05
			)

### Sub-plot 4
y_max <- max(sub.dc.info$synergy_loewe, na.rm = TRUE)
y_limit <- y_max * 1.75
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE', ylim = c(min(sub.dc.info$synergy_loewe), y_limit) )
draw_bracket(	1, 2, 
				y_pos = y_max + 1, 
				label = paste0('p-value = ', print.p.value.4 ),
				text.offset = 0.25
			)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index

pdf(file='Fig6a_heatmap_dasatinib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

############################################################

############################################################
# All Paclitaxel pairs cut into 2 groups using matrix H

my.drug <- 'Paclitaxel'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Paclitaxel.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Paclitaxel.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='greater' )

############################################################

############################################################
# All Quinacrine pairs cut into 2 groups using matrix H

my.drug <- 'Quinacrine'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Quinacrine.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Quinacrine.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='less' )

############################################################

############################################################
# All Pazopanib pairs cut into 2 groups using matrix H

my.drug <- 'Pazopanib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Pazopanib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Pazopanib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='greater' )

############################################################

############################################################
# All Crizotinib pairs cut into 2 groups using matrix H

my.drug <- 'Crizotinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Crizotinib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Crizotinib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='greater' )

############################################################

############################################################
# All Erlotinib pairs cut into 2 groups using matrix H

my.drug <- 'Erlotinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Erlotinib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Erlotinib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='greater' )

############################################################

############################################################
# All Gefitinib pairs cut into 2 groups using matrix H

my.drug <- 'Gefitinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Gefitinib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Gefitinib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='less' )

############################################################

############################################################
# All Ruxolitinib pairs cut into 2 groups using matrix H

my.drug <- 'Ruxolitinib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Ruxolitinib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
c(index.g1, index.g2) -> index
pdf(file='heatmap_Ruxolitinib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2)) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='greater' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='greater' )

############################################################

############################################################
# All Mitotane pairs cut into 2 groups using matrix H

my.drug <- 'Mitotane'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Mitotane.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
which( sub.dc.info$new.group == 3 ) -> index.g3
c(index.g1, index.g2, index.g3) -> index
pdf(file='heatmap_Mitotane.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2), length(index.g1)+length(index.g2)+length(index.g3) ) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='less' )

############################################################

############################################################
# All Vemurafenib pairs cut into 2 groups using matrix H

my.drug <- 'Vemurafenib'
c( 	which( dc.info$drug_row == my.drug ), which( dc.info$drug_col == my.drug ) ) -> index
H[,index] -> sub.H
dc.info[index,] -> sub.dc.info

swap.drug.names(my.drug, 1:607)
paste( dc.info$drug_row, '+', dc.info$drug_col ) -> dc.names
dc.names[index] -> sub.dc.names

sub.hc <- hclust( dist(t(sub.H), method = 'euclidean'), method='ward.D2')
cutree( sub.hc, k=2 ) -> ct
sub.dc.info$new.group <- ct

pdf(file='boxplot_Vemurafenib.pdf')
par(mfrow = c(2, 2), mar = c(3, 4.1, 3, 2.1))
boxplot( synergy_hsa ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='HSA' )
boxplot( synergy_zip ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='ZIP' )
boxplot( synergy_bliss ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='BLISS' )
boxplot( synergy_loewe ~ new.group, data=sub.dc.info, names=c('Group 1','Group 2'), xlab=NA, ylab='LOEWE' )
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Always reset plotting parameters after you're done
dev.off()

which( sub.dc.info$new.group == 1 ) -> index.g1
which( sub.dc.info$new.group == 2 ) -> index.g2
which( sub.dc.info$new.group == 3 ) -> index.g3
which( sub.dc.info$new.group == 4 ) -> index.g4
c(index.g1, index.g2, index.g3, index.g4) -> index
pdf(file='heatmap_Vemurafenib.pdf')
pheatmap( t(sub.H)[index,], cellwidth=8, cellheight=8, scale='none', cluster_rows=F, cluster_cols=F, legend=F, color=rev(heat.colors(10)), labels_row=sub.dc.names[index], fontsize_col=8, fontsize_row=8,
gaps_row = c(length(index.g1), length(index.g1)+length(index.g2), length(index.g1)+length(index.g2)+length(index.g3), length(index.g1)+length(index.g2)+length(index.g3)+length(index.g4) ) )
dev.off()

aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_hsa ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_hsa ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_zip ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_zip ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_bliss ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_bliss ~ new.group, data=sub.dc.info, alternative='less' )

aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=mean) |> round(digits=2)
aggregate(synergy_loewe ~ new.group, data=sub.dc.info, FUN=sd) |> round(digits=2)
t.test( synergy_loewe ~ new.group, data=sub.dc.info, alternative='less' )

############################################################

############################################################
# Analyze drug partners of the seven highly synergistic clusters with a strong Metagene 2 signal.

dc.info[ c(		index1, 
				index2,
				index4,
				index6,
				index7,
				index8,
				index9
				), c('drug_row','drug_col')] -> dc.info.7

c( dc.info.7$drug_row, dc.info.7$drug_col ) |> table() |> sort()

## highly frequent partners include: Lapatinib, Vemurafenib, Sunitinib, Axitinib, Vismodegib, Sorafenib, Lenalidomide, Vandetanib, and Imatinib

############################################################

