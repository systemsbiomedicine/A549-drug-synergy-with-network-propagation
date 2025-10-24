
############################################################
## Data preparation
############################################################

X <- read.csv('filesData/Network_Popagation_Result_alpha0.5(-log).csv', header=T, row.names=1)
# This is the -log10 of the network propagation result with alpha = 0.5
# Equivalent to -log10 of data in Supplementary Data S3.csv

W <- read.csv('filesData/W.csv', header=T, row.names=1)
colnames(W) <- paste0('Metagene ', 1:26)
rownames(W) <- rownames(X) 

gene.names <- row.names(X)

############################################################

gene.names[ order( W[, 'Metagene 2'], decreasing=F )[1:200] ]
write.table( data.frame( gene.names[ order( W[, 'Metagene 2'], decreasing=F )[1:200] ] ), file='200metagene2.csv', row.names=F, col.names=F )
