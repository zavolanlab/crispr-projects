args = commandArgs(trailingOnly=TRUE)

# test if there are all the necessary arguments
if (length(args)!=3) {
  stop("3 arguments must be supplied: counts_file_path, metadata_file_path, outdir_path", call.=FALSE)
}

require(dplyr)
require(stringr)
require(data.table)
library(BiocParallel)
library("DEXSeq")

count_file_path = args[1]
metadata_file_path = args[2]
outdir_path = args[3]

# helper function to parse input paths
trim_quotes_from_path<-function(path_string){
path_string = gsub("'", '', path_string)
path_string = gsub('"', '', path_string)
return(path_string)
}

# loading count matrix
cts <- read.csv(trim_quotes_from_path(count_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)

# loading metadata
sampleData <- read.csv(trim_quotes_from_path(metadata_file_path), sep='\t',header = TRUE,stringsAsFactors=FALSE,check.names = FALSE)
colnames(sampleData) <- c('sample','condition')
rownames(sampleData) <- sampleData$sample

groupID <- cts$dexseq_groupID
featureID <- cts$dexseq_featureID

countData <- data.matrix(cts[ ,rownames(sampleData)])

design <- formula( ~ sample + exon + condition:exon )

BPPARAM = MulticoreParam(workers=1)

dxd = DEXSeqDataSet( countData, sampleData, design, featureID, groupID )
dxd = estimateSizeFactors( dxd )

colData_df <- data.frame(colData(dxd))
colData_df$colname <- paste(colData_df$sample,':',colData_df$exon)

normalized_counts <- data.frame(counts(dxd, normalized=TRUE))
colnames(normalized_counts) <- (colData_df$colname)

write.table(normalized_counts,file=paste(trim_quotes_from_path(outdir_path),"/DEXseq_norm_counts.tsv",sep = ""),row.names = TRUE, quote=FALSE, sep='\t')