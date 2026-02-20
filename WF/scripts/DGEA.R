# DIFFERENTIAL GENE EXPRESSION ANALYSIS (DESeq2) #

#### DEPENDENCIES ####
library(Ckmeans.1d.dp, quietly = T)
library(RColorBrewer, quietly = T)
library(vroom, quietly = T)
library(tidyverse, quietly = T)
library(edgeR, quietly = T)
library(DESeq2, quietly = T)
library(ggplot2, quietly = T)
library(ggbeeswarm, quietly = T)

#### INPUT ####
args <- commandArgs(trailingOnly = T)

outdir_path <- args[1]
expr_matrix_path <- args[2]
EXP_samples <- args[3]
CTL_samples <- args[4]

#### INPUT DATA HANDLING ####
expr_matrix <- read.delim(expr_matrix_path, header = T, sep = "\t")

## MEGADEPTH COUNTS FILES READING ##
# creation of temp list of runs for each experiment
for (i in sample_data$experiment_id) {
  j <- as.list(unlist(as.list(list_df_withGTF[i])))
  names(j) <- sample_data$sample[sample_data$experiment_id==i]
  assign(i, list(j))
}

# aggregation of those lists into one single list (one for GTF and another for noGTF (later))
final_list_withGTF <- list()
for (i in ls(pattern='SRP')) {
  final_list_withGTF <- append(final_list_withGTF, get(i))
}
names(final_list_withGTF) <- unique(sample_data$experiment_id)

# rm of temp lists
for (i in ls(pattern='SRP')) {
  rm(i)
}

#### GTF LOOP ####
setwd(outdir)

# loop in final_list_withGTF for merge - norm - plot
# with gtf
for (i in final_list_withGTF) {
  cat('Analysis of ', unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  
  subDir <- as.character(unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  dir.create(file.path(outdir, subDir), showWarnings = FALSE)
  dir.create(file.path(subDir, 'withGTF'), showWarnings = FALSE)
  
  design <- data.frame('sample'= names(i), 'condition'= as.factor(sample_data$EXP_CTL[sample_data$sample==names(i)]))

  merged_df <- as.data.frame(map_dfc(i,
                                     ~vroom(.x,
                                            col_select=c(4),
                                            delim = '\t',
                                            col_names = FALSE,
                                            show_col_types = FALSE)))
  #renaming cols
  suppressWarnings(colnames(merged_df) <- names(i))
  suppressWarnings(merged_df <- head(merged_df, -2)) #remove last two lines
  rownames(merged_df) <- genenames$gene_id #gene ids as rownames
  
  # clustering for target gene norm counts #
  cpm_counts <- cpm(merged_df)
  target_counts <- cpm_counts[rownames(cpm_counts)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]
  
  k <- 2
  result <- Ckmeans.1d.dp(target_counts, k)
  colors <- brewer.pal(k, "Dark2")
  midpoints <- ahist(result, style="midpoints", data=target_counts, plot=FALSE)$breaks[2:k]
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'withGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'clustering', 'withGTF.pdf', sep = '_')), height = 2.5)
  plot(result,
       col.clusters = colors,
       main = paste('K-means clustering of', 
                    unique(unique(sample_data$experiment_id[sample_data$sample==names(i)]))),
       ylab = '',
       yaxt = 'n',
       xlab = paste(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), 'CPM', '(GTF)'),
       ylim=c(0,1))
  abline(v=midpoints, col="RoyalBlue", lwd=3, lty=3,)
  text(x = target_counts, y = par("usr")[3]-0.15, labels = paste(names(i), design$EXP_CTL[design$sample==names(i)]), xpd = NA, srt=45, cex=0.4, adj = 0.9)
  dev.off()
  
  ## DESEQ2 ##
  dds <- DESeqDataSetFromMatrix(countData = merged_df, 
                                colData = design, 
                                design= ~ condition)
  # FILTER 
  nrow(dds)
  keep <- rowSums(counts(dds) >= 1) >= length(design$sample)/2
  dds <- dds[keep,]
  nrow(dds)
  
  # VST COUNT TRANSFORMATION n PLOT
  vsd <- vst(dds, blind = FALSE)
  
  # PCA 
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'withGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'PCA', 'withGTF.pdf', sep = '_')))
  print(plotPCA(vsd, intgroup = c('condition')))
  dev.off()
  
  # DEA
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  summary(res)
  
  write.csv(as.data.frame(res), file = paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'withGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'withGTF', 'DE.tsv', sep = '_')), quote = F, sep = '\t')
  
  ## TARGET GENE COUNTS PLOT ##
  targetGene <- rownames(res)[rownames(res)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]
  geneCounts <- plotCounts(dds, 
                           gene = targetGene, 
                           intgroup = c('condition'),
                           transform = F,
                           returnData = TRUE)
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'withGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), targetGene, 'withGTF.pdf', sep = '_')))
  print(ggplot(geneCounts, aes(x = condition, y = count, color = condition, size=4)) + 
    theme_bw() +
    ggtitle(label = paste(targetGene, 'raw counts')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_beeswarm(cex = 3) +
    annotate(geom="text", 
             x=Inf, 
             y=res[rownames(res)==targetGene,]$baseMean, 
             label=paste('log2FC:', 
                         round(res[rownames(res)==targetGene,]$log2FoldChange, digits=3), '\n',
                         'pvalue:', 
                         round(res[rownames(res)==targetGene,]$pvalue, digits=3), '\n',
                         'padj:', 
                         round(res[rownames(res)==targetGene,]$padj, digits=3), '\n'),
             hjust=1, vjust=1,
             color="black") + 
    guides(condition = "legend", size = "none"))
  dev.off()
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'withGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'MAplot', 'withGTF.pdf', sep = '_')))
  plotMA(res,ylim = c(-max(res$log2FoldChange), max(res$log2FoldChange)))
  with(res[targetGene, ], {
    points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
    text(baseMean, log2FoldChange, targetGene, pos=2, col="red")
  })
  dev.off()
}

###############
#### NOGTF ####

# add gene count files paths_withGTF as column
sample_data$paths_noGTF <- file.path2(sample_dir, sample_data$LibraryLayout, sample_data$sample, path_to_noGTF_auc, paste0(sample_data$sample, '_auc.out'))
# prepare list for reading
list_df_noGTF <- with(sample_data, split(paths_noGTF, experiment_id))

for (i in sample_data$experiment_id) {
  j <- as.list(unlist(as.list(list_df_noGTF[i])))
  names(j) <- sample_data$sample[sample_data$experiment_id==i]
  assign(i, list(j))
}

final_list_noGTF <- list()
for (i in ls(pattern='SRP')) {
  final_list_noGTF <- append(final_list_noGTF, get(i))
}
names(final_list_noGTF) <- unique(sample_data$experiment_id)

#### no GTF loop ####
for (i in final_list_noGTF) {
  cat('Analysis of ', unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  
  subDir <- as.character(unique(sample_data$experiment_id[sample_data$sample==names(i)]))
  dir.create(file.path(subDir, 'noGTF'), showWarnings = FALSE)
  
  design <- data.frame('sample'= names(i), 'condition'= as.factor(sample_data$EXP_CTL[sample_data$sample==names(i)]))
  
  merged_df <- as.data.frame(map_dfc(i,
                                     ~vroom(.x,
                                            col_select=c(4),
                                            delim = '\t',
                                            col_names = FALSE,
                                            show_col_types = FALSE)))
  #renaming cols
  suppressWarnings(colnames(merged_df) <- names(i))
  suppressWarnings(merged_df <- head(merged_df, -2)) #remove last two lines
  rownames(merged_df) <- genenames$gene_id #gene ids as rownames
  
  cpm_counts <- cpm(merged_df)
  target_counts <- cpm_counts[rownames(cpm_counts)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]

  k <- 2
  result <- Ckmeans.1d.dp(target_counts, k)
  colors <- brewer.pal(k, "Dark2")
  midpoints <- ahist(result, style="midpoints", data=target_counts, plot=FALSE)$breaks[2:k]
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'noGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'clustering', 'noGTF.pdf', sep = '_')), height = 2.5)
  plot(result, 
            col.clusters = colors, 
            main = paste('K-means clustering of', unique(unique(sample_data$experiment_id[sample_data$paths_noGTF==i]))),
            ylab = '',
            yaxt = 'n',
            xlab = paste(unique(sample_data$targeted_gene[sample_data$sample==names(i)]), 'CPM', '(noGTF)'),
            ylim=c(0,1))
  abline(v=midpoints, col="RoyalBlue", lwd=3, lty=3,)
  text(x = target_counts, y = par("usr")[3]-0.15, labels = paste(names(i), design$EXP_CTL[design$sample==names(i)]), xpd = NA, srt=45, cex=0.4, adj = 0.9)
  dev.off()
  
  ################
  #### DESEQ2 ####
  ################
  
  dds <- DESeqDataSetFromMatrix(countData = merged_df, 
                                colData = design, 
                                design= ~ condition)
  ## FILTER ##
  nrow(dds)
  keep <- rowSums(counts(dds) >= 1) >= length(design$sample)/2
  dds <- dds[keep,]
  nrow(dds)
  
  ## VST COUNT TRANSFORMATION n PLOT ##
  vsd <- vst(dds, blind = FALSE)
  
  ## PCA ##
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'noGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'PCA',  'noGTF.pdf', sep = '_')))
  print(plotPCA(vsd, intgroup = c('condition')))
  dev.off()
  
  ## DEA ##
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  summary(res)
  
  write.csv(as.data.frame(res), file = paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'noGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'noGTF', 'DE.tsv', sep = '_')), quote = F, sep = '\t')
  
  ## TARGET GENE COUNTS PLOT ##
  targetGene <- rownames(res)[rownames(res)==unique(sample_data$targeted_gene[sample_data$sample==names(i)])]
  geneCounts <- plotCounts(dds, 
                           gene = targetGene, 
                           intgroup = c('condition'),
                           transform = F,
                           returnData = TRUE)
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'noGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), targetGene, 'noGTF.pdf', sep = '_')))
  print(ggplot(geneCounts, aes(x = condition, y = count, color = condition, size=4)) +  
    theme_bw() +
    ggtitle(label = paste(targetGene, 'raw counts')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_beeswarm(cex = 3) +
    annotate(geom="text", 
             x=Inf, 
             y=res[rownames(res)==targetGene,]$baseMean, 
             label=paste('log2FC:', 
                         round(res[rownames(res)==targetGene,]$log2FoldChange, digits=3), '\n',
                         'pvalue:', 
                         round(res[rownames(res)==targetGene,]$pvalue, digits=3), '\n',
                         'padj:', 
                         round(res[rownames(res)==targetGene,]$padj, digits=3), '\n'),
             hjust=1, vjust=1,
             color="black") + 
    guides(condition = "legend", size = "none"))
  dev.off()
  
  pdf(paste0(unique(sample_data$experiment_id[sample_data$sample==names(i)]), '/', 'noGTF/', paste(unique(sample_data$experiment_id[sample_data$sample==names(i)]), 'MAplot', 'noGTF.pdf', sep = '_')))
  plotMA(res,ylim = c(-max(res$log2FoldChange), max(res$log2FoldChange)))
  with(res[targetGene, ], {
    points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
    text(baseMean, log2FoldChange, targetGene, pos=2, col="red")
  })
  dev.off()
  
}

