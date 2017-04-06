library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(DESeq2)
library("scales")
library(ggdendro)

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

args = c('sample_key_all_kinase.txt', 'count_matrix.txt', 'global_alignment_report.txt')

##read in data
sample_key = read.table(args[1], header = T, sep = '\t', stringsAsFactors = F)
sample_key$group = sapply(1:nrow(sample_key), function(x) {
	paste(sample_key[x,2:5], collapse = '_')
})
counts = read.table(args[2], header = T, sep = '\t',row.names = 1, stringsAsFactors = F)
align_rep = read.table(args[3], header = T, sep = '\t', stringsAsFactors = F)

##stacked barplot of perfect/good/poor/unmapped reads per sample
align_plot = stack(align_rep[,c('perfect','good','poor','unmapped')])
align_plot$sample = rep(align_rep$sample, times = 4)
colnames(align_plot) = c('count','quality','sample')
align_plot$quality = factor(align_plot$quality, levels = c('perfect', 'good', 'poor', 'unmapped'))
tiff(file = 'alignment_quality.tiff', res = 300, width = 3000, height = 2000)
ggplot(align_plot, aes(sample, fill = quality, weight = count)) + 
	geom_bar(color = 'black') + 
	scale_fill_brewer(name = 'Quality', type = 'seq', palette = 1, direction = -1) +
	scale_y_continuous(expand = c(.01,0)) +
	labs(title = 'Alignment Quality Summary', x = 'Sample', y = 'Number of Reads') +
	theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()
	
##boxplot of alignment quality distribution across samples
align_plot$fraction = mapply(function(x,y) signif(x / sum(align_plot$count[align_plot$sample == y]), 3), align_plot$count, align_plot$sample)
# define the functions to calculate boxplot quantiles
boxplot_stats = function(x) {
	quantiles = quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
	names(quantiles) = c("ymin", "lower", "middle", "upper", "ymax")
	quantiles
}
outlier_stats = function(x) {
	outliers = quantile(x, probs = c(.1, .9))
	subset(x, x < outliers[1] | outliers[2] < x)
}
png(file = 'alignment_quality_distribution.png', width = 1000, height = 750)
ggplot(align_plot, aes(quality, fraction)) + 
	stat_summary(fun.data = boxplot_stats, geom = "boxplot", fill = 'lightcoral') + 
	stat_summary(fun.y = outlier_stats, geom = "point", alpha = 1/2) +
	coord_cartesian(ylim = c(0, 1)) +
	scale_y_continuous(expand = c(.01,0)) +
	labs(title = 'Distribution of Read Alignment Quality', x = 'Quality', y = 'Fraction of Reads per Sample') +
	theme_bw()
dev.off()

##boxplot of reads/hairpin/sample
hairpin_plot = stack(counts)
hairpin_plot$shrna = rep(rownames(counts), times = ncol(counts))
colnames(hairpin_plot) = c('count','sample','shrna')
png(file = 'counts_per_hp_raw.png', width = 1000, height = 750)
ggplot(hairpin_plot, aes(sample, count)) + 
	stat_summary(fun.data = boxplot_stats, geom = "boxplot", fill = 'lightblue') + 
	stat_summary(fun.y = outlier_stats, geom = "point", alpha = 1/5) +
	coord_cartesian(ylim = c(0, 3000)) +
	scale_y_continuous(expand = c(.01,0)) +
	labs(title = 'Distribution of Aligned Read Count per Hairpin', x = 'Sample', y = 'Number of Aligned Reads') +
	theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##boxplot of normalized reads/hairpin/sample
dds = DESeqDataSetFromMatrix(counts, colData = data.frame(rep(1, times = ncol(counts))), design = ~ 1)
dds = estimateSizeFactors(dds)
counts_norm = as.data.frame(round(counts(dds, normalized = T)))
sample_names = paste(sapply(colnames(counts), function(x) {strsplit(x, '_')[[1]][2]}), sample_key$group, sep = '_')
colnames(counts_norm) = sample_names
hairpin_norm = stack(counts_norm)
hairpin_norm$shrna = rep(rownames(counts_norm), times = ncol(counts_norm))
colnames(hairpin_norm) = c('count','sample','shrna')
png(file = 'counts_per_hp_norm.png', width = 1000, height = 750)
ggplot(hairpin_norm, aes(sample, count)) + 
	stat_summary(fun.data = boxplot_stats, geom = "boxplot", fill = 'lightblue') + 
	stat_summary(fun.y = outlier_stats, geom = "point", alpha = 1/5) +
	coord_cartesian(ylim = c(0, 3000)) +
	scale_y_continuous(expand = c(.01,0)) +
	labs(title = 'Distribution of Normalized Read Count per Hairpin', x = 'Sample', y = 'Number of Aligned Reads') +
	theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()
	
##number of hairpins below various count thresholds
thresholds = c(0,1,5,10,25,50,100,250,500,1000, 2500, 5000, 10000, max(counts))
nc = ncol(counts_norm)
nh = nrow(counts_norm)
N_samples = c(1, nc / 2, nc)
tbl = sapply(N_samples, function(x) sapply(thresholds, function(y) signif(sum(rowSums(counts_norm <= y) >= x) / nh, 5)))
colnames(tbl) = N_samples
rownames(tbl) = thresholds
tbl = as.data.frame(tbl)
tbl_stack = stack(tbl)
tbl_stack$threshold = factor(rep(thresholds, times = length(N_samples)))
colnames(tbl_stack) = c('fraction', 'samples', 'threshold')
tbl_stack$samples = factor(tbl_stack$samples, levels = N_samples)
png(file = 'counts_at_thresholds.png', width = 1000, height = 750)
ggplot(tbl_stack, aes(threshold, fraction, group = samples, color = samples)) +
	geom_line(size = 1.1) +
	geom_point(size = 3) +
	scale_y_continuous(expand = c(.01, 0)) +
	scale_color_brewer(name = 'Samples', type = 'qual', palette = 6) +
	labs(title = 'Fraction of Hairpins Below Count Threshold in N Samples', x = 'Count Threshold', y = 'Fraction of Hairpins') +
	theme_bw()
dev.off()

##normalized read count correlation matrix
corr_mat = signif(cor(counts_norm, counts_norm, method = 'p'), 2)
corr_stack = stack(split(corr_mat, col(corr_mat)))
corr_stack$y = rep(rownames(corr_mat), times = ncol(corr_mat))
colnames(corr_stack) = c('PCC', 'x' ,'y')
corr_stack$x = rep(colnames(corr_mat), each = nrow(corr_mat))
png(file = 'count_correlation_matrix.png', width = 1000, height = 750)
ggplot(corr_stack, aes(x, y)) + 
	geom_tile(aes(fill = PCC), colour = "white") +
	scale_fill_gradient(low = "white", high = "firebrick3") +
	labs(title = 'Pearson Correlation Coefficients of Normalized Counts', x = '', y = '') +
	theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

##normalized read count dendrogram
hc = hclust(dist(t(counts_norm)), method = 'ward.D2')
hc_data <- dendro_data(hc)
png(file = 'count_dendrogram.png', width = 1000, height = 750)
ggdendrogram(hc_data, rotate = F, size = 2) + 
	labs(title = "Dendrogram of Normalized Counts")
dev.off()

##read and format sample metadata
sample_key = read.table(args[1], header = T, sep = '\t', colClasses = c(rep('character',5), 'factor'))
sample_key$group = as.factor(sapply(1:nrow(sample_key), function(x) paste(sample_key[x,2:5],collapse = '_')))
sample_key = sample_key[1:80,]
sample_key = sample_key[c(-30,-38),]
unique_groups = unique(as.character(sample_key$group))

##read count data
counts = as.matrix(read.table(args[2], header = T, sep = '\t', stringsAsFactors = F))
counts = counts[,1:80]
counts = counts[,c(-30,-38)]
colnames(counts) = paste(colnames(counts), sample_key$group, sep = '_')

shrna = sapply(rownames(counts), function(x) {
	strsplit(as.character(x), '-', fixed = T)[[1]][1]
})
gene = sapply(rownames(counts), function(x) {
	strsplit(as.character(x), '-', fixed = T)[[1]][2]
})

## perform DESeq and extract results
gata3_veh_t0_counts = round(counts[,c(65:68,1:4)])
gata3_veh_t0_dds = DESeqDataSetFromMatrix(gata3_veh_t0_counts, colData = sample_key[c(65:68,1:4),], design = ~0 + group + replicate)
gata3_veh_t0_dds = DESeq(gata3_veh_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_veh_t0_counts_norm = cbind(shrna, gene, round(counts(gata3_veh_t0_dds, normalized = T)))
gata3_veh_t0_results = results(gata3_veh_t0_dds, contrast = c('group', unique_groups[1], unique_groups[17]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_dox_t0_counts = round(counts[,c(65:68,5:8)])
gata3_dox_t0_dds = DESeqDataSetFromMatrix(gata3_dox_t0_counts, colData = sample_key[c(65:68,5:8),], design = ~0 + group + replicate)
gata3_dox_t0_dds = DESeq(gata3_dox_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_dox_t0_counts_norm = cbind(shrna, gene, round(counts(gata3_dox_t0_dds, normalized = T)))
gata3_dox_t0_results = results(gata3_dox_t0_dds, contrast = c('group', unique_groups[2], unique_groups[17]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_tam_t0_counts = round(counts[,c(65:68,9:12)])
gata3_tam_t0_dds = DESeqDataSetFromMatrix(gata3_tam_t0_counts, colData = sample_key[c(65:68,9:12),], design = ~0 + group + replicate)
gata3_tam_t0_dds = DESeq(gata3_tam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_tam_t0_counts_norm = cbind(shrna, gene, round(counts(gata3_tam_t0_dds, normalized = T)))
gata3_tam_t0_results = results(gata3_tam_t0_dds, contrast = c('group', unique_groups[3], unique_groups[17]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_doxtam_t0_counts = round(counts[,c(65:68,13:16)])
gata3_doxtam_t0_dds = DESeqDataSetFromMatrix(gata3_doxtam_t0_counts, colData = sample_key[c(65:68,13:16),], design = ~0 + group + replicate)
gata3_doxtam_t0_dds = DESeq(gata3_doxtam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_doxtam_t0_counts_norm = cbind(shrna, gene, round(counts(gata3_doxtam_t0_dds, normalized = T)))
gata3_doxtam_t0_results = results(gata3_doxtam_t0_dds, contrast = c('group', unique_groups[4], unique_groups[17]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_dox_unt_counts = round(counts[,c(1:4,5:8)])
gata3_dox_unt_dds = DESeqDataSetFromMatrix(gata3_dox_unt_counts, colData = sample_key[c(1:4,5:8),], design = ~0 + group + replicate)
gata3_dox_unt_dds = DESeq(gata3_dox_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_dox_unt_counts_norm = cbind(shrna, gene, round(counts(gata3_dox_unt_dds, normalized = T)))
gata3_dox_unt_results = results(gata3_dox_unt_dds, contrast = c('group', unique_groups[2], unique_groups[1]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_tam_unt_counts = round(counts[,c(1:4,9:12)])
gata3_tam_unt_dds = DESeqDataSetFromMatrix(gata3_tam_unt_counts, colData = sample_key[c(1:4,9:12),], design = ~0 + group + replicate)
gata3_tam_unt_dds = DESeq(gata3_tam_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_tam_unt_counts_norm = cbind(shrna, gene, round(counts(gata3_tam_unt_dds, normalized = T)))
gata3_tam_unt_results = results(gata3_tam_unt_dds, contrast = c('group', unique_groups[3], unique_groups[1]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

gata3_doxtam_dox_counts = round(counts[,c(5:8,13:16)])
gata3_doxtam_dox_dds = DESeqDataSetFromMatrix(gata3_doxtam_dox_counts, colData = sample_key[c(5:8,13:16),], design = ~0 + group + replicate)
gata3_doxtam_dox_dds = DESeq(gata3_doxtam_dox_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
gata3_doxtam_dox_counts_norm = cbind(shrna, gene, round(counts(gata3_doxtam_dox_dds, normalized = T)))
gata3_doxtam_dox_results = results(gata3_doxtam_dox_dds, contrast = c('group', unique_groups[4], unique_groups[2]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_veh_t0_counts = round(counts[,c(67:70,17:20)])
ptpn12_veh_t0_dds = DESeqDataSetFromMatrix(ptpn12_veh_t0_counts, colData = sample_key[c(67:70,17:20),], design = ~0 + group + replicate)
ptpn12_veh_t0_dds = DESeq(ptpn12_veh_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_veh_t0_counts_norm = cbind(shrna, gene, round(counts(ptpn12_veh_t0_dds, normalized = T)))
ptpn12_veh_t0_results = results(ptpn12_veh_t0_dds, contrast = c('group', unique_groups[5], unique_groups[18]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_dox_t0_counts = round(counts[,c(67:70,21:24)])
ptpn12_dox_t0_dds = DESeqDataSetFromMatrix(ptpn12_dox_t0_counts, colData = sample_key[c(67:70,21:24),], design = ~0 + group + replicate)
ptpn12_dox_t0_dds = DESeq(ptpn12_dox_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_dox_t0_counts_norm = cbind(shrna, gene, round(counts(ptpn12_dox_t0_dds, normalized = T)))
ptpn12_dox_t0_results = results(ptpn12_dox_t0_dds, contrast = c('group', unique_groups[6], unique_groups[18]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_tam_t0_counts = round(counts[,c(67:70,25:28)])
ptpn12_tam_t0_dds = DESeqDataSetFromMatrix(ptpn12_tam_t0_counts, colData = sample_key[c(67:70,25:28),], design = ~0 + group + replicate)
ptpn12_tam_t0_dds = DESeq(ptpn12_tam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_tam_t0_counts_norm = cbind(shrna, gene, round(counts(ptpn12_tam_t0_dds, normalized = T)))
ptpn12_tam_t0_results = results(ptpn12_tam_t0_dds, contrast = c('group', unique_groups[7], unique_groups[18]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_doxtam_t0_counts = round(counts[,c(67:70,29:31)])
ptpn12_doxtam_t0_dds = DESeqDataSetFromMatrix(ptpn12_doxtam_t0_counts, colData = sample_key[c(67:70,29:31),], design = ~0 + group + replicate)
ptpn12_doxtam_t0_dds = DESeq(ptpn12_doxtam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_doxtam_t0_counts_norm = cbind(shrna, gene, round(counts(ptpn12_doxtam_t0_dds, normalized = T)))
ptpn12_doxtam_t0_results = results(ptpn12_doxtam_t0_dds, contrast = c('group', unique_groups[8], unique_groups[18]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_dox_unt_counts = round(counts[,c(17:20,21:24)])
ptpn12_dox_unt_dds = DESeqDataSetFromMatrix(ptpn12_dox_unt_counts, colData = sample_key[c(17:20,21:24),], design = ~0 + group + replicate)
ptpn12_dox_unt_dds = DESeq(ptpn12_dox_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_dox_unt_counts_norm = cbind(shrna, gene, round(counts(ptpn12_dox_unt_dds, normalized = T)))
ptpn12_dox_unt_results = results(ptpn12_dox_unt_dds, contrast = c('group', unique_groups[6], unique_groups[5]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_tam_unt_counts = round(counts[,c(17:20,25:28)])
ptpn12_tam_unt_dds = DESeqDataSetFromMatrix(ptpn12_tam_unt_counts, colData = sample_key[c(17:20,25:28),], design = ~0 + group + replicate)
ptpn12_tam_unt_dds = DESeq(ptpn12_tam_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_tam_unt_counts_norm = cbind(shrna, gene, round(counts(ptpn12_tam_unt_dds, normalized = T)))
ptpn12_tam_unt_results = results(ptpn12_tam_unt_dds, contrast = c('group', unique_groups[7], unique_groups[5]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

ptpn12_doxtam_dox_counts = round(counts[,c(21:24,29:31)])
ptpn12_doxtam_dox_dds = DESeqDataSetFromMatrix(ptpn12_doxtam_dox_counts, colData = sample_key[c(21:24,29:31),], design = ~0 + group + replicate)
ptpn12_doxtam_dox_dds = DESeq(ptpn12_doxtam_dox_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
ptpn12_doxtam_dox_counts_norm = cbind(shrna, gene, round(counts(ptpn12_doxtam_dox_dds, normalized = T)))
ptpn12_doxtam_dox_results = results(ptpn12_doxtam_dox_dds, contrast = c('group', unique_groups[8], unique_groups[6]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_veh_t0_counts = round(counts[,c(71:74,32:35)])
inpp4b_veh_t0_dds = DESeqDataSetFromMatrix(inpp4b_veh_t0_counts, colData = sample_key[c(71:74,32:35),], design = ~0 + group + replicate)
inpp4b_veh_t0_dds = DESeq(inpp4b_veh_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_veh_t0_counts_norm = cbind(shrna, gene, round(counts(inpp4b_veh_t0_dds, normalized = T)))
inpp4b_veh_t0_results = results(inpp4b_veh_t0_dds, contrast = c('group', unique_groups[9], unique_groups[19]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_dox_t0_counts = round(counts[,c(71:74,36:38)])
inpp4b_dox_t0_dds = DESeqDataSetFromMatrix(inpp4b_dox_t0_counts, colData = sample_key[c(71:74,36:38),], design = ~0 + group + replicate)
inpp4b_dox_t0_dds = DESeq(inpp4b_dox_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_dox_t0_counts_norm = cbind(shrna, gene, round(counts(inpp4b_dox_t0_dds, normalized = T)))
inpp4b_dox_t0_results = results(inpp4b_dox_t0_dds, contrast = c('group', unique_groups[10], unique_groups[19]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_tam_t0_counts = round(counts[,c(71:74,39:42)])
inpp4b_tam_t0_dds = DESeqDataSetFromMatrix(inpp4b_tam_t0_counts, colData = sample_key[c(71:74,39:42),], design = ~0 + group + replicate)
inpp4b_tam_t0_dds = DESeq(inpp4b_tam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_tam_t0_counts_norm = cbind(shrna, gene, round(counts(inpp4b_tam_t0_dds, normalized = T)))
inpp4b_tam_t0_results = results(inpp4b_tam_t0_dds, contrast = c('group', unique_groups[11], unique_groups[19]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_doxtam_t0_counts = round(counts[,c(71:74,43:46)])
inpp4b_doxtam_t0_dds = DESeqDataSetFromMatrix(inpp4b_doxtam_t0_counts, colData = sample_key[c(71:74,43:46),], design = ~0 + group + replicate)
inpp4b_doxtam_t0_dds = DESeq(inpp4b_doxtam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_doxtam_t0_counts_norm = cbind(shrna, gene, round(counts(inpp4b_doxtam_t0_dds, normalized = T)))
inpp4b_doxtam_t0_results = results(inpp4b_doxtam_t0_dds, contrast = c('group', unique_groups[12], unique_groups[19]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_dox_unt_counts = round(counts[,c(32:35,36:38)])
inpp4b_dox_unt_dds = DESeqDataSetFromMatrix(inpp4b_dox_unt_counts, colData = sample_key[c(32:35,36:38),], design = ~0 + group + replicate)
inpp4b_dox_unt_dds = DESeq(inpp4b_dox_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_dox_unt_counts_norm = cbind(shrna, gene, round(counts(inpp4b_dox_unt_dds, normalized = T)))
inpp4b_dox_unt_results = results(inpp4b_dox_unt_dds, contrast = c('group', unique_groups[10], unique_groups[9]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_tam_unt_counts = round(counts[,c(32:35,39:42)])
inpp4b_tam_unt_dds = DESeqDataSetFromMatrix(inpp4b_tam_unt_counts, colData = sample_key[c(32:35,39:42),], design = ~0 + group + replicate)
inpp4b_tam_unt_dds = DESeq(inpp4b_tam_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_tam_unt_counts_norm = cbind(shrna, gene, round(counts(inpp4b_tam_unt_dds, normalized = T)))
inpp4b_tam_unt_results = results(inpp4b_tam_unt_dds, contrast = c('group', unique_groups[11], unique_groups[9]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

inpp4b_doxtam_dox_counts = round(counts[,c(36:38,43:46)])
inpp4b_doxtam_dox_dds = DESeqDataSetFromMatrix(inpp4b_doxtam_dox_counts, colData = sample_key[c(36:38,43:46),], design = ~0 + group + replicate)
inpp4b_doxtam_dox_dds = DESeq(inpp4b_doxtam_dox_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
inpp4b_doxtam_dox_counts_norm = cbind(shrna, gene, round(counts(inpp4b_doxtam_dox_dds, normalized = T)))
inpp4b_doxtam_dox_results = results(inpp4b_doxtam_dox_dds, contrast = c('group', unique_groups[12], unique_groups[10]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_veh_t0_counts = round(counts[,c(75:78,47:50)])
map3k1_veh_t0_dds = DESeqDataSetFromMatrix(map3k1_veh_t0_counts, colData = sample_key[c(75:78,47:50),], design = ~0 + group + replicate)
map3k1_veh_t0_dds = DESeq(map3k1_veh_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_veh_t0_counts_norm = cbind(shrna, gene, round(counts(map3k1_veh_t0_dds, normalized = T)))
map3k1_veh_t0_results = results(map3k1_veh_t0_dds, contrast = c('group', unique_groups[13], unique_groups[20]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_dox_t0_counts = round(counts[,c(75:78,51:54)])
map3k1_dox_t0_dds = DESeqDataSetFromMatrix(map3k1_dox_t0_counts, colData = sample_key[c(75:78,51:54),], design = ~0 + group + replicate)
map3k1_dox_t0_dds = DESeq(map3k1_dox_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_dox_t0_counts_norm = cbind(shrna, gene, round(counts(map3k1_dox_t0_dds, normalized = T)))
map3k1_dox_t0_results = results(map3k1_dox_t0_dds, contrast = c('group', unique_groups[14], unique_groups[20]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_tam_t0_counts = round(counts[,c(75:78,55:58)])
map3k1_tam_t0_dds = DESeqDataSetFromMatrix(map3k1_tam_t0_counts, colData = sample_key[c(75:78,55:58),], design = ~0 + group + replicate)
map3k1_tam_t0_dds = DESeq(map3k1_tam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_tam_t0_counts_norm = cbind(shrna, gene, round(counts(map3k1_tam_t0_dds, normalized = T)))
map3k1_tam_t0_results = results(map3k1_tam_t0_dds, contrast = c('group', unique_groups[15], unique_groups[20]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_doxtam_t0_counts = round(counts[,c(75:78,59:62)])
map3k1_doxtam_t0_dds = DESeqDataSetFromMatrix(map3k1_doxtam_t0_counts, colData = sample_key[c(75:78,59:62),], design = ~0 + group + replicate)
map3k1_doxtam_t0_dds = DESeq(map3k1_doxtam_t0_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_doxtam_t0_counts_norm = cbind(shrna, gene, round(counts(map3k1_doxtam_t0_dds, normalized = T)))
map3k1_doxtam_t0_results = results(map3k1_doxtam_t0_dds, contrast = c('group', unique_groups[16], unique_groups[20]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_dox_unt_counts = round(counts[,c(47:50,51:54)])
map3k1_dox_unt_dds = DESeqDataSetFromMatrix(map3k1_dox_unt_counts, colData = sample_key[c(47:50,51:54),], design = ~0 + group + replicate)
map3k1_dox_unt_dds = DESeq(map3k1_dox_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_dox_unt_counts_norm = cbind(shrna, gene, round(counts(map3k1_dox_unt_dds, normalized = T)))
map3k1_dox_unt_results = results(map3k1_dox_unt_dds, contrast = c('group', unique_groups[14], unique_groups[13]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_tam_unt_counts = round(counts[,c(47:50,55:58)])
map3k1_tam_unt_dds = DESeqDataSetFromMatrix(map3k1_tam_unt_counts, colData = sample_key[c(47:50,55:58),], design = ~0 + group + replicate)
map3k1_tam_unt_dds = DESeq(map3k1_tam_unt_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_tam_unt_counts_norm = cbind(shrna, gene, round(counts(map3k1_tam_unt_dds, normalized = T)))
map3k1_tam_unt_results = results(map3k1_tam_unt_dds, contrast = c('group', unique_groups[15], unique_groups[13]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

map3k1_doxtam_dox_counts = round(counts[,c(51:54,59:62)])
map3k1_doxtam_dox_dds = DESeqDataSetFromMatrix(map3k1_doxtam_dox_counts, colData = sample_key[c(51:54,59:62),], design = ~0 + group + replicate)
map3k1_doxtam_dox_dds = DESeq(map3k1_doxtam_dox_dds, fitType = 'local', minReplicatesForReplace = Inf, betaPrior = F)
map3k1_doxtam_dox_counts_norm = cbind(shrna, gene, round(counts(map3k1_doxtam_dox_dds, normalized = T)))
map3k1_doxtam_dox_results = results(map3k1_doxtam_dox_dds, contrast = c('group', unique_groups[16], unique_groups[14]), cooksCutoff = Inf, independentFiltering = F, pAdjustMethod = "none")

results = list(
	gata3_veh_t0 = gata3_veh_t0_results,
	gata3_dox_t0 = gata3_dox_t0_results,
	gata3_tam_t0 = gata3_tam_t0_results,
	gata3_doxtam_t0 = gata3_doxtam_t0_results,
	gata3_dox_unt = gata3_dox_unt_results,
	gata3_tam_unt = gata3_tam_unt_results,
	gata3_doxtam_dox = gata3_doxtam_dox_results,
	ptpn12_veh_t0 = ptpn12_veh_t0_results,
	ptpn12_dox_t0 = ptpn12_dox_t0_results,
	ptpn12_tam_t0 = ptpn12_tam_t0_results,
	ptpn12_doxtam_t0 = ptpn12_doxtam_t0_results,
	ptpn12_dox_unt = ptpn12_dox_unt_results,
	ptpn12_tam_unt = ptpn12_tam_unt_results,
	ptpn12_doxtam_dox = ptpn12_doxtam_dox_results,
	inpp4b_veh_t0 = inpp4b_veh_t0_results,
	inpp4b_dox_t0 = inpp4b_dox_t0_results,
	inpp4b_tam_t0 = inpp4b_tam_t0_results,
	inpp4b_doxtam_t0 = inpp4b_doxtam_t0_results,
	inpp4b_dox_unt = inpp4b_dox_unt_results,
	inpp4b_tam_unt = inpp4b_tam_unt_results,
	inpp4b_doxtam_dox = inpp4b_doxtam_dox_results,
	map3k1_veh_t0 = map3k1_veh_t0_results,
	map3k1_dox_t0 = map3k1_dox_t0_results,
	map3k1_tam_t0 = map3k1_tam_t0_results,
	map3k1_doxtam_t0 = map3k1_doxtam_t0_results,
	map3k1_dox_unt = map3k1_dox_unt_results,
	map3k1_tam_unt = map3k1_tam_unt_results,
	map3k1_doxtam_dox = map3k1_doxtam_dox_results
)

names = names(results)
int_results = do.call(cbind, as.matrix((unlist(results))))
colnames(int_results) = gsub('log2FoldChange', 'l2fc', colnames(int_results))
colnames(int_results) = sapply(colnames(int_results), function(x) strsplit(as.character(x), '.', fixed = T)[[1]][1])
colnames(int_results) = paste(colnames(int_results), rep(names, each = 6), sep = '_')
cols = unlist(lapply(6 * seq(0,length(names) - 1), function(x) x + c(2,3,6)))
int_results = as.data.frame(int_results[,cols])
int_results = apply(int_results, 2, function(x) signif(x, 3))
final_results = data.frame(shRNA = shrna, gene = gene, int_results)
write.table(final_results, file = 'deseq_results.txt', row.names = F, sep = '\t')

###################################################################
##counting the number of lethal/beneficial hairpins in each class##
###################################################################

na_rows = which(rowSums(is.na(final_results)) != 0)
final_results = final_results[-na_rows,]

l_thr = c(.5, .5, .5, .5, .5, .5)
p_thr = c(.1, .1, .1, .1, .1, .1)

final_results$gata3_veh_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_veh_t0']
	pval = final_results[x,'padj_gata3_veh_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_dox_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_dox_t0']
	pval = final_results[x,'padj_gata3_dox_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_tam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_tam_t0']
	pval = final_results[x,'padj_gata3_tam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_doxtam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_doxtam_t0']
	pval = final_results[x,'padj_gata3_doxtam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_dox_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_dox_unt']
	pval = final_results[x,'padj_gata3_dox_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_tam_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_tam_unt']
	pval = final_results[x,'padj_gata3_tam_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$gata3_doxtam_dox_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_gata3_doxtam_dox']
	pval = final_results[x,'padj_gata3_doxtam_dox']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_veh_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_veh_t0']
	pval = final_results[x,'padj_ptpn12_veh_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_dox_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_dox_t0']
	pval = final_results[x,'padj_ptpn12_dox_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_tam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_tam_t0']
	pval = final_results[x,'padj_ptpn12_tam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_doxtam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_doxtam_t0']
	pval = final_results[x,'padj_ptpn12_doxtam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_dox_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_dox_unt']
	pval = final_results[x,'padj_ptpn12_dox_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_tam_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_tam_unt']
	pval = final_results[x,'padj_ptpn12_tam_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$ptpn12_doxtam_dox_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_ptpn12_doxtam_dox']
	pval = final_results[x,'padj_ptpn12_doxtam_dox']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_veh_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_veh_t0']
	pval = final_results[x,'padj_inpp4b_veh_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_dox_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_dox_t0']
	pval = final_results[x,'padj_inpp4b_dox_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_tam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_tam_t0']
	pval = final_results[x,'padj_inpp4b_tam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_doxtam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_doxtam_t0']
	pval = final_results[x,'padj_inpp4b_doxtam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_dox_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_dox_unt']
	pval = final_results[x,'padj_inpp4b_dox_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_tam_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_tam_unt']
	pval = final_results[x,'padj_inpp4b_tam_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$inpp4b_doxtam_dox_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_inpp4b_doxtam_dox']
	pval = final_results[x,'padj_inpp4b_doxtam_dox']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_veh_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_veh_t0']
	pval = final_results[x,'padj_map3k1_veh_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_dox_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_dox_t0']
	pval = final_results[x,'padj_map3k1_dox_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_tam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_tam_t0']
	pval = final_results[x,'padj_map3k1_tam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_doxtam_t0_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_doxtam_t0']
	pval = final_results[x,'padj_map3k1_doxtam_t0']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_dox_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_dox_unt']
	pval = final_results[x,'padj_map3k1_dox_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_tam_unt_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_tam_unt']
	pval = final_results[x,'padj_map3k1_tam_unt']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$map3k1_doxtam_dox_class = sapply(1:nrow(final_results), function(x) {
	l2fc = final_results[x,'l2fc_map3k1_doxtam_dox']
	pval = final_results[x,'padj_map3k1_doxtam_dox']
	if((l2fc < -l_thr[1]) & (pval < p_thr[1])) {
		'lethal'
	} else if((l2fc > l_thr[1]) & (pval < p_thr[1])) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$overall_veh_t0_class = sapply(1:nrow(final_results), function(x) {
	cls = final_results[x,grep('_veh_t0_class', colnames(final_results))]
	if(sum(cls == 'lethal') >= 3) {
		'lethal'
	} else if(sum(cls == 'beneficial') >= 3) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$overall_tam_t0_class = sapply(1:nrow(final_results), function(x) {
	cls = final_results[x,grep('_tam_t0_class', colnames(final_results))]
	if(sum(cls == 'lethal') >= 3) {
		'lethal'
	} else if(sum(cls == 'beneficial') >= 3) {
		'beneficial'
	} else {
		'neutral'
	}
})
final_results$overall_tam_unt_class = sapply(1:nrow(final_results), function(x) {
	cls = final_results[x,grep('_tam_unt_class', colnames(final_results))]
	if(sum(cls == 'lethal') >= 3) {
		'lethal'
	} else if(sum(cls == 'beneficial') >= 3) {
		'beneficial'
	} else {
		'neutral'
	}
})

class_count_table = t(sapply(final_results$gene, function(x) {
	gata3_veh_t0_lethal = sum(final_results$gene == x & final_results$gata3_veh_t0_class == 'lethal')
	gata3_veh_t0_beneficial = sum(final_results$gene == x & final_results$gata3_veh_t0_class == 'beneficial')
	gata3_dox_t0_lethal = sum(final_results$gene == x & final_results$gata3_dox_t0_class == 'lethal')
	gata3_dox_t0_beneficial = sum(final_results$gene == x & final_results$gata3_dox_t0_class == 'beneficial')
	gata3_tam_t0_lethal = sum(final_results$gene == x & final_results$gata3_tam_t0_class == 'lethal')
	gata3_tam_t0_beneficial = sum(final_results$gene == x & final_results$gata3_tam_t0_class == 'beneficial')
	gata3_doxtam_t0_lethal = sum(final_results$gene == x & final_results$gata3_doxtam_t0_class == 'lethal')
	gata3_doxtam_t0_beneficial = sum(final_results$gene == x & final_results$gata3_doxtam_t0_class == 'beneficial')
	gata3_dox_unt_lethal = sum(final_results$gene == x & final_results$gata3_dox_unt_class == 'lethal')
	gata3_dox_unt_beneficial = sum(final_results$gene == x & final_results$gata3_dox_unt_class == 'beneficial')
	gata3_tam_unt_lethal = sum(final_results$gene == x & final_results$gata3_tam_unt_class == 'lethal')
	gata3_tam_unt_beneficial = sum(final_results$gene == x & final_results$gata3_tam_unt_class == 'beneficial')
	gata3_doxtam_dox_lethal = sum(final_results$gene == x & final_results$gata3_doxtam_dox_class == 'lethal')
	gata3_doxtam_dox_beneficial = sum(final_results$gene == x & final_results$gata3_doxtam_dox_class == 'beneficial')
	ptpn12_veh_t0_lethal = sum(final_results$gene == x & final_results$ptpn12_veh_t0_class == 'lethal')
	ptpn12_veh_t0_beneficial = sum(final_results$gene == x & final_results$ptpn12_veh_t0_class == 'beneficial')
	ptpn12_dox_t0_lethal = sum(final_results$gene == x & final_results$ptpn12_dox_t0_class == 'lethal')
	ptpn12_dox_t0_beneficial = sum(final_results$gene == x & final_results$ptpn12_dox_t0_class == 'beneficial')
	ptpn12_tam_t0_lethal = sum(final_results$gene == x & final_results$ptpn12_tam_t0_class == 'lethal')
	ptpn12_tam_t0_beneficial = sum(final_results$gene == x & final_results$ptpn12_tam_t0_class == 'beneficial')
	ptpn12_doxtam_t0_lethal = sum(final_results$gene == x & final_results$ptpn12_doxtam_t0_class == 'lethal')
	ptpn12_doxtam_t0_beneficial = sum(final_results$gene == x & final_results$ptpn12_doxtam_t0_class == 'beneficial')
	ptpn12_dox_unt_lethal = sum(final_results$gene == x & final_results$ptpn12_dox_unt_class == 'lethal')
	ptpn12_dox_unt_beneficial = sum(final_results$gene == x & final_results$ptpn12_dox_unt_class == 'beneficial')
	ptpn12_tam_unt_lethal = sum(final_results$gene == x & final_results$ptpn12_tam_unt_class == 'lethal')
	ptpn12_tam_unt_beneficial = sum(final_results$gene == x & final_results$ptpn12_tam_unt_class == 'beneficial')
	ptpn12_doxtam_dox_lethal = sum(final_results$gene == x & final_results$ptpn12_doxtam_dox_class == 'lethal')
	ptpn12_doxtam_dox_beneficial = sum(final_results$gene == x & final_results$ptpn12_doxtam_dox_class == 'beneficial')
	inpp4b_veh_t0_lethal = sum(final_results$gene == x & final_results$inpp4b_veh_t0_class == 'lethal')
	inpp4b_veh_t0_beneficial = sum(final_results$gene == x & final_results$inpp4b_veh_t0_class == 'beneficial')
	inpp4b_dox_t0_lethal = sum(final_results$gene == x & final_results$inpp4b_dox_t0_class == 'lethal')
	inpp4b_dox_t0_beneficial = sum(final_results$gene == x & final_results$inpp4b_dox_t0_class == 'beneficial')
	inpp4b_tam_t0_lethal = sum(final_results$gene == x & final_results$inpp4b_tam_t0_class == 'lethal')
	inpp4b_tam_t0_beneficial = sum(final_results$gene == x & final_results$inpp4b_tam_t0_class == 'beneficial')
	inpp4b_doxtam_t0_lethal = sum(final_results$gene == x & final_results$inpp4b_doxtam_t0_class == 'lethal')
	inpp4b_doxtam_t0_beneficial = sum(final_results$gene == x & final_results$inpp4b_doxtam_t0_class == 'beneficial')
	inpp4b_dox_unt_lethal = sum(final_results$gene == x & final_results$inpp4b_dox_unt_class == 'lethal')
	inpp4b_dox_unt_beneficial = sum(final_results$gene == x & final_results$inpp4b_dox_unt_class == 'beneficial')
	inpp4b_tam_unt_lethal = sum(final_results$gene == x & final_results$inpp4b_tam_unt_class == 'lethal')
	inpp4b_tam_unt_beneficial = sum(final_results$gene == x & final_results$inpp4b_tam_unt_class == 'beneficial')
	inpp4b_doxtam_dox_lethal = sum(final_results$gene == x & final_results$inpp4b_doxtam_dox_class == 'lethal')
	inpp4b_doxtam_dox_beneficial = sum(final_results$gene == x & final_results$inpp4b_doxtam_dox_class == 'beneficial')
	map3k1_veh_t0_lethal = sum(final_results$gene == x & final_results$map3k1_veh_t0_class == 'lethal')
	map3k1_veh_t0_beneficial = sum(final_results$gene == x & final_results$map3k1_veh_t0_class == 'beneficial')
	map3k1_dox_t0_lethal = sum(final_results$gene == x & final_results$map3k1_dox_t0_class == 'lethal')
	map3k1_dox_t0_beneficial = sum(final_results$gene == x & final_results$map3k1_dox_t0_class == 'beneficial')
	map3k1_tam_t0_lethal = sum(final_results$gene == x & final_results$map3k1_tam_t0_class == 'lethal')
	map3k1_tam_t0_beneficial = sum(final_results$gene == x & final_results$map3k1_tam_t0_class == 'beneficial')
	map3k1_doxtam_t0_lethal = sum(final_results$gene == x & final_results$map3k1_doxtam_t0_class == 'lethal')
	map3k1_doxtam_t0_beneficial = sum(final_results$gene == x & final_results$map3k1_doxtam_t0_class == 'beneficial')
	map3k1_dox_unt_lethal = sum(final_results$gene == x & final_results$map3k1_dox_unt_class == 'lethal')
	map3k1_dox_unt_beneficial = sum(final_results$gene == x & final_results$map3k1_dox_unt_class == 'beneficial')
	map3k1_tam_unt_lethal = sum(final_results$gene == x & final_results$map3k1_tam_unt_class == 'lethal')
	map3k1_tam_unt_beneficial = sum(final_results$gene == x & final_results$map3k1_tam_unt_class == 'beneficial')
	map3k1_doxtam_dox_lethal = sum(final_results$gene == x & final_results$map3k1_doxtam_dox_class == 'lethal')
	map3k1_doxtam_dox_beneficial = sum(final_results$gene == x & final_results$map3k1_doxtam_dox_class == 'beneficial')
	overall_veh_t0_lethal = sum(final_results$gene == x & final_results$overall_veh_t0_class == 'lethal')
	overall_veh_t0_beneficial = sum(final_results$gene == x & final_results$overall_veh_t0_class == 'beneficial')
	overall_tam_t0_lethal = sum(final_results$gene == x & final_results$overall_tam_t0_class == 'lethal')
	overall_tam_t0_beneficial = sum(final_results$gene == x & final_results$overall_tam_t0_class == 'beneficial')
	overall_tam_unt_lethal = sum(final_results$gene == x & final_results$overall_tam_unt_class == 'lethal')
	overall_tam_unt_beneficial = sum(final_results$gene == x & final_results$overall_tam_unt_class == 'beneficial')
	c(gata3_veh_t0_lethal,gata3_veh_t0_beneficial,gata3_dox_t0_lethal,gata3_dox_t0_beneficial,
	gata3_tam_t0_lethal,gata3_tam_t0_beneficial,gata3_doxtam_t0_lethal,gata3_doxtam_t0_beneficial,
	gata3_dox_unt_lethal,gata3_dox_unt_beneficial,gata3_tam_unt_lethal,gata3_tam_unt_beneficial,
	gata3_doxtam_dox_lethal,gata3_doxtam_dox_beneficial,ptpn12_veh_t0_lethal,ptpn12_veh_t0_beneficial,
	ptpn12_dox_t0_lethal,ptpn12_dox_t0_beneficial,ptpn12_tam_t0_lethal,ptpn12_tam_t0_beneficial,
	ptpn12_doxtam_t0_lethal,ptpn12_doxtam_t0_beneficial,ptpn12_dox_unt_lethal,ptpn12_dox_unt_beneficial,
	ptpn12_tam_unt_lethal,ptpn12_tam_unt_beneficial,ptpn12_doxtam_dox_lethal,ptpn12_doxtam_dox_beneficial,
	inpp4b_veh_t0_lethal,inpp4b_veh_t0_beneficial,inpp4b_dox_t0_lethal,inpp4b_dox_t0_beneficial,
	inpp4b_tam_t0_lethal,inpp4b_tam_t0_beneficial,inpp4b_doxtam_t0_lethal,inpp4b_doxtam_t0_beneficial,
	inpp4b_dox_unt_lethal,inpp4b_dox_unt_beneficial,inpp4b_tam_unt_lethal,inpp4b_tam_unt_beneficial,
	inpp4b_doxtam_dox_lethal,inpp4b_doxtam_dox_beneficial,map3k1_veh_t0_lethal,map3k1_veh_t0_beneficial,
	map3k1_dox_t0_lethal,map3k1_dox_t0_beneficial,map3k1_tam_t0_lethal,map3k1_tam_t0_beneficial,
	map3k1_doxtam_t0_lethal,map3k1_doxtam_t0_beneficial,map3k1_dox_unt_lethal,map3k1_dox_unt_beneficial,
	map3k1_tam_unt_lethal,map3k1_tam_unt_beneficial,map3k1_doxtam_dox_lethal,map3k1_doxtam_dox_beneficial,
	overall_veh_t0_lethal,overall_veh_t0_beneficial,
	overall_tam_t0_lethal,overall_tam_t0_beneficial,
	overall_tam_unt_lethal,overall_tam_unt_beneficial
	)
}))
colnames(class_count_table) = c(
	'gata3_veh_t0_lethal','gata3_veh_t0_beneficial','gata3_dox_t0_lethal','gata3_dox_t0_beneficial',
	'gata3_tam_t0_lethal','gata3_tam_t0_beneficial','gata3_doxtam_t0_lethal','gata3_doxtam_t0_beneficial',
	'gata3_dox_unt_lethal','gata3_dox_unt_beneficial','gata3_tam_unt_lethal','gata3_tam_unt_beneficial',
	'gata3_doxtam_dox_lethal','gata3_doxtam_dox_beneficial','ptpn12_veh_t0_lethal','ptpn12_veh_t0_beneficial',
	'ptpn12_dox_t0_lethal','ptpn12_dox_t0_beneficial','ptpn12_tam_t0_lethal','ptpn12_tam_t0_beneficial',
	'ptpn12_doxtam_t0_lethal','ptpn12_doxtam_t0_beneficial','ptpn12_dox_unt_lethal','ptpn12_dox_unt_beneficial',
	'ptpn12_tam_unt_lethal','ptpn12_tam_unt_beneficial','ptpn12_doxtam_dox_lethal','ptpn12_doxtam_dox_beneficial',
	'inpp4b_veh_t0_lethal','inpp4b_veh_t0_beneficial','inpp4b_dox_t0_lethal','inpp4b_dox_t0_beneficial',
	'inpp4b_tam_t0_lethal','inpp4b_tam_t0_beneficial','inpp4b_doxtam_t0_lethal','inpp4b_doxtam_t0_beneficial',
	'inpp4b_dox_unt_lethal','inpp4b_dox_unt_beneficial','inpp4b_tam_unt_lethal','inpp4b_tam_unt_beneficial',
	'inpp4b_doxtam_dox_lethal','inpp4b_doxtam_dox_beneficial','map3k1_veh_t0_lethal','map3k1_veh_t0_beneficial',
	'map3k1_dox_t0_lethal','map3k1_dox_t0_beneficial','map3k1_tam_t0_lethal','map3k1_tam_t0_beneficial',
	'map3k1_doxtam_t0_lethal','map3k1_doxtam_t0_beneficial','map3k1_dox_unt_lethal','map3k1_dox_unt_beneficial',
	'map3k1_tam_unt_lethal','map3k1_tam_unt_beneficial','map3k1_doxtam_dox_lethal','map3k1_doxtam_dox_beneficial',
	'overall_veh_t0_lethal','overall_veh_t0_beneficial',
	'overall_tam_t0_lethal','overall_tam_t0_beneficial',
	'overall_tam_unt_lethal','overall_tam_unt_beneficial'
)

final_table = cbind(final_results, class_count_table)	

final_table$pool = sapply(all_final_table$shRNA, function(x) {
	strsplit(as.character(x), '_', fixed = T)[[1]][2]
})
final_table$gene_pool = sapply(1:nrow(all_final_table), function(x) {
	paste(all_final_table$gene[x], all_final_table$pool[x], sep = '_')
})

write.table(final_table, paste(l_thr[1], p_thr[1], 'deseq_table.txt', sep = '_'), sep = '\t', row.names = F)

###################################

gene_table = data.frame(gene_pool = unique(final_table$gene_pool))

dat = subset(final_table, overall_veh_t0_class != 'neutral')
gene_table$gata3_veh_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_veh_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$ptpn12_veh_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_veh_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$inpp4b_veh_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_veh_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$map3k1_veh_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_veh_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$overall_veh_t0_l2fc = sapply(1:nrow(gene_table), function(x) {
	sum(gene_table$gata3_veh_t0_l2fc[x], gene_table$ptpn12_veh_t0_l2fc[x], gene_table$inpp4b_veh_t0_l2fc[x], gene_table$map3k1_veh_t0_l2fc[x]) 
})
dat = subset(final_table, overall_tam_t0_class != 'neutral')
gene_table$gata3_tam_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_tam_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$ptpn12_tam_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_tam_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$inpp4b_tam_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_tam_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$map3k1_tam_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_tam_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$overall_tam_t0_l2fc = sapply(1:nrow(gene_table), function(x) {
	sum(gene_table$gata3_tam_t0_l2fc[x], gene_table$ptpn12_tam_t0_l2fc[x], gene_table$inpp4b_tam_t0_l2fc[x], gene_table$map3k1_tam_t0_l2fc[x]) 
})
dat = subset(final_table, overall_tam_unt_class != 'neutral')
gene_table$gata3_tam_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_tam_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$ptpn12_tam_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_tam_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$inpp4b_tam_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_tam_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$map3k1_tam_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_tam_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$overall_tam_unt_l2fc = sapply(1:nrow(gene_table), function(x) {
	sum(gene_table$gata3_tam_unt_l2fc[x], gene_table$ptpn12_tam_unt_l2fc[x], gene_table$inpp4b_tam_unt_l2fc[x], gene_table$map3k1_tam_unt_l2fc[x]) 
})
gene_table$gata3_doxtam_dox_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_doxtam_dox[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$ptpn12_doxtam_dox_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_doxtam_dox[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$inpp4b_doxtam_dox_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_doxtam_dox[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
gene_table$map3k1_doxtam_dox_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_doxtam_dox[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})

#################### grubbs #####################

gene_table$grubbs = sapply(1:nrow(gene_table), function(x) {
	effects = c(
		gene_table$gata3_tam_unt_l2fc[x],
		gene_table$ptpn12_tam_unt_l2fc[x],
		gene_table$inpp4b_tam_unt_l2fc[x],
		gene_table$map3k1_tam_unt_l2fc[x],
		gene_table$gata3_doxtam_dox_l2fc[x],
		gene_table$ptpn12_doxtam_dox_l2fc[x],
		gene_table$inpp4b_doxtam_dox_l2fc[x],
		gene_table$map3k1_doxtam_dox_l2fc[x]
	)
	grubbs.test(effects, two.sided = T)$p.value
})

#####################################################

dat = subset(final_table, gata3_dox_t0_class != 'neutral')
gene_table$gata3_dox_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_dox_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, ptpn12_dox_t0_class != 'neutral')
gene_table$ptpn12_dox_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_dox_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, inpp4b_dox_t0_class != 'neutral')
gene_table$inpp4b_dox_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_dox_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, map3k1_dox_t0_class != 'neutral')
gene_table$map3k1_dox_t0_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_dox_t0[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, gata3_dox_unt_class != 'neutral')
gene_table$gata3_dox_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_gata3_dox_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, ptpn12_dox_unt_class != 'neutral')
gene_table$ptpn12_dox_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_ptpn12_dox_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, inpp4b_dox_unt_class != 'neutral')
gene_table$inpp4b_dox_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_inpp4b_dox_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})
dat = subset(final_table, map3k1_dox_unt_class != 'neutral')
gene_table$map3k1_dox_unt_l2fc = sapply(gene_table$gene_pool, function(x) {
	effects = dat$l2fc_map3k1_dox_unt[dat$gene_pool == x]
	neg_effects = effects[effects < 0]
	if(length(neg_effects) > 0) {
	} else {
		neg_effects = 0
	}
	pos_effects = effects[effects >= 0]
	if(length(pos_effects) > 0) {
	} else {
		pos_effects = 0
	}
	l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
	l2fc
})


###################################



sig_results = lapply(results, function(x) x[!is.na(x$log2FoldChange),])
sig_results = lapply(sig_results, function(x) x[!is.na(x$padj),])

sig_results = lapply(seq_along(sig_results), function(x) sig_results[[x]][abs(sig_results[[x]]$log2FoldChange) > l_thr[x],])
sig_results = lapply(seq_along(sig_results), function(x) sig_results[[x]][sig_results[[x]]$padj < p_thr[x],])

shrna = lapply(seq_along(sig_results), function(x) {
	sapply(rownames(sig_results[[x]]), function(y) {
		strsplit(as.character(y), '-', fixed = T)[[1]][1]
	})
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], shrna = shrna[[x]]))
gene = lapply(seq_along(sig_results), function(x) {
	sapply(rownames(sig_results[[x]]), function(y) {
		strsplit(as.character(y), '-', fixed = T)[[1]][2]
	})
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], gene = gene[[x]]))
gene_l2fc = lapply(seq_along(sig_results), function(x) {
	sapply(sig_results[[x]]$gene, function(y) {
		effects = sig_results[[x]]$log2FoldChange[sig_results[[x]]$gene == y]
		neg_effects = effects[effects < 0]
		if(length(neg_effects) > 0) {
		} else {
			neg_effects = 0
		}
		pos_effects = effects[effects >= 0]
		if(length(pos_effects) > 0) {
		} else {
			pos_effects = 0
		}
		l2fc  = sum(mean(neg_effects) * length(neg_effects), mean(pos_effects) * length(pos_effects))
		if(l2fc > 10) {10} else if(l2fc < -10) {-10} else {l2fc}
	})
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], gene_l2fc = gene_l2fc[[x]]))
gene_pval = lapply(seq_along(sig_results), function(i) {
	mapply(function(j,k) {
		if(j >= 0) {
			pvals = sig_results[[i]]$padj[sig_results[[i]]$log2FoldChange >= 0 & sig_results[[i]]$gene == k]
			Xsq = -2 * sum(log(pvals))
			p = pchisq(Xsq, df = 2 * length(pvals), lower.tail = F)
			if(p < 1e-20) {1e-20} else {p}
		} else {		
			pvals = sig_results[[i]]$padj[sig_results[[i]]$log2FoldChange < 0 & sig_results[[i]]$gene == k]
			Xsq = -2 * sum(log(pvals))
			p = pchisq(Xsq, df = 2 * length(pvals), lower.tail = F)
			if(p < 1e-20) {1e-20} else {p}
		}
	}, sig_results[[i]]$gene_l2fc, sig_results[[i]]$gene)
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], gene_pval = gene_pval[[x]]))
hit_hps = lapply(seq_along(sig_results), function(i) {
	mapply(function(j,k) {
		if(j > 0) {
			count = sum(sig_results[[i]]$log2FoldChange >= 0 & sig_results[[i]]$gene == k)
		} else {
			count = sum(sig_results[[i]]$log2FoldChange < 0 & sig_results[[i]]$gene == k)
		}
	}, sig_results[[i]]$gene_l2fc, sig_results[[i]]$gene)
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], hit_hps = hit_hps[[x]]))
gene_res = lapply(seq_along(sig_results), function(i) {
	mapply(function(j,k) {
		res = sqrt(j^2 + (-log(k,10))^2) * sign(j)
		all_l2fc = sig_results[[i]]$gene_l2fc
		all_pvals = -log(sig_results[[i]]$gene_pval, 10)
		sign_indicator = sign(all_l2fc)
		all_res = sqrt(all_l2fc^2 + (all_pvals - 0)^2) * sign_indicator
		sd_res = sd(all_res)
		res / sd_res
	}, sig_results[[i]]$gene_l2fc, sig_results[[i]]$gene_pval)
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], gene_res = gene_res[[x]]))
gene_sig = lapply(seq_along(sig_results), function(i) {
	mapply(function(k) {
		ifelse(k >= 4, 'Yes', 'No')
	}, sig_results[[i]]$hit_hps)
})
sig_results = lapply(seq_along(sig_results), function(x) data.frame(sig_results[[x]], gene_sig = gene_sig[[x]]))
names(sig_results) = names(results)

draw_plot = function(data, name) {
	tiff(filename = paste(name, '_volcano_plot.png', sep = ''), res = 300, width = 3000, height = 3000)
	print(
		ggplot(data, aes(gene_l2fc, -log(gene_pval, 10), size = factor(hit_hps))) +
		geom_point(data = subset(data, gene_sig == 'No'), shape = 21, fill = 'grey', color = "black", stroke = 1, alpha = 0.75) +
		geom_point(data = subset(data, gene_sig == 'Yes'), shape = 21, fill = 'red', color = "black", stroke = 1, alpha = 0.75) +
		geom_label_repel(data = subset(data, gene_sig == 'Yes'), aes(label = gene), 
			fill = 'white', color = 'black', segment.size = .5, segment.color = 'black', size = 3, force = 5, point.padding = unit(1, "lines")) +
		scale_size_manual(
			breaks = 1:10,
			values = c(1, 2*2:10)
		) +
		geom_hline(yintercept = 2, color = 'black', linetype = 2, alpha = 0.5) +
		geom_vline(xintercept = 1, color = 'black', linetype = 2, alpha = 0.5) +
		geom_vline(xintercept = -1, color = 'black', linetype = 2, alpha = 0.5) +
		scale_y_continuous(expand = c(.1,0)) + scale_x_continuous(expand = c(.1,0), limits = c(-10,10)) +
		labs(title = paste(name, 'Volcano Plot', sep = ' '), x = 'l2FC Sum', y = '-log(Fisher P-value)') +
		theme_bw()
	)
	dev.off()
}
volcano_data = lapply(sig_results, function(x) unique(x[,c(8:13)]))
names(volcano_data) = names(sig_results)
invisible(lapply(seq_along(volcano_data), function(x) write.table(volcano_data[[x]], paste(names(volcano_data)[x], 'volcano_plot_data.txt', sep = '_'), row.names = F, sep = '\t')))
invisible(lapply(seq_along(volcano_data), function(x) draw_plot(data = volcano_data[[x]], name = names(volcano_data)[x])))

############ finding background suppressed/enhanced modifiers ################

library(outliers)
all_final_table$pool = sapply(all_final_table$shRNA, function(x) {
	strsplit(as.character(x), '_', fixed = T)[[1]][2]
})
all_final_table$gene_pool = sapply(1:nrow(all_final_table), function(x) {
	paste(all_final_table$gene[x], all_final_table$pool[x], sep = '_')
})


