#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")
library(DESeq2)
#install.packeges("tidyverse")
library(tidyverse)
#install.packages("ggplot2")
library(ggplot2)
library(ggrepel)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(ggtext)
library(pheatmap)

# rep1
setwd("--/countsdata/rep1")

# data preparation for DESeq
## bring in counts files
countfiles <- list.files(pattern = "htseq copy.tsv$")
col_name <- sub("_counts_htseq copy.tsv", "", countfiles)

## empty list to store dataframes
dfs <- list()

for (i in seq_along(countfiles)) {
  df <- read.delim(countfiles[i], 
                   sep = "\t",
                   header = FALSE, 
                   col.names = c("gene_ID", col_name[i]))
  dfs[[i]] <- df
}

# counts df
rep1_counts = Reduce(function(x, y) merge(x, y, by = "gene_ID", all = TRUE), dfs)

## replace AFUB with appropriate pattern
rep1_counts <- rep1_counts[grepl("^AFUB", rep1_counts$gene_ID), ]
row_names <- rep1_counts$gene_ID
rep1_counts <- rep1_counts[, -1]
rownames(rep1_counts) <- row_names

## reorder columns in dataframe
rep1_counts <- rep1_counts %>% select("AL.C1AR_S28_L002", "AL.C2AR_S29_L002", "AL.C3AR_S30_L002",
                                      "AL.01AR_S34_L002", "AL.02AR_S35_L002", "AL.03AR_S36_L002",
                                      "AL.C1FR_S31_L002", "AL.C2FR_S32_L002", "AL.C3FR_S33_L002",
                                      "AL.01FR_S37_L002", "AL.02FR_S38_L002", "AL.03FR_S39_L002")

## rename columns to match sample info
colnames(rep1_counts) <- c("Euploid_nodrug_1", "Euploid_nodrug_2", "Euploid_nodrug_3",
                            "Aneu7A_nodrug_1", "Aneu7A_nodrug_2", "Aneu7A_nodrug_3",
                            "Euploid_FK506_1", "Euploid_FK506_2", "Euploid_FK506_3",
                            "Aneu7A_FK506_1", "Aneu7A_FK506_2", "Aneu7A_FK506_3")

# filter out genes with less than 1 TPM in any sample (only keep those that have 1TPM or more average in at least 1 condition/strain)
## calculate gene lengths from gff
gff <- read.delim("A1163_CEA10T2T_liftoff.gff", header = FALSE, col.names = c("contig", "source", "cat", "start", "end", ".", "strand", ".", "ID"), sep = "\t", skip = 6)
gff_proteincoding <- gff[gff$cat == "protein_coding_gene", ]
gff_proteincoding$gene <- sub(".*(AFUB_\\d+).*", "\\1", gff_proteincoding$ID)
gene_lengths <- gff_proteincoding[, c("gene", "start", "end")]
gene_lengths$length <- gene_lengths$end - gene_lengths$start
row.names(gene_lengths) <- gene_lengths$gene
rep1_counts$gene <- row.names(rep1_counts)
tpm_counts_rep1 <- merge(rep1_counts, gene_lengths[ , c("gene", "length")], by = "gene")
## calculate tpms
### convert length from bp to kbp
tpm_counts_rep1$length <- tpm_counts_rep1$length/1000
### df/character vector to hold tpm values
tpm_calc_rep1 <- tpm_counts_rep1[, c("gene"), drop = FALSE]
### get cols except gene
tpm_cols_rep1 <- setdiff(colnames(tpm_counts_rep1), c("gene", "length"))
### get tpm
for (gene in tpm_cols_rep1) {
  rpk <- tpm_counts_rep1[[gene]] / tpm_counts_rep1$length
  tpm <- rpk / sum(rpk) * 1e6
  tpm_calc_rep1[[gene]] <- tpm
}
## average tpm for technical reps
tpm_calc_rep1 <- tpm_calc_rep1 %>%
  pivot_longer(
    cols = -gene,
    names_to = "replicate",
    values_to = "tpm"
  )

tpm_calc_rep1_long <- tpm_calc_rep1 %>%
  mutate(sample = str_extract(replicate, "^[^_]+_[^_]+"))

tpm_avg_rep1 <- tpm_calc_rep1_long %>%
  group_by(gene, sample) %>%
  summarise(avg = mean(tpm), .groups = "drop") %>%
  pivot_wider(names_from = "sample", values_from = "avg")

tpm_filtered_rep1 <- tpm_avg_rep1 %>%
  filter(if_any(-gene, ~ . >= 1)) %>%
  pull(gene)

counts_rep1_tpmfiltered <- rep1_counts[rownames(rep1_counts) %in% tpm_filtered_rep1, ]
counts_rep1_tpmfiltered <- counts_rep1_tpmfiltered[ ,1:12]

# make df with appropriate sample info
sample <- c("Euploid_nodrug_1", "Euploid_nodrug_2", "Euploid_nodrug_3", 
            "Aneu7_nodrug_1", "Aneu7_nodrug_2", "Aneu7_nodrug_3",
            "Euploid_FK506_1", "Euploid_FK506_2", "Euploid_FK506_3", 
            "Aneu7_FK506_1", "Aneu7_FK506_2", "Aneu7_FK506_3")

condition <- c("Euploid_nodrug", "Euploid_nodrug", "Euploid_nodrug", 
               "Aneu7_nodrug", "Aneu7_nodrug", "Aneu7_nodrug",
               "Euploid_FK506", "Euploid_FK506", "Euploid_FK506", 
               "Aneu7_FK506", "Aneu7_FK506", "Aneu7_FK506")

colData_rep1 <- data.frame(sample, condition)

rownames(colData_rep1) <- colnames(counts_rep1_tpmfiltered)

# DESeq
dds_rep1 <- DESeqDataSetFromMatrix(countData = counts_rep1_tpmfiltered,
                              colData = colData_rep1,
                              design = ~ condition) 
dds_rep1 <- DESeq(dds_rep1)

#### PCA
trans_rep1 <- vst(dds_rep1, blind = TRUE)
##### extract matrix and get variances
genevar_rep1 <- as.data.frame(sort(rowVars(assay(trans_rep1)), decreasing = TRUE))
colnames(genevar_rep1) <- "Variance"
genevar_rep1$Rank <- seq_len(nrow(genevar_rep1))

genevar_rep1_plot <- ggplot(genevar_rep1, aes(x = Rank, y = Variance)) +
  geom_point(shape = 1, size = 0.5) +
  geom_vline(xintercept = 1000, color = "red", linetype = "solid") +
  xlim(0,10000) +
  ylim(0,8) +
  theme_bw()
ggsave(file="--/genevar_rep1_plot.svg", plot=genevar_rep1_plot, width=3, height=2)

##### manually plot PCs to get variance in all components
genevar_rep1_top1000 <- rownames(genevar_rep1)[1:1000]
pca_rep1_manual_vst <- assay(trans_rep1)[genevar_rep1_top1000, ]
pca_rep1_manual_data <- prcomp(t(pca_rep1_manual_vst), scale = TRUE)
pca_rep1_variance <- pca_rep1_manual_data$sdev^2
pca_rep1_var_explained <- (pca_rep1_variance / sum(pca_rep1_variance))*100

pca_rep1_var_df <- as.data.frame(pca_rep1_var_explained)
pca_rep1_var_df$PC <- factor(paste0("PC", seq_along(pca_rep1_var_explained)), levels = paste0("PC", seq_along(pca_rep1_var_explained)))

rep1_pcs_plot <- ggplot(pca_rep1_var_df, aes(x = PC, y = pca_rep1_var_explained, fill = PC)) +
  geom_bar(stat = "identity", fill = "#004c4c") +
  theme_bw()
ggsave(file="--/rep1_pcs_plot.svg", plot=rep1_pcs_plot, width=3, height=2)

pca_rep1_df <- data.frame(
  Sample = rownames(pca_rep1_manual_data$x),
  PC1 = pca_rep1_manual_data$x[ , 1],
  PC2 = pca_rep1_manual_data$x[ , 2]       
)

pca_rep1_df <- cbind(pca_rep1_df, colData_rep1[rownames(pca_rep1_df), , drop = FALSE])

#pca_rep1_data <- plotPCA(trans_rep1, intgroup = c("condition"), returnData=TRUE) # PC1: 64%; PC2: 24%
pca_rep1 <- ggplot(pca_rep1_df, 
                           aes(PC1, PC2, color=condition, shape=condition)) + 
  geom_richtext(data = pca_rep1_data, 
                aes(label = condition),
                color = "black") +
  geom_point(size = 3, shape = c(23, 23, 23, 21, 21, 21, 23, 23, 23, 21, 21, 21),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#444545", "#bcbcbc", "#6850a1", "#b4a8d2")) +
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1), 
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 7)) +
  xlab("PC1 (52.2% variance)") +
  ylab("PC2 (21.6% variance)")
ggsave(file="--/pca_rep1.svg", plot=pca_rep1, width=2.5, height=2.5, units = "in")

# DEG comparisons
res_aneu_vs_euploid_FK506_rep1 <- as.data.frame(results(dds_rep1,
                                                     contrast = c("condition", "Aneu7_FK506",
                                                                  "Euploid_FK506"),
                                                     alpha = 0.01))

# save data
write.csv(res_aneu_vs_euploid_FK506_rep1, "res_aneu_vs_euploid_FK506_rep1.csv")

volcano_7A <- EnhancedVolcano(res_aneu_vs_euploid_FK506_rep1, 
                              lab = row.names(res_aneu_vs_euploid_FK506_rep1),
                              x = 'log2FoldChange',
                              y = 'padj',
                              title = NULL,
                              pCutoff = 0.05,
                              FCcutoff = 1,
                              xlab = 'log2FC',
                              ylab = expression(-log[10](padj)),
                              subtitle = NULL,
                              caption = NULL,
                              selectLab = c("AFUB_086660", "AFUB_086680", "AFUB_086700", "AFUB_086710", "AFUB_086720", "AFUB_086750"),
                              legendLabels = NULL,
                              legendPosition = NULL,
                              legendLabSize = 0,
                              legendIconSize = 0,
                              drawConnectors = TRUE,
                              titleLabSize = 8,
                              axisLabSize = 7,
                              pointSize = 1,
                              labSize = 1, 
                              #boxedLabels = TRUE,
                              colAlpha = 0.6,
                              col=c('#dad9d9', '#dad9d9', '#999999', '#684fa1'),
                              max.overlaps = 100,
                              xlim = c(-7, 7),
                              ylim = c(-0.5, 65),
                              gridlines.major = 0.3,
                              gridlines.minor = 0.1,
                              border = "partial",
                              borderWidth = 0.3,
                              borderColour = "black")

ggsave(file="--/volcano_rep1.svg", plot=volcano_7A, width=3.2, height=3)

# rep2
setwd("--/countsdata/rep2")

# data preparation for DESeq
## bring in counts files
countfiles_rep2<- list.files(pattern = "htseq.tsv$")
col_name_2 <- sub("_counts_htseq.tsv", "", countfiles_rep2)

## empty list to store dataframes
dfs_2 <- list()

for (i in seq_along(countfiles_rep2)) {
  df <- read.delim(countfiles_rep2[i], 
                   sep = "\t",
                   header = FALSE, 
                   col.names = c("gene_ID", col_name_2[i]))
  dfs_2[[i]] <- df
}

# counts df
rep2_counts = Reduce(function(x, y) merge(x, y, by = "gene_ID", all = TRUE), dfs_2)

## replace AFUB with appropriate pattern
rep2_counts <- rep2_counts[grepl("^AFUB", rep2_counts$gene_ID), ]
row_names2 <- rep2_counts$gene_ID
rep2_counts <- rep2_counts[, 2:13]
rownames(rep2_counts) <- row_names2

## reorder columns in dataframe
rep2_counts <- rep2_counts %>% select("AL_WT_A1_S387_L006", "AL_WT_A2_S388_L006", "AL_WT_A3_S389_L006",
                                      "AL_6BR_A1_S405_L006", "AL_6BR_A2_S406_L006", "AL_6BR_A3_S407_L006",
                                      "AL_WT_F1_S390_L006", "AL_WT_F2_S391_L006", "AL_WT_F3_S392_L006",
                                      "AL_6BR_F1_S408_L006", "AL_6BR_F2_S409_L006", "AL_6BR_F3_S410_L006")

## rename columns to match sample info
colnames(rep2_counts) <- c("Euploid2_nodrug_1", "Euploid2_nodrug_2", "Euploid2_nodrug_3",
                           "Aneu7B_nodrug_1", "Aneu7B_nodrug_2", "Aneu7B_nodrug_3",
                           "Euploid2_FK506_1", "Euploid2_FK506_2", "Euploid2_FK506_3",
                           "Aneu7B_FK506_1", "Aneu7B_FK506_2", "Aneu7B_FK506_3")

# filter out genes with less than 1 TPM in any sample (only keep those that have 1TPM or more average in at least 1 condition/strain)
## calculate gene lengths from gff
gff <- read.delim("A1163_CEA10T2T_liftoff.gff", header = FALSE, col.names = c("contig", "source", "cat", "start", "end", ".", "strand", ".", "ID"), sep = "\t", skip = 6)
gff_proteincoding <- gff[gff$cat == "protein_coding_gene", ]
gff_proteincoding$gene <- sub(".*(AFUB_\\d+).*", "\\1", gff_proteincoding$ID)
gene_lengths <- gff_proteincoding[, c("gene", "start", "end")]
gene_lengths$length <- gene_lengths$end - gene_lengths$start
row.names(gene_lengths) <- gene_lengths$gene
rep2_counts$gene <- row.names(rep2_counts)
tpm_counts_rep2 <- merge(rep2_counts, gene_lengths[ , c("gene", "length")], by = "gene")
## calculate tpms
### convert length from bp to kbp
tpm_counts_rep2$length <- tpm_counts_rep2$length/1000
### df/character vector to hold tpm values
tpm_calc_rep2 <- tpm_counts_rep2[, c("gene"), drop = FALSE]
### get cols except gene
tpm_cols_rep2 <- setdiff(colnames(tpm_counts_rep2), c("gene", "length"))
### get tpm
for (gene in tpm_cols_rep2) {
  rpk_2 <- tpm_counts_rep2[[gene]] / tpm_counts_rep2$length
  tpm_2 <- rpk_2 / sum(rpk_2) * 1e6
  tpm_calc_rep2[[gene]] <- tpm_2
}
## average tpm for technical reps
tpm_calc_rep2 <- tpm_calc_rep2 %>%
  pivot_longer(
    cols = -gene,
    names_to = "replicate",
    values_to = "tpm_2"
  )

tpm_calc_rep2_long <- tpm_calc_rep2 %>%
  mutate(sample = str_extract(replicate, "^[^_]+_[^_]+"))

tpm_avg_rep2 <- tpm_calc_rep2_long %>%
  group_by(gene, sample) %>%
  summarise(avg = mean(tpm_2), .groups = "drop") %>%
  pivot_wider(names_from = "sample", values_from = "avg")

tpm_filtered_rep2 <- tpm_avg_rep2 %>%
  filter(if_any(-gene, ~ . >= 1)) %>%
  pull(gene)

counts_rep2_tpmfiltered <- rep2_counts[rownames(rep2_counts) %in% tpm_filtered_rep2, ]
counts_rep2_tpmfiltered <- counts_rep2_tpmfiltered[ ,1:12]

# make df with appropriate sample info
sample <- c("Euploid_nodrug_1", "Euploid_nodrug_2", "Euploid_nodrug_3", 
            "Aneu7_nodrug_1", "Aneu7_nodrug_2", "Aneu7_nodrug_3",
            "Euploid_FK506_1", "Euploid_FK506_2", "Euploid_FK506_3", 
            "Aneu7_FK506_1", "Aneu7_FK506_2", "Aneu7_FK506_3")

condition <- c("Euploid_nodrug", "Euploid_nodrug", "Euploid_nodrug", 
               "Aneu7_nodrug", "Aneu7_nodrug", "Aneu7_nodrug",
               "Euploid_FK506", "Euploid_FK506", "Euploid_FK506", 
               "Aneu7_FK506", "Aneu7_FK506", "Aneu7_FK506")

colData_rep2 <- data.frame(sample, condition)

rownames(colData_rep2) <- colnames(counts_rep2_tpmfiltered)

# DESeq
dds_rep2 <- DESeqDataSetFromMatrix(countData = counts_rep2_tpmfiltered,
                                   colData = colData_rep2,
                                   design = ~ condition) 
dds_rep2 <- DESeq(dds_rep2)

#### PCA
trans_rep2 <- vst(dds_rep2, blind = TRUE)
##### extract matrix and get variances
genevar_rep2 <- as.data.frame(sort(rowVars(assay(trans_rep2)), decreasing = TRUE))
colnames(genevar_rep2) <- "Variance"
genevar_rep2$Rank <- seq_len(nrow(genevar_rep2))

genevar_rep2_plot <- ggplot(genevar_rep2, aes(x = Rank, y = Variance)) +
  geom_point(shape = 1, size = 0.5) +
  geom_vline(xintercept = 1000, color = "red", linetype = "solid") +
  xlim(0,10000) +
  ylim(0, 8) +
  theme_bw()
ggsave(file="--/genevar_rep2_plot.svg", plot=genevar_rep2_plot, width=3, height=2)


##### manually plot PCs to get variance in all components
genevar_rep2_top1000 <- rownames(genevar_rep2)[1:1000]
pca_rep2_manual_vst <- assay(trans_rep2)[genevar_rep2_top1000, ]
pca_rep2_manual_data <- prcomp(t(pca_rep2_manual_vst), scale = TRUE)
pca_rep2_variance <- pca_rep2_manual_data$sdev^2
pca_rep2_var_explained <- (pca_rep2_variance / sum(pca_rep2_variance))*100

pca_rep2_var_df <- as.data.frame(pca_rep2_var_explained)
pca_rep2_var_df$PC <- factor(paste0("PC", seq_along(pca_rep2_var_explained)), levels = paste0("PC", seq_along(pca_rep2_var_explained)))

rep2_pcs_plot <- ggplot(pca_rep2_var_df, aes(x = PC, y = pca_rep2_var_explained, fill = PC)) +
  geom_bar(stat = "identity", fill = "#004c4c") +
  theme_bw()
ggsave(file="--/rep2_pcs_plot.svg", plot=rep2_pcs_plot, width=3, height=2)

pca_rep2_df <- data.frame(
  Sample = rownames(pca_rep2_manual_data$x),
  PC1 = pca_rep2_manual_data$x[ , 1],
  PC2 = pca_rep2_manual_data$x[ , 2]       
)

pca_rep2_df <- cbind(pca_rep2_df, colData_rep2[rownames(pca_rep2_df), , drop = FALSE])

#pca_rep2_data <- plotPCA(trans_rep2, intgroup = c("condition"), returnData=TRUE) # PC1: 65%; PC2: 25%
pca_rep2 <- ggplot(pca_rep2_df, 
                   aes(PC1, PC2, color=condition, shape=condition)) + 
  geom_richtext(data = pca_rep2_data, 
                aes(label = condition),
                color = "black") +
  geom_point(size = 3, shape = c( 21, 21, 21, 24, 24, 24, 21, 21, 21, 24, 24, 24),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#397739","#b6d7a8", "#444545", "#bcbcbc")) +
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1), 
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 7)) +
  xlab("PC1 (60.5% variance)") +
  ylab("PC2 (19.3% variance")
ggsave(file="pca_rep2.svg", plot=pca_rep2, width=3, height=3)

ggsave(file="--/pca_rep2.svg", plot=pca_rep2, width=2.5, height=2.5)

# DEG comparisons
res_aneu_vs_euploid_FK506_rep2 <- as.data.frame(results(dds_rep2,
                                                        contrast = c("condition", "Aneu7_FK506",
                                                                     "Euploid_FK506"),
                                                        alpha = 0.01))
# save data
write.csv(res_aneu_vs_euploid_FK506_rep2, "res_aneu_vs_euploid_FK506_rep2.csv")

## get results for aneu vs euploid no drug comparisons
res_aneu_vs_euploid_nodrug_rep1 <- as.data.frame(results(dds_rep1,
                                                        contrast = c("condition", "Aneu7_nodrug",
                                                                     "Euploid_nodrug"),
                                                        alpha = 0.01))
write.csv(res_aneu_vs_euploid_nodrug_rep1, "res_aneu_vs_euploid_nodrug_rep1.csv")

res_aneu_vs_euploid_nodrug_rep2 <- as.data.frame(results(dds_rep2,
                                                         contrast = c("condition", "Aneu7_nodrug",
                                                                      "Euploid_nodrug"),
                                                         alpha = 0.01))
write.csv(res_aneu_vs_euploid_nodrug_rep2, "res_aneu_vs_euploid_nodrug_rep2.csv")

volcano_7B <- EnhancedVolcano(res_aneu_vs_euploid_FK506_rep2, 
                              lab = row.names(res_aneu_vs_euploid_FK506_rep2),
                              x = 'log2FoldChange',
                              y = 'padj',
                              title = NULL,
                              pCutoff = 0.05,
                              FCcutoff = 1,
                              xlab = 'log2FC',
                              ylab = expression(-log[10](padj)),
                              subtitle = NULL,
                              caption = NULL,
                              selectLab = c("AFUB_086660", "AFUB_086680", "AFUB_086700", "AFUB_086710", "AFUB_086720", "AFUB_086750"),
                              legendLabels = NULL,
                              legendPosition = NULL,
                              legendLabSize = 0,
                              legendIconSize = 0,
                              drawConnectors = TRUE,
                              titleLabSize = 8,
                              axisLabSize = 7,
                              pointSize = 1,
                              labSize = 1, 
                              #boxedLabels = TRUE,
                              colAlpha = 0.6,
                              col=c('#dad9d9', '#dad9d9', '#999999', '#397739'),
                              max.overlaps = 200,
                              xlim = c(-7, 7),
                              ylim = c(-0.5, 65),
                              gridlines.major = 0.3,
                              gridlines.minor = 0.1,
                              border = "partial",
                              borderWidth = 0.3,
                              borderColour = "black")

ggsave(file="--/volcano_rep2.svg", plot=volcano_7B, width=3.2, height=3)

# make heatmap of tpms in all BGCs (Fig. S4)
SM_genes <- c(#DHN melanin
  "AFUB_033230", "AFUB_033220", "AFUB_033240", "AFUB_033270", "AFUB_033250", "AFUB_033290",
  # gliotoxin
  "AFUB_075680", "AFUB_075690", "AFUB_075700", "AFUB_075710", "AFUB_075720", "AFUB_075730", "AFUB_075740", "AFUB_075750", "AFUB_075760", "AFUB_075770", "AFUB_075780", "AFUB_075790", 
  "AFUB_026880",
  # fumigaclanive
  "AFUB_033750", "AFUB_033740", "AFUB_033730", "AFUB_033720", "AFUB_033710", "AFUB_033700", "AFUB_033690", "AFUB_033680", "AFUB_033670", "AFUB_033660", "AFUB_033650",
  #fumitremorgin
  "AFUB_086360", "AFUB_086350", "AFUB_086340", "AFUB_086330", "AFUB_086320", "AFUB_086310", "AFUB_086300", "AFUB_086290", "AFUB_086280",
  # pseurotin
  "AFUB_086150", "AFUB_086130", "AFUB_086040", "AFUB_086030", "AFUB_086020", "AFUB_086010", "AFUB_085990",
  # helvolic acid
  "AFUB_072030", "AFUB_072040", "AFUB_072050", "AFUB_072060", "AFUB_072070", "AFUB_072080", "AFUB_072090", "AFUB_072100", "AFUB_072110",
  # fumiquinazoline
  "AFUB_078070", "AFUB_078050", "AFUB_078040", "AFUB_078060", "AFUB_078030",
  # pyripyropene
  "AFUB_000830", "AFUB_000820", "AFUB_000810", "AFUB_000800", "AFUB_000790", "AFUB_000780", "AFUB_000770", "AFUB_000760", "AFUB_000750",
  # endocrocin
  "AFUB_100730", "AFUB_100740", "AFUB_100750", "AFUB_100760",
  # fumagillin
  "AFUB_086200", "AFUB_086190", "AFUB_086180", "AFUB_086150", "AFUB_086100", "AFUB_086090", "AFUB_086060", "AFUB_086050",
  # hexadehdroastechrome
  "AFUB_036300", "AFUB_036290", "AFUB_036280", "AFUB_036270", "AFUB_036260", "AFUB_036250", "AFUB_036240", "AFUB_036230",
  # neosartoricin
  "AFUB_086720", "AFUB_086710", "AFUB_086700", "AFUB_086690", "AFUB_086680", "AFUB_086670",
  # trypacidin
  "AFUB_071820", "AFUB_071810", "AFUB_071800", "AFUB_071790", "AFUB_071780", "AFUB_071770", "AFUB_071760", "AFUB_071750", "AFUB_071740", "AFUB_071730", 
  "AFUB_071720", "AFUB_071710",
  # fumisoquin
  "AFUB_094860", "AFUB_094850", "AFUB_094840", "AFUB_094830", "AFUB_094820", "AFUB_094810", "AFUB_094800",
  # xanthocillin
  "AFUB_051210", "AFUB_051200", "AFUB_051190", "AFUB_051180", "AFUB_051170", "AFUB_051160", "AFUB_051150",
  # fumihopaside
  "AFUB_071550", "AFUB_071560", "AFUB_071569", "AFUB_071570",
  # fumivaline
  "AFUB_035500", "AFUB_035510", "AFUB_035520", "AFUB_035530",
  # sphingofungin
  "AFUB_034530", "AFUB_034520", "AFUB_034510", "AFUB_034500", "AFUB_034490", "AFUB_034480", "AFUB_034470", "AFUB_034460", "AFUB_034450",
  # satorypyrone
  "AFUB_084240", "AFUB_084230", "AFUB_084220", "AFUB_084210", "AFUB_084200", "0841090")

tpms_merged <- merge(tpm_avg_rep1, tpm_avg_rep2, by = "gene")

SM_tpms <- subset(tpms_merged, tpms_merged$gene %in% SM_genes)
rownames(SM_tpms) <- SM_tpms$gene
SM_tpms <- SM_tpms[ , 2:9]
SM_tpms <- SM_tpms %>% select(c("Euploid_nodrug", "Euploid2_nodrug", "Aneu7A_nodrug", "Aneu7B_nodrug", 
                                "Euploid_FK506", "Euploid2_FK506", "Aneu7A_FK506", "Aneu7B_FK506"))
SM_tpms <- SM_tpms[SM_genes, ]

breaks <- seq(0, 500, length.out = 101)
colors_heatmap <- colorRampPalette(c("white", "#f4cccc", "#ea9999", "#e06666", "#cc0000", "#990000", "#660000"))(length(breaks) - 1)

S4_heatmap <- pheatmap(SM_tpms,
color = colors_heatmap,
breaks = breaks,
cluster_rows = FALSE,
cluster_cols = FALSE,
fontsize = 5)

ggsave(file="--/S4_heatmap.svg", plot=S4_heatmap, width=6, height=9, units = "in")
