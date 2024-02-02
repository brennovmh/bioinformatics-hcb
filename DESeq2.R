library(DESeq2)
library(ggplot2)
library(ComplexHeatMap)
library(circlize)
library(tibble)


count -> Counts
Counts <- Counts[, c("Geneid", "X.BM1075.1", "X.BM1077.1", "X.BM1078.1", "X.BM1079.1", "X.BM1092.1", "X.BM1115.1", "X.BM1116.1", "X.BM1141.1", "X.BM1146.1", "X.BM1152.1", "X.BM1153.1", "X.BM1158.1", "X.BM1161.1", "X.BM1170.2", "X.BM1179.1", "X.BM1181.1", "X.BM1182.1", "X.BM1185.2", "X.BM1186.1", "X.BM1189.1", "X.BM1206.1", "X.BM1260.1", "X.BM1277.1", "X.BM1284.1", "X.BM1284.5", "X.BM1292.1", "X.BM1314.1", "X.BM1324.1", "X.BM1395.1", "X.BM1471.1", "X.BM1489.1", "X.BM1490.1", "X.BM1498.1", "X.BM1500.1", "X.BM1501.1", "X.BM1503.1", "X.BM1505.1", "X.BM1507.1", "X.BM1519.1", "X.BM1523.1", "X.BM1526.1", "X.BM1532.1", "X.BM770.13", "X.BM810.1", "X.BM821.4", "X.BM943.4", "X.BM999.2"
)]


#Transformar primeira coluna em index
Counts %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Geneid') -> Counts

## FIltrar linhas com count menor que 50
Counts <- Counts[rowSums(Counts) >= 50, ]


condition <- factor(c('1and19', 'C', 'C', 'C', 'C', '1and19', 'C', '1and19', 'C', 'C', 'C', 'C', '12and21', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', '12and21', '12and21', 'C', 'C', 'C', 'C', 'C', '12and21', 'IKZ', '12and21', 'C', '12and21', 'C', '12and21', 'BCRABL', '12and21', 'BRD9andNUT', '12and21', 'C', '12and21', 'IGH', 'BCRABL', '12and21', 'C'))

coldata <- data.frame(row.names = colnames(Counts), condition)

dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)

dds <- DESeq(dds)

vsdata <- vst(dds, blind = FALSE)

plotPCA(vsdata, intgroup = 'condition')

plotDispEsts(dds)

res <- results(dds, contrast = c('condition', '12and21', 'C'))

# Significant genes
sigs <- na.omit(res)
sigs <- sigs[sigs$padj<0.05,]
sigs

library('org.Hs.eg.db')

sigs_df <- as.data.frame(sigs)
sigs_df$symbol <- mapIds(org.Hs.eg.db, keys = row.names(sigs_df), keytype = 'ENSEMBL', column = 'SYMBOL')
sigs_df

df.top <- sigs_df[ sigs_df$baseMean >30 & (abs(sigs_df$log2FoldChange) >2),]
df.top

df.top <- df.top[order(df.top$log2FoldChange, decreasing = TRUE),]
df.top

rlog_out <- rlog(dds, blind =FALSE) #get normalized count data from dds object mat<-assay(rlog_out)[rownanes(df.top), rownames(coldata)] #sig genes × samples
mat <- assay(rlog_out)[row.names(df.top), row.names(coldata)]
colnames (mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, function(x) scale(x))) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(mat)

num_keep <- 10
rows_keep <- c(seq(1:num_keep), seq(nrow(mat.scaled)-num_keep), nrow(mat.scaled))

l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"

mean <- as.matrix(df.top[rows_keep,]$baseMean)
colnames(mean) <- "AveExpr"

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red"))
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2), 
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
                column_labels = colnames(mat.scaled), name="Z-score",
                cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = df.top$symbol[rows_keep], 
                cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                  grid.text(round(l2_val[i, j],2), x, y)
                })
h3 <- Heatmap(mean, row_labels = df.top$symbol[rows_keep], 
                cluster_rows = F, name = "AveExpr", col=col_AveExpr,
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                  grid.text(round(mean[i, j],2), x, y)
                })

h<-h1+h2+h3
h
