library(tximport)
library(ensembldb)
library(AnnotationHub)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(genefilter)
library(org.Hs.eg.db)
BiocManager::install("M3C")
library('M3C')

hub = AnnotationHub()

#make sure to use the right species
ensdb_query <- query(hub, c("EnsDb", "sapiens", "110"))
ensdb_query

ensdb_109 <- ensdb_query[['AH113665']]

# Extract transcript and gene information
tx_data <- transcripts(ensdb_109, return.type = "DataFrame")

# Create the tx2gene data.frame
tx2gene <- tx_data[, c("tx_id", "gene_id")]

# Especifique o caminho absoluto para o diretório
dir_path <- "/home/bvmh/Salmon_Trusight/sf"

# Lista todos os arquivos dentro do diretório com extensão .sf
files <- list.files(dir_path, pattern = "\\.sf$", full.names = TRUE)
sample_names <- gsub("_quant$", "", basename(files))
names(files) <- sample_names
files

# Importa os arquivos usando tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

df_tx <- as.data.frame(txi)

tx_abundance <- data.frame(
  GeneName = row.names(df_tx),
  '4833BM2' = df_tx$abundance.4833BM2.sf,
  'BM1014-1' = df_tx$abundance.BM1284.1.sf,
  'BM1126-1' = df_tx$abundance.BM1126.1.sf,
  'BM1284-5' = df_tx$abundance.BM1284.5.sf,
  'BM374-4' = df_tx$abundance.BM374.4.sf,
  'BM943-1' = df_tx$abundance.BM943.1.sf,
  '4959BM1' = df_tx$abundance.4959BM1.sf,
  'BM1072-1' = df_tx$abundance.BM1072.1.sf,
  'BM1284-1' = df_tx$abundance.BM1284.1.sf,
  'BM374-1' = df_tx$abundance.BM374.1.sf, 
  'BM911-1' = df_tx$abundance.BM911.1.sf,
  'BM943-4' = df_tx$abundance.BM943.4.sf
)

tx_abundance$symbol <- mapIds(org.Hs.eg.db, keys = tx_abundance$GeneName, keytype = 'ENSEMBL', column = 'SYMBOL')

# Carregar os símbolos do arquivo de texto para uma lista
symbols <- readLines("/home/bvmh/trusight-rna-fusion-gene-list.txt")

# Filtrar o dataframe usando os símbolos da lista
filtered_df <- tx_abundance %>%
  filter(symbol %in% symbols)

tx_abundance_unique <- filtered_df[!duplicated(filtered_df$symbol), ]

row.names(tx_abundance_unique) <- tx_abundance_unique$symbol

tx_abundance_unique <- tx_abundance_unique[ , -ncol(tx_abundance_unique)]

tx_abundance_unique <- tx_abundance_unique[, -1]
tx_abundance_unique <- tx_abundance_unique[, -45]

tx_abundance_unique_log2 <- log2(tx_abundance_unique)

write.csv(tx_abundance_unique_log2, file = "tx_abundance.csv")

tx_abundance_unique_log2
