#!/usr/bin/env Rscript

###############################################
### STEP1. Construct rsem to edgeR input
###############################################
BiocManager::install("tximportData")
library(tximport)
path.case = "/home/hykim/RNA-seq/"

files_N <- list.files(path = path.case, pattern = "N.genes.results$", full.names = TRUE)
files_T <- list.files(path = path.case, pattern = "T.genes.results$", full.names = TRUE)
files <- c(files_N, files_T)
print(length(files))

N.list<-dir(path=path.case, pattern = "N.genes.results")
sample.N.list <- list()
for( i in 1:length(N.list)){
        sample.N.list[i] <- unlist(strsplit(N.list[i],split=".genes.results"))[1]
        }
T.list<-dir(path=path.case, pattern = "T.genes.results")
sample.T.list <- list()
for( i in 1:length(T.list)){
        sample.T.list[i] <- unlist(strsplit(T.list[i],split=".genes.results"))[1]
        }
sample.list <- c(sample.N.list, sample.T.list)
names(files) <- sample.list

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
a <- round(txi.rsem$counts[,], digit=0)
write.table(a, "edgeR_input.txt",col.names=NA, row.names=T, quote=F,sep='\t')
##write.table(txi.rsem$abundance, "edgeR_input.txt",col.names=NA, row.names=T, quote=F,sep='\t')

###############################################
### STEP2. Convert Ensembl to Entrezgene id
###############################################
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library('biomaRt')

df <- read.table("edgeR_input.txt",sep='\t',header=T,stringsAsFactors=F)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

for (i in 1:length(df$X)){     
	df$X[i] <- unlist(strsplit(df$X[i], fixed=TRUE, split="."))[1]
}

genes <- df$X
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","entrezgene_id", "description"), values=genes, mart= mart)
df$entrezgene_id=""
df["entrezgene_id"] = lapply("entrezgene_id", function(x) G_list[[x]][match(df$X, G_list$ensembl_gene_id)])

#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description"),values=genes,mart= mart)
#df$hgnc_symbol = "" 
#df["hgnc_symbol"] = lapply("hgnc_symbol", function(x) G_list[[x]][match(df$X, G_list$ensembl_gene_id)]) 

write.table(df, "edgeR_input_ENTREZ.txt", col.names=NA, row.names=T, quote=F,sep='\t')

###############################################
### STEP3. CMScaller
###############################################
library(CMScaller)
df1 <- read.table("edgeR_input_ENTREZ.txt", sep='\t',header=T, stringsAsFactors=F)
A<-df1$entrezgene_id
B <- df1[ ,229:454]
#B <- df1[ ,3:454]
C <- cbind(A,B)
C_nondup = C[-which(duplicated(C$A)),]
C_nondup2 = C_nondup[!(is.na(C_nondup$A) | C_nondup$A==""), ]
counts <- C_nondup2[,-1]       #Counts per sample
rownames(counts) <- as.character(C_nondup2[,1]) #entrezgene_id (num -> cha)

res <- CMScaller(emat=counts, RNAseq=TRUE, FDR=0.05)
write.table(df, "res_prediction.txt", col.names=NA, row.names=T, quote=F,sep='\t')

subPCA(emat = counts, class = res$prediction, normMethod = "quantile", pch=16, frame=FALSE)

cam <- CMSgsa(emat=counts, class=res$prediction, RNAseq=TRUE)

deg_all <- subDEG(emat=counts, class=res$prediction, doVoom=TRUE)
for (deg in deg_all) {
	deg$symbol <- fromTo(rownames(deg), id.in = "entrez", id.out="symbol", rough=TRUE)
#	subVolcano(deg, geneID="symbol", lfc=1, cex=0.66)
}
#head(deg2)
##logFC  AveExpr        t      P.Value    adj.P.Val        B symbol



