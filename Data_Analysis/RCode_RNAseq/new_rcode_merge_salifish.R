setwd("RNASeq_data")
dir()
rm(list=ls())


######data merging##################

GSM1616876_fv1 <- read.table(file="SRR1813698.txt", header= TRUE)
dim(GSM1616876_fv1)
GSM1616876_fv1 <- GSM1616876_fv1[order(GSM1616876_fv1$Name),]
head(GSM1616876_fv1)

GSM1616877_fv2 <- read.table(file="SRR1813699.txt", header= TRUE)
dim(GSM1616877_fv2)
GSM1616877_fv2 <- GSM1616877_fv2[order(GSM1616877_fv2$Name),]
head(GSM1616877_fv2)

GSM1616878_fv3 <- read.table(file="SRR1813700.txt", header= TRUE)
dim(GSM1616878_fv3)
GSM1616878_fv3 <- GSM1616878_fv3[order(GSM1616878_fv3$Name),]
head(GSM1616878_fv3)

GSM1616879_fz1 <- read.table(file="SRR1813701.txt", header= TRUE)
dim(GSM1616879_fz1)
GSM1616879_fz1 <- GSM1616879_fz1[order(GSM1616879_fz1$Name),]
head(GSM1616879_fz1)


GSM1616880_fz2 <- read.table(file="SRR1813702.txt", header= TRUE)
dim(GSM1616880_fz2)
GSM1616880_fz2 <- GSM1616880_fz2[order(GSM1616880_fz2$Name),]
head(GSM1616880_fz2)

GSM1616881_fz3 <- read.table(file="SRR1813703.txt", header= TRUE)
dim(GSM1616881_fz3)
GSM1616881_fz3 <- GSM1616881_fz3[order(GSM1616881_fz3$Name),]
head(GSM1616881_fz3)

GSM1616882_mv1 <- read.table(file="SRR1813704.txt", header= TRUE)
dim(GSM1616882_mv1)
GSM1616882_mv1 <- GSM1616882_mv1[order(GSM1616882_mv1$Name),]
head(GSM1616882_mv1)

GSM1616883_mv2 <- read.table(file="SRR1813705.txt", header= TRUE)
dim(GSM1616883_mv2)
GSM1616883_mv2 <- GSM1616883_mv2[order(GSM1616883_mv2$Name),]
head(GSM1616883_mv2)

GSM1616884_mv3 <- read.table(file="SRR1813706.txt", header= TRUE)
dim(GSM1616884_mv3)
GSM1616884_mv3 <- GSM1616884_mv3[order(GSM1616884_mv3$Name),]
head(GSM1616884_mv3)

GSM1616885_mz1 <- read.table(file="SRR1813707.txt", header= TRUE)
dim(GSM1616885_mz1)
GSM1616885_mz1 <- GSM1616885_mz1[order(GSM1616885_mz1$Name),]
head(GSM1616885_mz1)

GSM1616886_mz2 <- read.table(file="SRR1813708.txt", header= TRUE)
dim(GSM1616886_mz2)
GSM1616886_mz2 <- GSM1616886_mz2[order(GSM1616886_mz2$Name),]
head(GSM1616886_mz2)

GSM1616887_mz3 <- read.table(file="SRR1813709.txt", header= TRUE)
dim(GSM1616887_mz3)
GSM1616887_mz3 <- GSM1616887_mz3[order(GSM1616887_mz3$Name),]
head(GSM1616887_mz3)

#data file combination  
merged <- data.frame(cbind(GSM1616876_fv1, GSM1616877_fv2$TPM, GSM1616878_fv3$TPM, 
                           GSM1616879_fz1$TPM, GSM1616880_fz2$TPM, 
                           GSM1616881_fz3$TPM, GSM1616882_mv1$TPM, 
                           GSM1616883_mv2$TPM, GSM1616884_mv3$TPM, 
                           GSM1616885_mz1$TPM, GSM1616886_mz2$TPM, 
                           GSM1616887_mz3$TPM))
#merged$gene.no <- seq(1:30897)
colnames(merged) <- c("genes", "GSM1616876", "GSM1616877", 
                      "GSM1616878", "GSM1616879", "GSM1616880", 
                      "GSM1616881", "GSM1616882", "GSM1616883", 
                      "GSM1616884", "GSM1616885", "GSM1616886", 
                      "GSM1616887")

row.names(merged) <- merged$genes
merged <- data.frame(merged[,c(2:13)])
names(merged)
dim(merged)
head(merged)

####writing into aggregated file
write.table(merged, "RNAseq_new_merged_TPM.txt", sep="\t") 

###meta data: design matrix

meta_dat <- read.table(file="sailfish_file_table.txt")
colnames(meta_dat) <- c("SRR ID", "sample.no", "gender", "treatment")
row.names(meta_dat) <- meta_dat$sample.no
meta_dat <- data.frame(meta_dat[,c(1,3,4)])
head(meta_dat)

###another aggregation

allGene <- row.names(merged)

prepareData <- function(myGenes) {
  miniDat <- t(merged[myGenes, ])
  miniDat <- suppressWarnings(data.frame(gExp = as.vector(miniDat),
                                         gene = rep(colnames(miniDat), each = nrow(miniDat))))
  miniDat <- suppressWarnings(data.frame(meta_dat, miniDat))
  miniDat
}

#combine expression data and design
pDat <- prepareData(allGene)
head(pDat)
str(pDat)
dim(pDat)

####writing into aggregated file
write.table(pDat, "row.merged_TPM.txt", sep="\t")
