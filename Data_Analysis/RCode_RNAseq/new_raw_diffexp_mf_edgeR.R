############### R code by Rashed ################

setwd("../RNASeq_data/new_data_Tony_TPM")
dir()
rm(list = ls())

###Raw counts
datn <- read.table(file="RNAseq_new_merged_raw.txt", header = TRUE)
#input counts need to be integer so that rounding all entries into the next integer
datn <- round(datn, digits = 0)
head(datn)
str(datn)

###comparing male and female vehicle data preparation

###count data
countsfm <- datn[,c(1:3, 7:9)]
head(countsfm)
str(countsfm)
dim(countsfm)

#row_sub = apply(countsfm, 1, function(row) all(row !=0 ))
#countsfmn <- countsfm[row_sub,]
#dim(countsfmn)

##densityplot
require(dplyr)
require(ggplot2)
require(tidyr)

expression_mf_long <- countsfm %>% 
  add_rownames() %>%
  gather(key = sample, value = raw_counts, ... = -rowname)

expression_mf_long %>%
  ggplot(aes(raw_counts, color = sample)) +
  geom_density() +
  scale_x_continuous(trans = "log10") +
  geom_vline(xintercept = 1.5)


##Calculate edgeR object DGEList data class
library(edgeR)

#factor variable group (1 means female 2 means male)
group <- factor(rep(1:2, each = 3))
y <- DGEList(counts=countsfm, group=group)
###check y
y

###filtering the low counts

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

####differential expression analysis consideration

# The most obvious technical factor that affects the read counts, other than gene expression
#levels, is the sequencing depth of each RNA sample. 
# RNA composition effect is adjusted for, the remaining genes may falsely
#appear to be down-regulated in that sample
# For example, users should not enter RPKM or FPKM values to edgeR in place of read counts.

# The classic edgeR functions estimateCommonDisp and exactTest produce a matrix of pseudocounts 
#as part of the output object. The pseudo-counts are used internally to speed up
#computation of the conditional likelihood used for dispersion estimation and exact tests
#in the classic edgeR pipeline.

# The classic edgeR functions estimateCommonDisp and exactTest produce a matrix of pseudocounts 
#as part of the output object. The pseudo-counts are used internally to speed up
#computation of the conditional likelihood used for dispersion estimation and exact tests in the 
#classic edgeR pipeline. The qCML method is only applicable on datasets with a single factor design

# Testing for DE genes: The exact test is based on the qCML methods. Knowing the conditional 
#distribution for the sum of counts in a group, we can compute exact p-values by summing over all 
#sums of counts that have a probability less than the probability under the null hypothesis of the
#observed sum of counts. The exact test for the negative binomial distribution has strong
#parallels with Fisher’s exact test. The testing can be done by using the function exactTest(), 
#and the function allows both common dispersion and tagwise dispersion approaches. 


#######################  Negative Binomial Model  ############################
#design
design <- model.matrix(~group)

#Normalization
y <- calcNormFactors(y, method= "TMM")

#Estimating dispersions (without design goes to the classic version)
#To estimate common dispersion and tagwise dispersions in one run (recommended):
 y <- estimateDisp(y, design = design)
#Alternatively, to estimate common dispersion:
  y <- estimateCommonDisp(y)
#Then to estimate tagwise dispersions:
  y <- estimateTagwiseDisp(y)

#Testing for DE genes

et <- exactTest(y)
genes.de <- topTags(et)  #10 most differential genes

###genes sig with fdr 0.05

genes.de.fdr.05 <- topTags(et, n = Inf, p.value=0.05)
genes.de.fdr.05
(genecounts <- nrow(genes.de.fdr.05))

#########DE genes for different normalizations (Negative Binomial distribution):

#1. Without Design Matrix: TMM, RLE, upperquartile, none: 33 genes
#2. With Design Matrix: TMM, RLE, upperquartile, none: 21 genes


################ edgeR GLM DE #############################

rm(list = ls())

###Raw counts
datn <- read.table(file="../RNASeq_data/new_data_Tony_TPM/RNAseq_new_merged_raw.txt", header = TRUE)[,c(1:3, 7:9)]
#input counts need to be integer so that rounding all entries into the next integer
datn <- round(datn, digits = 0)

###comparing male and female vehicle data preparation
countsfm <- datn[,c(1:3, 7:9)]

##Calculate edgeR object DGEList data class for GLM edgeR
library(edgeR)

#factor variable group (1 means female 2 means male)
group <- factor(rep(1:2, each = 3))

dge.glm <- DGEList(counts=countsfm, group=group)
###check dge.glm
dge.glm
str(dge.glm)
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[1]])

#design
design <- model.matrix(~group)
design

#dispersion paramete estimate
# to estimate common dispersion:
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm,design, verbose=TRUE)
# to estimate trendwise dispersion:
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, design)
# to estimate tagwise dispersion:
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
#check plot
plotBCV(dge.glm.tag.disp)

###GLM model fitting
fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))

#DE analysis
lrt <- glmLRT(fit,coef=2) #on specific coefficent
topTags(lrt)  #top genes

#finding DE genes with a specific FDA: 0.05 here

tt.glm <- topTags(lrt, n=Inf)
class(tt.glm)
tt.glm.fdr0.05 <- tt.glm$table[tt.glm$table$FDR < 0.05,]
tt.glm.fdr0.05

#number of DE genes at FDR 0.05
nrow(tt.glm.fdr0.05)

#interesting genes (need to specify the criterion)
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50,])
cpm(dge.glm.tag.disp)[interestingSamples,]

#summary results
summary(de.glm <- decideTestsDGE(lrt, p=0.05, adjust="BH"))

#plotting the tagwise log fold changes against log-cpm
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags=tags.glm)
abline(h=c(-2,2),col="blue")

#makes more sense in glm fitting: found 62 DE genes with down regulated 36 genes 
#and up regulated 26 genes


################## Work with DESeq Package in R ####################

rm(list = ls())

###Raw counts
datn <- read.table(file="RNAseq_new_merged_raw.txt", header = TRUE)
#input counts need to be integer so that rounding all entries into the next integer
datn <- round(datn, digits = 0)

###comparing male and female vehicle data preparation
countsfm <- datn[,c(1:3, 7:9)]

#factor variable group (1 means female 2 means male)
group <- factor(rep(1:2, each = 3))
library(DESeq)

#reading in the same count table data and grouping information for input of DESeq DEA
deSeqDat <- newCountDataSet(countsfm, group)
head(counts(deSeqDat))

#Next, we estimate the size factors to account for differences in library coverage 
#and estimate the variance:

deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)

#estimate dispersion parameter
deSeqDat <- estimateDispersions(deSeqDat)
#plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)

#fit the model and examine the results

results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)
dim(results)

#finding DE genes with a specific FDA: 0.05 here
results1 <- na.omit(results)
deseq.nbinom.fdr0.05 <- results1[results1$padj <0.05,]
deseq.nbinom.fdr0.05

#number of DE genes at FDR 0.05
nrow(deseq.nbinom.fdr0.05)

#plot the results
plotMA(results)

#DESeq finds 21 DE genes like as edgeR with design matrix provided


####################### Voom & limma  ########################

rm(list = ls())

###Raw counts
datn <- read.table(file="RNAseq_new_merged_raw.txt", header = TRUE)
#input counts need to be integer so that rounding all entries into the next integer
datn <- round(datn, digits = 0)

###comparing male and female vehicle data preparation
countsfm <- datn[,c(1:3, 7:9)]

#factor variable group (1 means female 2 means male)
group <- factor(rep(1:2, each = 3))

#design
design <- model.matrix(~group)

library(limma)

#normalizing factor for raw count data matrix
norm.factor <- calcNormFactors(countsfm)

#voom transformation
dat.voomed <- voom(countsfm, design, plot=TRUE,lib.size=colSums(countsfm)*norm.factor)
dat.voomed

##fit the model, for all genes at once, and use eBayes() to moderate 
#the estimated error variances:

limmafit <- lmFit(dat.voomed, design)
limmafiteB <- eBayes(limmafit)

#top 10 performing genes
topTable(limmafiteB, coef="group2", sort.by="p")

##find out genes with a specific fdr 0.05

limmafitHits <- topTable(limmafiteB, coef=2, n = Inf, p.value= 0.05)
limmafitHits 
(voomcounts.de <- nrow(limmafitHits))

#This gives only 12 DE genes
