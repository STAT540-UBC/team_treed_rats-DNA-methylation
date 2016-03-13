setwd("RNASeq_data")
dir()

datn <- read.table(file="RNAseq_all_merged.txt", header = TRUE)
rownames(datn) <- datn$genes
datn <- datn[,c(2:14)]
head(datn)
str(datn)

countsfm <- datn[,c(1:3, 7:9, 13)]

##distribution of all female and male samples
par(mfrow = c(2,3))
hist(log(countsfm$GSM1616876_FVEH_3_1), xlab = "log RPKM count", main ="female replicate 1")
hist(log(countsfm$GSM1616877_FVEH_5_1), xlab = "log RPKM count", main ="female replicate 2")
hist(log(countsfm$GSM1616878_FVEH_6_1), xlab = "log RPKM count", main ="female replicate 3")
hist(log(countsfm$GSM1616882_MVEH_1_1), xlab = "log RPKM count", main ="male replicate 1")
hist(log(countsfm$GSM1616883_MVEH_3_1), xlab = "log RPKM count", main ="male replicate 2")
hist(log(countsfm$GSM1616884_MVEH_6_1), xlab = "log RPKM count", main ="male replicate 3")
dev.off()

###possible decision is to remove rpkm counts less than 1 (seems so from the histograms)
##revised rpkm counts after revision

#countsfmnz <- countsfm[countsfm[,1] >= 1 | countsfm[,2] >= 1 | countsfm[,3] >= 1 | countsfm[,4] >= 1 | countsfm[,5] >= 1 | countsfm[,6] >= 1 , ]
countsfm$fm.mean <- (countsfm$GSM1616876_FVEH_3_1 + countsfm$GSM1616877_FVEH_5_1 + countsfm$GSM1616878_FVEH_6_1)/3
countsfm$m.mean <- (countsfm$GSM1616882_MVEH_1_1 + countsfm$GSM1616883_MVEH_3_1 + countsfm$GSM1616884_MVEH_6_1)/3
head(countsfm)
str(countsfm)
dim(countsfm)

####distribution of female and male mean log rpkm counts (raw)
par(mfrow = c(1,2))
hist(log(countsfm$fm.mean), xlab = "mean log RPKM count", main ="female sample")
hist(log(countsfm$m.mean), xlab = "mean log RPKM count", main ="male sample")
dev.off()

countsfmnz <- countsfm[countsfm[,8] >= 1 & countsfm[,9] >= 1, ]
head(countsfmnz)
str(countsfmnz)
dim(countsfmnz)
#reduced to 14046 genes

####distribution of female and male mean log rpkm counts after reduction
par(mfrow = c(1,2))
hist(log(countsfmnz$fm.mean), xlab = "mean log RPKM count", main ="female sample")
hist(log(countsfmnz$m.mean), xlab = "mean log RPKM count", main ="male sample")
dev.off()

####fold change
##consider male in the denominator
countsfmnz$fold.change <- countsfmnz$fm.mean/countsfmnz$m.mean
#countsfmnz1 <- countsfmnz[(countsfmnz$fold.change >= 2 | countsfmnz$fold.change <= 0.5), ]
countsfmnz<- countsfmnz[order(countsfmnz$gene.no),]
countsfmnz$new.gene.order <- seq(1: dim(countsfmnz)[1])
head(countsfmnz)
str(countsfmnz)

dat_diffexp <- read.table(file="row.merged.txt", header = TRUE)
head(dat_diffexp)

library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)

###work with revised rpkm counts after removing genes below cutoff
rev_gene <- countsfmnz$gene.no
dat_diffexpn <- dat_diffexp[dat_diffexp$gene.no %in% rev_gene,]
dat_diffexpn<- dat_diffexpn[order(dat_diffexpn$gene.no),]
dat_diffexpn$new.gene.order <- rep(1:length(rev_gene), each = 12)


library(gdata)
dat_diff.fvmv <- dat_diffexpn[dat_diffexpn$group.fct == "FVEH" | dat_diffexpn$group.fct == "MVEH",]
dim(dat_diff.fvmv)

###individual genes
dat_diff.fvmv<- dat_diff.fvmv[order(dat_diff.fvmv$new.gene.order),]
resfvmv <- matrix(0, nrow = length(unique(dat_diff.fvmv$gene.no)), ncol = 3)
resfvmv[,1] <- unique(dat_diff.fvmv$gene.no)
resfvmv[,2] <- unique(dat_diff.fvmv$new.gene.order)

for (i in 1: dim(resfvmv)[1]) {
  resfvmv[i,3] <- t.test(dat_diff.fvmv[dat_diff.fvmv$new.gene.order==i,]$rpkm_value ~ dat_diff.fvmv[dat_diff.fvmv$new.gene.order==i,]$group)$p.value
}

resfvmv <- resfvmv
resfvmv <- data.frame(resfvmv)
#resfvmv$genes <- factor(levels(droplevels(dat_fvmv$genes)))
colnames(resfvmv) <- c("gene.no", "new.gene.order", "fvmv.pval")

#dim(resfvmv[resfvmv$fvmv.pval < 0.05,])  ###251 de genes

###differentially expression result sheet preparation

diff.exp.fm <- data.frame(countsfmnz[,7:10], resfvmv[,3])
str(diff.exp.fm)
head(diff.exp.fm) 
colnames(diff.exp.fm) <- c("gene.no", "fm.mean", "m.mean", "fold.change", "raw.pvalue") 

####adjustment in p-value for multiple comparison test: three methods

diff.exp.fm$pvalue.fdr.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "BH"), 4)
diff.exp.fm$pvalue.hb.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "BY"), 4)
#diff.exp.fm$pvalue.bn.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "bonferroni"), 4)

##discard the NA values
#fvmvna <- na.omit(diff.exp.fm)
#dim(fvmvna)

###checking the genes signifcantly different at 5% level of significance between male anf female

sig.fvmv1 <- diff.exp.fm[diff.exp.fm$pvalue.fdr.adj < 0.05,]
dim(sig.fvmv1)
View(sig.fvmv1)

sig.fvmv2 <- diff.exp.fm[diff.exp.fm$pvalue.hb.adj < 0.05,]
dim(sig.fvmv2)
#View(sig.fvmv2)

