setwd("data")
rm(list=ls())

###read the merged data
asd <- read.table(file="RNAseq_all_merged.txt", header=T)
head(asd)
dim(asd)
str(asd)

#read the row.merged data
asd1 <- read.table(file="row.merged.txt", header=T)
str(asd1)
head(asd1)

####read the p-value data
bdat <- read.table(file = "RHgroups_anova.pval.txt", header = T)
head(bdat)
View(bdat)

##discard the NA values
bdat <- na.omit(bdat)
dim(bdat)

###checking the genes signifcant at 5% level of significance for further analysis

bdat1 <- bdat[bdat$p.value < 0.05,]
dim(bdat1)
View(bdat1)
###found 1325 genes as differential ones between any of the 4 groups
### next step is to reduce the workload and work on these specific genes only

####reduce the size of the merged datasets according to only the significant ones

library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)

###significant genes at 5% level of significance
sig.gene.no <- bdat1$gene.no
tbl_anovdat <- tbl_df(bdat1)

###reduction of the both data file 

sig_anovdat1 <- asd[asd$gene.no %in% sig.gene.no,]
dim(sig_anovdat1)  ###check the reduced data for number of genes significant
head(sig_anovdat1)

sig_anovdat2 <- asd1[asd1$gene.no %in% sig.gene.no,]
length(unique(sig_anovdat2$gene.no)) ###check the reduced data for number of genes significant
head(sig_anovdat2)

####now t-tests among the groups among the reduced result

###seems not important
anova(lm(sig_anovdat2$rpkm_value ~ sig_anovdat2$group.fct))
pairwise.t.test(sig_anovdat2$rpkm_value, sig_anovdat2$group.fct, 
                pool.sd=FALSE, p.adjust.method="none")
t.test(rpkm_value ~ group.fct, data=sig_anovdat2[sig_anovdat2$group.fct == "FVEH" | sig_anovdat2$group.fct == "MVEH",])


####check across the genes 
sig_anovdat2 <- sig_anovdat2[order(sig_anovdat2$gene.no),]
tbl_sig_anovdat2 <- tbl_df(sig_anovdat2)
sig_anovdat2$new.gene.order <- rep(1:1325, each = 12)
head(sig_anovdat2)

###FVEH vs. MVEH (Male vs Female)
library(gdata)
dat_fvmv <- sig_anovdat2[sig_anovdat2$group.fct == "FVEH" | sig_anovdat2$group.fct == "MVEH",]
dim(dat_fvmv)

t.test(dat_fvmv[dat_fvmv$new.gene.order==1,]$rpkm_value ~ dat_fvmv[dat_fvmv$new.gene.order==1,]$group)

t.test(dat_fvmv[dat_fvmv$new.gene.order==1,]$rpkm_value ~ dat_fvmv[dat_fvmv$new.gene.order==1,]$group)$p.value

###individual genes
dat_fvmv <- dat_fvmv[order(dat_fvmv$gene.no),]
resfvmv <- matrix(0, nrow = length(unique(dat_fvmv$gene.no)), ncol = 3)
resfvmv[,1] <- unique(dat_fvmv$gene.no)
resfvmv[,2] <- unique(dat_fvmv$new.gene.order)

for (i in 1: dim(resfvmv)[1]) {
  resfvmv[i,3] <- t.test(dat_fvmv[dat_fvmv$new.gene.order==i,]$rpkm_value ~ dat_fvmv[dat_fvmv$new.gene.order==i,]$group)$p.value
}

resfvmv <- resfvmv
resfvmv <- data.frame(resfvmv)
resfvmv$genes <- factor(levels(droplevels(dat_fvmv$genes)))
colnames(resfvmv) <- c("gene.no", "new.gene.order", "fvmv.pval")

##discard the NA values
fvmvna <- na.omit(resfvmv)
dim(fvmvna)

###checking the genes signifcantly different at 5% level of significance between male anf female

sig.fvmv <- fvmvna[fvmvna$fvmv.pval < 0.05,]
dim(sig.fvmv)
View(sig.fvmv)

###found 154 genes as differential ones between between male anf female

genes.fvmv <- sig.fvmv$gene.no

###reduction of the both data file to analyze 

sigt.fvmv1 <- asd[asd$gene.no %in% genes.fvmv,]
dim(sigt.fvmv1)  ###check the reduced data for number of genes significant
head(sigt.fvmv1)

sigt.fvmv2 <- asd1[asd1$gene.no %in% genes.fvmv,]
length(unique(sigt.fvmv2$gene.no)) ###check the reduced data for number of genes significant
head(sigt.fvmv2)


###FVEH vs. FZEB (Female vs Female Treatment)

dat_fvfz <- sig_anovdat2[sig_anovdat2$group.fct == "FVEH" | 
                           sig_anovdat2$group.fct == "FZEB",]
dim(dat_fvfz)

t.test(dat_fvfz[dat_fvfz$new.gene.order==1,]$rpkm_value ~ dat_fvfz[dat_fvfz$new.gene.order==1,]$group)

###individual genes
dat_fvfz <- dat_fvfz[order(dat_fvfz$gene.no),]
resfvfz <- matrix(0, nrow = length(unique(dat_fvfz$gene.no)), ncol = 3)
resfvfz[,1] <- unique(dat_fvfz$gene.no)
resfvfz[,2] <- unique(dat_fvfz$new.gene.order)

for (i in 1: dim(resfvfz)[1]) {
  resfvfz[i,3] <- t.test(dat_fvfz[dat_fvfz$new.gene.order==i,]$rpkm_value ~ dat_fvfz[dat_fvfz$new.gene.order==i,]$group)$p.value
}

resfvfz <- resfvfz
resfvfz <- data.frame(resfvfz)
resfvfz$genes <- factor(levels(droplevels(dat_fvfz$genes)))
colnames(resfvfz) <- c("gene.no", "new.gene.order", "fvfz.pval", "genes")

##discard the NA values
fvfzna <- na.omit(resfvfz)
dim(fvfzna)

###checking the genes signifcantly different at 5% level of significance between female and female treatment

sig.fvfz <- fvfzna[fvfzna$fvfz.pval < 0.05,]
dim(sig.fvfz)
View(sig.fvfz)

###found 232 genes as differential ones between between female and female treatment

genes.fvfz <- sig.fvfz$gene.no

###reduction of the both data file to analyze 

sigt.fvfz1 <- asd[asd$gene.no %in% genes.fvfz,]
dim(sigt.fvfz1)  ###check the reduced data for number of genes significant
head(sigt.fvfz1)

sigt.fvfz2 <- asd1[asd1$gene.no %in% genes.fvfz,]
length(unique(sigt.fvfz2$gene.no)) ###check the reduced data for number of genes significant
head(sigt.fvfz2)


###MVEH vs. MZEB

dat_mvmz <- sig_anovdat2[sig_anovdat2$group.fct == "MZEB" | sig_anovdat2$group.fct == "MVEH",]
dim(dat_mvmz)

###MVEH vs. FZEB

dat_mvfz <- sig_anovdat2[sig_anovdat2$group.fct == "FZEB" | sig_anovdat2$group.fct == "MVEH",]
dim(dat_mvfz)

###MZEB vs. FZEB

dat_fzmz <- sig_anovdat2[sig_anovdat2$group.fct == "FZEB" | sig_anovdat2$group.fct == "MZEB",]
dim(dat_fzmz)

