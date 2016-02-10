setwd("data")
dir()
rm(list=ls())
dat1 <- read.table(file="GSM1616876_FVEH_3_1.rpkm.txt")
dim(dat1)
colnames(dat1) <- c("genes", "GSM1616876_FVEH_3_1")

dat2 <- read.table(file="GSM1616877_FVEH_5_1.rpkm.txt")
dim(dat2)
colnames(dat2) <- c("genes", "GSM1616877_FVEH_5_1")

dat3 <- read.table(file="GSM1616878_FVEH_6_1.rpkm.txt")
dim(dat3)
colnames(dat3) <- c("genes", "GSM1616878_FVEH_6_1")

dat4 <- read.table(file="GSM1616879_FZEB_2_1.rpkm.txt")
dim(dat4)
colnames(dat4) <- c("genes", "GSM1616879_FZEB_2_1")

dat5 <- read.table(file="GSM1616880_FZEB_3_1.rpkm.txt")
dim(dat5)
colnames(dat5) <- c("genes", "GSM1616880_FZEB_3_1")

dat6 <- read.table(file="GSM1616881_FZEB_5_1.rpkm.txt")
dim(dat6)
colnames(dat6) <- c("genes", "GSM1616881_FZEB_5_1")

dat7 <- read.table(file="GSM1616882_MVEH_1_1.rpkm.txt")
dim(dat7)
colnames(dat7) <- c("genes", "GSM1616882_MVEH_1_1")

dat8 <- read.table(file="GSM1616883_MVEH_3_1.rpkm.txt")
dim(dat8)
colnames(dat8) <- c("genes", "GSM1616883_MVEH_3_1")

dat9 <- read.table(file="GSM1616884_MVEH_6_1.rpkm.txt")
dim(dat9)
colnames(dat9) <- c("genes", "GSM1616884_MVEH_6_1")

dat10 <- read.table(file="GSM1616885_MZEB_3_1.rpkm.txt")
dim(dat10)
colnames(dat10) <- c("genes", "GSM1616885_MZEB_3_1")

dat11 <- read.table(file="GSM1616886_MZEB_5_1.rpkm.txt")
dim(dat11)
colnames(dat11) <- c("genes", "GSM1616886_MZEB_5_1")

dat12 <- read.table(file="GSM1616887_MZEB_6_1.rpkm.txt")
dim(dat12)
colnames(dat12) <- c("genes", "GSM1616887_MZEB_6_1")


###data_merge
###ordering genes before merging

dat1 <- dat1[order(dat1$genes),]
dat2 <- dat2[order(dat2$genes),]
dat3 <- dat3[order(dat3$genes),]
dat4 <- dat4[order(dat4$genes),]
dat5 <- dat5[order(dat5$genes),]
dat6 <- dat6[order(dat6$genes),]
dat7 <- dat7[order(dat7$genes),]
dat8 <- dat8[order(dat8$genes),]
dat9 <- dat9[order(dat9$genes),]
dat10 <- dat10[order(dat10$genes),]
dat11 <- dat11[order(dat11$genes),]
dat12 <- dat12[order(dat12$genes),]

merged <- data.frame(cbind(dat1, dat2$GSM1616877_FVEH_5_1, dat3$GSM1616878_FVEH_6_1, 
                           dat4$GSM1616879_FZEB_2_1, dat5$GSM1616880_FZEB_3_1, 
                           dat6$GSM1616881_FZEB_5_1, dat7$GSM1616882_MVEH_1_1, 
                           dat8$GSM1616883_MVEH_3_1, dat9$GSM1616884_MVEH_6_1, 
                           dat10$GSM1616885_MZEB_3_1, dat11$GSM1616886_MZEB_5_1, 
                           dat12$GSM1616887_MZEB_6_1))
merged$gene.no <- seq(1:29516)
colnames(merged) <- c("genes", "GSM1616876_FVEH_3_1", "GSM1616877_FVEH_5_1", 
                      "GSM1616878_FVEH_6_1", "GSM1616879_FZEB_2_1", "GSM1616880_FZEB_3_1", 
                      "GSM1616881_FZEB_5_1", "GSM1616882_MVEH_1_1", "GSM1616883_MVEH_3_1", 
                      "GSM1616884_MVEH_6_1", "GSM1616885_MZEB_3_1", "GSM1616886_MZEB_5_1", 
                      "GSM1616887_MZEB_6_1", "gene.no")
names(merged)
dim(merged)
head(merged)
str(merged)


####writing into aggregated file
write.table(merged, "merged.txt", sep="\t") 

##other format
library(xlsx)
write.xlsx(merged, "merged.xlsx")

library(foreign)
write.foreign(merged, "mydata1.txt", "mydata.sps",   package="SPSS")
write.foreign(merged, "mydata2.txt", "mydata.sas",   package="SAS")


###check the data
asd <- read.table(file="merged.txt", header=T)
head(asd)
dim(asd)
str(asd)



#####another way of merging for analysis ####row merging 
rm(list=ls())
dat1 <- read.table(file="GSM1616876_FVEH_3_1.rpkm.txt")
dat2 <- read.table(file="GSM1616877_FVEH_5_1.rpkm.txt")
dat3 <- read.table(file="GSM1616878_FVEH_6_1.rpkm.txt")
dat4 <- read.table(file="GSM1616879_FZEB_2_1.rpkm.txt")
dat5 <- read.table(file="GSM1616880_FZEB_3_1.rpkm.txt")
dat6 <- read.table(file="GSM1616881_FZEB_5_1.rpkm.txt")
dat7 <- read.table(file="GSM1616882_MVEH_1_1.rpkm.txt")
dat8 <- read.table(file="GSM1616883_MVEH_3_1.rpkm.txt")
dat9 <- read.table(file="GSM1616884_MVEH_6_1.rpkm.txt")
dat10 <- read.table(file="GSM1616885_MZEB_3_1.rpkm.txt")
dat11 <- read.table(file="GSM1616886_MZEB_5_1.rpkm.txt")
dat11 <- read.table(file="GSM1616886_MZEB_5_1.rpkm.txt")
dat12 <- read.table(file="GSM1616887_MZEB_6_1.rpkm.txt")

head(dat1)
str(dat1)

dat1 <- dat1[order(dat1$V1),]
dat2 <- dat2[order(dat2$V1),]
dat3 <- dat3[order(dat3$V1),]
dat4 <- dat4[order(dat4$V1),]
dat5 <- dat5[order(dat5$V1),]
dat6 <- dat6[order(dat6$V1),]
dat7 <- dat7[order(dat7$V1),]
dat8 <- dat8[order(dat8$V1),]
dat9 <- dat9[order(dat9$V1),]
dat10 <- dat10[order(dat10$V1),]
dat11 <- dat11[order(dat11$V1),]
dat12 <- dat12[order(dat12$V1),]

dat1$gene.no <- seq(1:29516)
dat2$gene.no <- seq(1:29516)
dat3$gene.no <- seq(1:29516)
dat4$gene.no <- seq(1:29516)
dat5$gene.no <- seq(1:29516)
dat6$gene.no <- seq(1:29516)
dat7$gene.no <- seq(1:29516)
dat8$gene.no <- seq(1:29516)
dat9$gene.no <- seq(1:29516)
dat10$gene.no <- seq(1:29516)
dat11$gene.no <- seq(1:29516)
dat12$gene.no <- seq(1:29516)

dat1$replication <- rep(1, 29516)
dat2$replication <- rep(2, 29516)
dat3$replication <- rep(3, 29516)
dat4$replication <- rep(1, 29516)
dat5$replication <- rep(2, 29516)
dat6$replication <- rep(3, 29516)
dat7$replication <- rep(1, 29516)
dat8$replication <- rep(2, 29516)
dat9$replication <- rep(3, 29516)
dat10$replication <- rep(1, 29516)
dat11$replication <- rep(2, 29516)
dat12$replication <- rep(3, 29516)

dat1$group <- rep(1, 29516)
dat2$group <- rep(1, 29516)
dat3$group <- rep(1, 29516)
dat4$group <- rep(2, 29516)
dat5$group <- rep(2, 29516)
dat6$group <- rep(2, 29516)
dat7$group <- rep(3, 29516)
dat8$group <- rep(3, 29516)
dat9$group <- rep(3, 29516)
dat10$group <- rep(4, 29516)
dat11$group <- rep(4, 29516)
dat12$group <- rep(4, 29516)

datbind <- data.frame(rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, 
                            dat9, dat10, dat11, dat12))
datbind$group1 <- factor(datbind$group, levels = c(1, 2, 3, 4), labels = 
                           c("FVEH", "FZEB", "MVEH", "MZEB"))

colnames(datbind) <- c("genes", "rpkm_value", "gene.no", "replication", 
                       "group", "group.fct")
str(datbind)
head(datbind)

###example of extracting each gene information for all 4 groups

datbind[datbind$gene.no==1,]

####writing into aggregated file
write.table(datbind, "row.merged.txt", sep="\t")

#check
asd1 <- read.table(file="row.merged.txt", header=T)
str(asd1)
head(asd1)


##################program to compute p-values for all genes in 
############### comparison of all 4 groups

asd1 <- asd1[order(asd1$gene.no),]
resmat <- matrix(0, nrow = length(unique(asd1$gene.no)), ncol = 2)

resmat[,1] <- unique(asd1$gene.no)

for (i in 1: 29516) {
  resmat[i,2] <- anova(lm(asd1[asd1$gene.no==i,]$rpkm_value ~ asd1[asd1$gene.no==i,]$group.fct))$`Pr(>F)`[1]
}

resmat <- resmat
colnames(resmat) <- c("gene.no", "p.value")
####writing the anova value result
write.table(resmat, "anova.pval.txt", sep="\t")


