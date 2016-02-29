# Finding differentially expressed genes for RNAseq `rpkm` read counts data using multiple `t` test and fold change
Rashed  
February 29, 2016  

# Steps to find out differentially expressed genes using multiple `t` test and fold change

The analysis pipeline consists of the following steps.

 1. Check the distribution of the log rpkm counts of individual samples (replicates) for male and female
 2. Compute the mean of rpkm counts over male and female sample and check the distribution of the mean log rpkm counts for male and female
 3. Filtering the rpkm counts based on discarding some rpkm counts after checking the distribution above
 4. Calculate the fold change and prepare the data for t-test
 5. Perform t-test and adjust p values using `FDR` (or other) for multiple testing

# Check the distribution of the log rpkm counts of individual samples (replicates) for male and female

## Reading Data

```r
#setwd("team_treed_rats-DNA-methylation/RNASeq_data/RNASeq_data")
#dir()
library(knitr)
library(rmarkdown)
library(xtable)

datn <- read.table(file=
                     "C:/Users/Rashed/Documents/all_github/team_treed_rats-DNA-methylation/RNASeq_data/RNAseq_all_merged.txt", 
                   header = TRUE)
rownames(datn) <- datn$genes
datn <- datn[,c(2:14)]
head(datn)
```

```
##                    GSM1616876_FVEH_3_1 GSM1616877_FVEH_5_1
## ENSRNOG00000000001           0.3067553          0.34937193
## ENSRNOG00000000007          76.5608522         72.84595259
## ENSRNOG00000000008           0.0557963          0.01308340
## ENSRNOG00000000009           0.2127083          0.07980306
## ENSRNOG00000000010           7.3180708         14.39158730
## ENSRNOG00000000012           0.5348094          0.68972645
##                    GSM1616878_FVEH_6_1 GSM1616879_FZEB_2_1
## ENSRNOG00000000001           0.2068158           0.3142311
## ENSRNOG00000000007          69.4314626          68.4164006
## ENSRNOG00000000008           0.5595697           0.0000000
## ENSRNOG00000000009           0.2509657           0.2542075
## ENSRNOG00000000010          26.5290822          20.2018857
## ENSRNOG00000000012           0.3943746           0.5592566
##                    GSM1616880_FZEB_3_1 GSM1616881_FZEB_5_1
## ENSRNOG00000000001          0.11537004           0.2086471
## ENSRNOG00000000007         74.40058336          69.0384660
## ENSRNOG00000000008          0.08813648           0.8799938
## ENSRNOG00000000009          0.00000000           0.0000000
## ENSRNOG00000000010         13.63452893           8.5958703
## ENSRNOG00000000012          0.63359340           0.7161597
##                    GSM1616882_MVEH_1_1 GSM1616883_MVEH_3_1
## ENSRNOG00000000001           0.2421934          0.22098735
## ENSRNOG00000000007          83.3136515         67.79241856
## ENSRNOG00000000008           0.0000000          0.02557915
## ENSRNOG00000000009           0.1410698          0.07801088
## ENSRNOG00000000010           8.0890022         18.30070229
## ENSRNOG00000000012           0.9606181          0.98070826
##                    GSM1616884_MVEH_6_1 GSM1616885_MZEB_3_1
## ENSRNOG00000000001          0.24217070          0.29265232
## ENSRNOG00000000007         74.78969651         69.51902744
## ENSRNOG00000000008          0.02371864          0.06986581
## ENSRNOG00000000009          0.07233670          0.07102520
## ENSRNOG00000000010         14.41252577         15.42670804
## ENSRNOG00000000012          0.79570374          0.33483308
##                    GSM1616886_MZEB_5_1 GSM1616887_MZEB_6_1 gene.no
## ENSRNOG00000000001           0.1406568          0.13247706       1
## ENSRNOG00000000007          80.7306331         76.36488333       2
## ENSRNOG00000000008           0.0000000          0.01686755       3
## ENSRNOG00000000009           0.1365467          0.15432717       4
## ENSRNOG00000000010           8.7244650         17.07601107       5
## ENSRNOG00000000012           0.6437203          1.05089455       6
```

```r
#str(datn)

###comparing male and female vehicle

####data preparation

###count data
#can input raw or norrmalized count. here we put rpkm values as count
countsfm <- datn[,c(1:3, 7:9, 13)]
#head(countsfm)
#str(countsfm)
#dim(countsfm)
```


## Check the distribution of individual samples



```r
##distribution of all female and male samples
par(mfrow = c(2,3))
hist(log(countsfm$GSM1616876_FVEH_3_1), xlab = "log RPKM count", main ="female replicate 1")
hist(log(countsfm$GSM1616877_FVEH_5_1), xlab = "log RPKM count", main ="female replicate 2")
hist(log(countsfm$GSM1616878_FVEH_6_1), xlab = "log RPKM count", main ="female replicate 3")
hist(log(countsfm$GSM1616882_MVEH_1_1), xlab = "log RPKM count", main ="male replicate 1")
hist(log(countsfm$GSM1616883_MVEH_3_1), xlab = "log RPKM count", main ="male replicate 2")
hist(log(countsfm$GSM1616884_MVEH_6_1), xlab = "log RPKM count", main ="male replicate 3")
```

![](1-RNA_Seq_t_test_differential_mf_files/figure-html/ch12-1.png)<!-- -->

```r
#dev.off()

###possible decision is to remove rpkm counts less than 1 (seems so from the histograms)
##revised rpkm counts after revision
```

#Compute the mean of rpkm counts over male and female sample and check the distribution of the mean log rpkm counts for male and female


##ompute the mean of rpkm counts over male and female sample



```r
countsfm$fm.mean <- (countsfm$GSM1616876_FVEH_3_1 + countsfm$GSM1616877_FVEH_5_1 + countsfm$GSM1616878_FVEH_6_1)/3
countsfm$m.mean <- (countsfm$GSM1616882_MVEH_1_1 + countsfm$GSM1616883_MVEH_3_1 + countsfm$GSM1616884_MVEH_6_1)/3
head(countsfm)
```

```
##                    GSM1616876_FVEH_3_1 GSM1616877_FVEH_5_1
## ENSRNOG00000000001           0.3067553          0.34937193
## ENSRNOG00000000007          76.5608522         72.84595259
## ENSRNOG00000000008           0.0557963          0.01308340
## ENSRNOG00000000009           0.2127083          0.07980306
## ENSRNOG00000000010           7.3180708         14.39158730
## ENSRNOG00000000012           0.5348094          0.68972645
##                    GSM1616878_FVEH_6_1 GSM1616882_MVEH_1_1
## ENSRNOG00000000001           0.2068158           0.2421934
## ENSRNOG00000000007          69.4314626          83.3136515
## ENSRNOG00000000008           0.5595697           0.0000000
## ENSRNOG00000000009           0.2509657           0.1410698
## ENSRNOG00000000010          26.5290822           8.0890022
## ENSRNOG00000000012           0.3943746           0.9606181
##                    GSM1616883_MVEH_3_1 GSM1616884_MVEH_6_1 gene.no
## ENSRNOG00000000001          0.22098735          0.24217070       1
## ENSRNOG00000000007         67.79241856         74.78969651       2
## ENSRNOG00000000008          0.02557915          0.02371864       3
## ENSRNOG00000000009          0.07801088          0.07233670       4
## ENSRNOG00000000010         18.30070229         14.41252577       5
## ENSRNOG00000000012          0.98070826          0.79570374       6
##                       fm.mean      m.mean
## ENSRNOG00000000001  0.2876477  0.23511715
## ENSRNOG00000000007 72.9460891 75.29858885
## ENSRNOG00000000008  0.2094831  0.01643260
## ENSRNOG00000000009  0.1811590  0.09713913
## ENSRNOG00000000010 16.0795801 13.60074341
## ENSRNOG00000000012  0.5396368  0.91234338
```

```r
#str(countsfm)
dim(countsfm)
```

```
## [1] 29516     9
```

##check the distribution of the mean log rpkm counts for male and female


```r
####distribution of female and male mean log rpkm counts (raw)
par(mfrow = c(1,2))
hist(log(countsfm$fm.mean), xlab = "mean log RPKM count", main ="female sample")
hist(log(countsfm$m.mean), xlab = "mean log RPKM count", main ="male sample")
```

![](1-RNA_Seq_t_test_differential_mf_files/figure-html/ch14-1.png)<!-- -->

```r
#dev.off()
```

#Filtering the rpkm counts based on discarding some rpkm counts after checking the distribution above

##Filtering the rpkm counts discarding genes with mean rpkm counts less than 1 



```r
countsfmnz <- countsfm[countsfm[,8] >= 1 & countsfm[,9] >= 1, ]
head(countsfmnz)
```

```
##                    GSM1616876_FVEH_3_1 GSM1616877_FVEH_5_1
## ENSRNOG00000000007           76.560852           72.845953
## ENSRNOG00000000010            7.318071           14.391587
## ENSRNOG00000000014            1.292220            1.837586
## ENSRNOG00000000021            8.038146            6.282750
## ENSRNOG00000000024           11.175178           13.266981
## ENSRNOG00000000028            1.382465            1.640918
##                    GSM1616878_FVEH_6_1 GSM1616882_MVEH_1_1
## ENSRNOG00000000007          69.4314626          83.3136515
## ENSRNOG00000000010          26.5290822           8.0890022
## ENSRNOG00000000014           0.3442727           0.4146824
## ENSRNOG00000000021           7.4405996           7.9819941
## ENSRNOG00000000024          11.7098678          11.2726827
## ENSRNOG00000000028           1.8043282           1.4604873
##                    GSM1616883_MVEH_3_1 GSM1616884_MVEH_6_1 gene.no
## ENSRNOG00000000007           67.792419           74.789697       2
## ENSRNOG00000000010           18.300702           14.412526       5
## ENSRNOG00000000014            1.796318            1.630222       7
## ENSRNOG00000000021            6.770800            7.861791       9
## ENSRNOG00000000024           12.037560           12.457588      10
## ENSRNOG00000000028            1.559198            1.528999      11
##                      fm.mean    m.mean
## ENSRNOG00000000007 72.946089 75.298589
## ENSRNOG00000000010 16.079580 13.600743
## ENSRNOG00000000014  1.158026  1.280408
## ENSRNOG00000000021  7.253832  7.538195
## ENSRNOG00000000024 12.050676 11.922610
## ENSRNOG00000000028  1.609237  1.516228
```

```r
#str(countsfmnz)
dim(countsfmnz)
```

```
## [1] 14046     9
```

```r
#reduced to 14046 genes
```

## Distribution of the revised results 


```r
####distribution of female and male mean log rpkm counts after reduction
par(mfrow = c(1,2))
hist(log(countsfmnz$fm.mean), xlab = "mean log RPKM count", main ="female sample")
hist(log(countsfmnz$m.mean), xlab = "mean log RPKM count", main ="male sample")
```

![](1-RNA_Seq_t_test_differential_mf_files/figure-html/ch16-1.png)<!-- -->

```r
#dev.off()
```

#Calculate the fold change and prepare the data for t-test


```r
####fold change
##consider male in the denominator
countsfmnz$fold.change <- countsfmnz$fm.mean/countsfmnz$m.mean
#countsfmnz1 <- countsfmnz[(countsfmnz$fold.change >= 2 | countsfmnz$fold.change <= 0.5), ]
countsfmnz<- countsfmnz[order(countsfmnz$gene.no),]
countsfmnz$new.gene.order <- seq(1: dim(countsfmnz)[1])
head(countsfmnz)
```

```
##                    GSM1616876_FVEH_3_1 GSM1616877_FVEH_5_1
## ENSRNOG00000000007           76.560852           72.845953
## ENSRNOG00000000010            7.318071           14.391587
## ENSRNOG00000000014            1.292220            1.837586
## ENSRNOG00000000021            8.038146            6.282750
## ENSRNOG00000000024           11.175178           13.266981
## ENSRNOG00000000028            1.382465            1.640918
##                    GSM1616878_FVEH_6_1 GSM1616882_MVEH_1_1
## ENSRNOG00000000007          69.4314626          83.3136515
## ENSRNOG00000000010          26.5290822           8.0890022
## ENSRNOG00000000014           0.3442727           0.4146824
## ENSRNOG00000000021           7.4405996           7.9819941
## ENSRNOG00000000024          11.7098678          11.2726827
## ENSRNOG00000000028           1.8043282           1.4604873
##                    GSM1616883_MVEH_3_1 GSM1616884_MVEH_6_1 gene.no
## ENSRNOG00000000007           67.792419           74.789697       2
## ENSRNOG00000000010           18.300702           14.412526       5
## ENSRNOG00000000014            1.796318            1.630222       7
## ENSRNOG00000000021            6.770800            7.861791       9
## ENSRNOG00000000024           12.037560           12.457588      10
## ENSRNOG00000000028            1.559198            1.528999      11
##                      fm.mean    m.mean fold.change new.gene.order
## ENSRNOG00000000007 72.946089 75.298589   0.9687577              1
## ENSRNOG00000000010 16.079580 13.600743   1.1822574              2
## ENSRNOG00000000014  1.158026  1.280408   0.9044199              3
## ENSRNOG00000000021  7.253832  7.538195   0.9622770              4
## ENSRNOG00000000024 12.050676 11.922610   1.0107414              5
## ENSRNOG00000000028  1.609237  1.516228   1.0613422              6
```

```r
str(countsfmnz)
```

```
## 'data.frame':	14046 obs. of  11 variables:
##  $ GSM1616876_FVEH_3_1: num  76.56 7.32 1.29 8.04 11.18 ...
##  $ GSM1616877_FVEH_5_1: num  72.85 14.39 1.84 6.28 13.27 ...
##  $ GSM1616878_FVEH_6_1: num  69.431 26.529 0.344 7.441 11.71 ...
##  $ GSM1616882_MVEH_1_1: num  83.314 8.089 0.415 7.982 11.273 ...
##  $ GSM1616883_MVEH_3_1: num  67.79 18.3 1.8 6.77 12.04 ...
##  $ GSM1616884_MVEH_6_1: num  74.79 14.41 1.63 7.86 12.46 ...
##  $ gene.no            : int  2 5 7 9 10 11 13 16 19 20 ...
##  $ fm.mean            : num  72.95 16.08 1.16 7.25 12.05 ...
##  $ m.mean             : num  75.3 13.6 1.28 7.54 11.92 ...
##  $ fold.change        : num  0.969 1.182 0.904 0.962 1.011 ...
##  $ new.gene.order     : int  1 2 3 4 5 6 7 8 9 10 ...
```

```r
dat_diffexp <- read.table(file=
            "C:/Users/Rashed/Documents/all_github/team_treed_rats-DNA-methylation/RNASeq_data/row.merged.txt", header = TRUE)
head(dat_diffexp)
```

```
##                genes rpkm_value gene.no replication group group.fct
## 1 ENSRNOG00000000001  0.3067553       1           1     1      FVEH
## 2 ENSRNOG00000000007 76.5608522       2           1     1      FVEH
## 3 ENSRNOG00000000008  0.0557963       3           1     1      FVEH
## 4 ENSRNOG00000000009  0.2127083       4           1     1      FVEH
## 5 ENSRNOG00000000010  7.3180708       5           1     1      FVEH
## 6 ENSRNOG00000000012  0.5348094       6           1     1      FVEH
```

```r
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))

###work with revised rpkm counts after removing genes below cutoff
rev_gene <- countsfmnz$gene.no
dat_diffexpn <- dat_diffexp[dat_diffexp$gene.no %in% rev_gene,]
dat_diffexpn<- dat_diffexpn[order(dat_diffexpn$gene.no),]
dat_diffexpn$new.gene.order <- rep(1:length(rev_gene), each = 12)

##Data preparation for t-test
suppressPackageStartupMessages(library(gdata))
##data with male and female samples only
dat_diff.fvmv <- dat_diffexpn[dat_diffexpn$group.fct == "FVEH" | dat_diffexpn$group.fct == "MVEH",]
dim(dat_diff.fvmv)
```

```
## [1] 84276     7
```


# Perform t-test and adjust p values using `FDR` (or other) for multiple testing

## Perform t-test and observe the results


```r
###individual genes
dat_diff.fvmv<- dat_diff.fvmv[order(dat_diff.fvmv$new.gene.order),]

##storing the raw p-values
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
```

Now we observe the results for raw p values (before adjustment for multiple testing) in the following table. Here we have found 251 differentially expressed genes (seems too low!!!) at 5% level of significance. 


```r
###differentially expression result sheet preparation

diff.exp.fmraw <- data.frame(countsfmnz[,7:10], resfvmv[,3])
diff.exp.fmraw <- diff.exp.fmraw[diff.exp.fmraw$resfvmv...3. < 0.05, ]
str(diff.exp.fmraw)
```

```
## 'data.frame':	251 obs. of  5 variables:
##  $ gene.no     : int  224 383 388 491 555 557 580 681 784 824 ...
##  $ fm.mean     : num  13.63 10.55 2.76 17.01 18.99 ...
##  $ m.mean      : num  11.83 12.29 3.21 17.69 20.07 ...
##  $ fold.change : num  1.152 0.859 0.859 0.962 0.946 ...
##  $ resfvmv...3.: num  0.00789 0.04334 0.03036 0.04885 0.01872 ...
```

```r
head(diff.exp.fmraw) 
```

```
##                    gene.no   fm.mean    m.mean fold.change resfvmv...3.
## ENSRNOG00000000466     224 13.630361 11.834890   1.1517101  0.007886797
## ENSRNOG00000000711     383 10.553780 12.287086   0.8589327  0.043340247
## ENSRNOG00000000720     388  2.755691  3.208746   0.8588061  0.030361329
## ENSRNOG00000000894     491 17.013327 17.689166   0.9617936  0.048849376
## ENSRNOG00000000989     555 18.985484 20.073723   0.9457879  0.018716211
## ENSRNOG00000000991     557  5.554510  4.285393   1.2961496  0.020799509
```

```r
colnames(diff.exp.fmraw) <- c("gene.no", "fm.mean", "m.mean", "fold.change", "raw.pvalue") 
diff.exp.fmraw1 <- diff.exp.fmraw[1:10,] #10 DE genes
knitr::kable(xtable(diff.exp.fmraw1), digits=3, caption = "Differentially expressed genes based on raw p-value")
```



Table: Differentially expressed genes based on raw p-value

                      gene.no   fm.mean   m.mean   fold.change   raw.pvalue
-------------------  --------  --------  -------  ------------  -----------
ENSRNOG00000000466        224    13.630   11.835         1.152        0.008
ENSRNOG00000000711        383    10.554   12.287         0.859        0.043
ENSRNOG00000000720        388     2.756    3.209         0.859        0.030
ENSRNOG00000000894        491    17.013   17.689         0.962        0.049
ENSRNOG00000000989        555    18.985   20.074         0.946        0.019
ENSRNOG00000000991        557     5.555    4.285         1.296        0.021
ENSRNOG00000001032        580    26.144   24.545         1.065        0.026
ENSRNOG00000001172        681    48.957   51.048         0.959        0.013
ENSRNOG00000001313        784    54.500   57.069         0.955        0.022
ENSRNOG00000001375        824     2.083    1.765         1.180        0.033


##Adjust p values using `FDR` (or other) for multiple testing and observe the results


```r
###differentially expression result sheet preparation considering all

diff.exp.fm <- data.frame(countsfmnz[,7:10], resfvmv[,3])
str(diff.exp.fm)
```

```
## 'data.frame':	14046 obs. of  5 variables:
##  $ gene.no     : int  2 5 7 9 10 11 13 16 19 20 ...
##  $ fm.mean     : num  72.95 16.08 1.16 7.25 12.05 ...
##  $ m.mean      : num  75.3 13.6 1.28 7.54 11.92 ...
##  $ fold.change : num  0.969 1.182 0.904 0.962 1.011 ...
##  $ resfvmv...3.: num  0.668 0.722 0.852 0.683 0.869 ...
```

```r
head(diff.exp.fm) 
```

```
##                    gene.no   fm.mean    m.mean fold.change resfvmv...3.
## ENSRNOG00000000007       2 72.946089 75.298589   0.9687577    0.6683702
## ENSRNOG00000000010       5 16.079580 13.600743   1.1822574    0.7219895
## ENSRNOG00000000014       7  1.158026  1.280408   0.9044199    0.8523140
## ENSRNOG00000000021       9  7.253832  7.538195   0.9622770    0.6830687
## ENSRNOG00000000024      10 12.050676 11.922610   1.0107414    0.8692072
## ENSRNOG00000000028      11  1.609237  1.516228   1.0613422    0.5310957
```

```r
colnames(diff.exp.fm) <- c("gene.no", "fm.mean", "m.mean", "fold.change", "raw.pvalue") 

####adjustment in p-value for multiple comparison test: three methods

diff.exp.fm$pvalue.fdr.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "BH"), 4)
diff.exp.fm$pvalue.hb.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "BY"), 4)
#diff.exp.fm$pvalue.bn.adj <- round(p.adjust(diff.exp.fm$raw.pvalue, "bonferroni"), 4)

###checking the genes signifcantly different at 5% level of significance between male anf female

sig.fvmv1 <- diff.exp.fm[diff.exp.fm$pvalue.fdr.adj < 0.05,]
dim(sig.fvmv1)
```

```
## [1] 0 7
```

```r
#View(sig.fvmv1)

sig.fvmv2 <- diff.exp.fm[diff.exp.fm$pvalue.hb.adj < 0.05,]
dim(sig.fvmv2)
```

```
## [1] 0 7
```

```r
#View(sig.fvmv2)
```



#Comments

I am confused at the output found!!! Results after adjusting for multiple comparisons seem horrible for differentially expressed genes. No differentially expressed genes are found!!! 



