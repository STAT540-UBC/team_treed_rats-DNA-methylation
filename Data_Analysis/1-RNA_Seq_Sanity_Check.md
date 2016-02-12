# 1-RNA_Seq_sanity_checks
Emma  
10 February 2016  



**Processing RNA Seq Data and Sanity Checks**
================================================

**Load Libraries**
-------------------


```r
require(ggplot2)
require(data.table)
require(knitr)
```

**Import Data**
----------------


```r
row.merged <- read.table("C:/Users/Emma/Documents/Masters/STAT540/team_treed_rats-DNA-methylation/RNASeq_data/row.merged.txt", header=TRUE, row.names = 1)
```

The data (rnadata) is arranged in the following table, displaying the gene, rpkm value, as well as the group and replication it belongs to - that is, for example, female vehicle 1, or female z 1. 

```r
head(row.merged)
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

Another data set was made (rpkmlog) with the same data but with log2 rpkm values, as this may be useful for some plots.

**Taking a peek at the data**
------------------------------

Getting a summary of the data set shows that there are 88,548 genes for each of the groups; FVEH, FZEB, MVEH, MZEB. This is good as each of the replicates should have 1/3 of this - 29,516 as expected! 

```r
summary(row.merged)
```

```
##                 genes          rpkm_value           gene.no     
##  ENSRNOG00000000001:    12   Min.   :    0.000   Min.   :    1  
##  ENSRNOG00000000007:    12   1st Qu.:    0.000   1st Qu.: 7380  
##  ENSRNOG00000000008:    12   Median :    0.789   Median :14758  
##  ENSRNOG00000000009:    12   Mean   :   15.194   Mean   :14758  
##  ENSRNOG00000000010:    12   3rd Qu.:    9.759   3rd Qu.:22137  
##  ENSRNOG00000000012:    12   Max.   :12039.940   Max.   :29516  
##  (Other)           :354120                                      
##   replication     group      group.fct   
##  Min.   :1    Min.   :1.00   FVEH:88548  
##  1st Qu.:1    1st Qu.:1.75   FZEB:88548  
##  Median :2    Median :2.50   MVEH:88548  
##  Mean   :2    Mean   :2.50   MZEB:88548  
##  3rd Qu.:3    3rd Qu.:3.25               
##  Max.   :3    Max.   :4.00               
## 
```

Looking at the rpkm values of the genes as a whole, there are a LOT that have a value of 0. Is this bad? Shouldn't every gene have some expression?  


```r
test <- subset(row.merged, rpkm_value == 0)
summary(test)
```

```
##                 genes          rpkm_value    gene.no       replication   
##  ENSRNOG00000000039:    12   Min.   :0    Min.   :    3   Min.   :1.000  
##  ENSRNOG00000000050:    12   1st Qu.:0    1st Qu.:17461   1st Qu.:1.000  
##  ENSRNOG00000000199:    12   Median :0    Median :21800   Median :2.000  
##  ENSRNOG00000000252:    12   Mean   :0    Mean   :20612   Mean   :2.001  
##  ENSRNOG00000000326:    12   3rd Qu.:0    3rd Qu.:25921   3rd Qu.:3.000  
##  ENSRNOG00000000401:    12   Max.   :0    Max.   :29516   Max.   :3.000  
##  (Other)           :116184                                               
##      group       group.fct   
##  Min.   :1.000   FVEH:28801  
##  1st Qu.:2.000   FZEB:29378  
##  Median :2.000   MVEH:28593  
##  Mean   :2.505   MZEB:29484  
##  3rd Qu.:4.000               
##  Max.   :4.000               
## 
```

It also needs to be checked for each group, how many have a readout of 0, as this may cause some problems / bias.

The mean values of rpkm are shown below:


```r
meanrpkms <- aggregate(row.merged[, 2], list(row.merged$group.fct), mean)
meanrpkms
```

```
##   Group.1        x
## 1    FVEH 15.24095
## 2    FZEB 15.17979
## 3    MVEH 15.19049
## 4    MZEB 15.16369
```


