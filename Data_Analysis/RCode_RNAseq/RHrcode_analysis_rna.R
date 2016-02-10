
rm(list=ls())

###read the merged data
asd <- read.table(file="merged.txt", header=T)
head(asd)
dim(asd)
str(asd)

#read the row.merged data
asd1 <- read.table(file="row.merged.txt", header=T)
str(asd1)
head(asd1)

####read the p-value data
bdat <- read.table(file = "anova.pval.txt", header = T)
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

tbl_anovdat <- tbl_df(bdat1)