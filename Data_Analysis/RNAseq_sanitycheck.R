## Load Libraries
require(ggplot2)
require(plyr)
require(reshape2)
require(dplyr)
require(knitr)

## Load data set
row.merged <- read.delim("C:/Users/David/team_treed_rats-DNA-methylation/RNASeq_data/row.merged.txt")

## Check data
head(row.merged)
summary(row.merged)

## Make log2 data set
rpkmlog <- row.merged
rpkmlog[, 2] <- log(rpkmlog[2], 2)

## Analyzing row.merged for no expression values
rpkm.null <- subset(row.merged, rpkm_value == 0)
summary(rpkm.null)

## data frame "rpkm.null" contains all gene observations with no expression
rpkm.null  <- rpkm.null[ order(rpkm.null[,1], rpkm.null[,5]),]

## Organized by genes by group and the frequency they occur, then subsetted freq=3
rpkm.summary <- table(rpkm.null[, c('genes', 'group.fct')], useNA='ifany')
rpkm.freq <- as.data.frame(rpkm.summary)
rpkm.freq <- subset(rpkm.freq, Freq==3)
rpkm.freq <- dcast(rpkm.freq, genes~group.fct)
  
## perhaps before here, need to setup so that the gene is followed by columns of its freq
## for each group condition for easier comparison
rpkm.comp <- dcast(rpkm.null, genes~group.fct)

## Analyzing expressed data subset
expressed <- subset(row.merged, rpkm_value > 0.00001)
summary(expressed)

## Taking the mean of the replicates and grouping data by gene+group
expressed.mean <- expressed %>% group_by(genes,group.fct) %>% summarize(rpkm.mean = mean(rpkm_value))

## Transposing data so columns are the group.fct's
expressed.mean <- dcast(expressed.mean, genes~group.fct)
