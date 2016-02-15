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
test <- subset(row.merged, rpkm_value == 0)
summary(test)

## data frame "test" contains all gene observations with no expression
test.sorted  <- test[ order(test[,1], test[,5]),]

## Organized by genes by group and the frequency they occur, then subsetted freq=3
test.summary <- table(test[, c('genes', 'group.fct')], useNA='ifany')
test.freq <- as.data.frame(test.summary)
test.freq3 <- subset(test.freq, Freq==3)

## perhaps before here, need to setup so that the gene is followed by columns of its freq
## for each group condition for easier comparison

