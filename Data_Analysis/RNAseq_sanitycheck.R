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

## Organize by group??