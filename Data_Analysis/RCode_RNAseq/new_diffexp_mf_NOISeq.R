setwd("RNASeq_data/new_data_Tony_TPM")
dir()

###TPM
datn <- read.table(file="RNAseq_new_merged_TPM.txt", header = TRUE)
head(datn)
str(datn)

###comparing male and female vehicle

####data preparation

###count data
#can input raw or norrmalized count. here we put rpkm values as count
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
  gather(key = sample, value = TPM, ... = -rowname)

expression_mf_long %>%
  ggplot(aes(TPM, color = sample)) +
  geom_density() +
  scale_x_continuous(trans = "log2") +
  geom_vline(xintercept = 0.01)


###removing rows with all zero counts

countsfmnz <- countsfm[countsfm[,1] != 0 | countsfm[,2] != 0 | countsfm[,3] != 0 | countsfm[,4] != 0 | countsfm[,5] != 0 | countsfm[,6] != 0 , ]
dim(countsfmnz)
#reduced to 26743 genes

###design (factor var) for NOISeq
meta_dat <- read.table(file="sailfish_file_table.txt")
colnames(meta_dat) <- c("SRR ID", "sample.no", "gender", "treatment")
row.names(meta_dat) <- meta_dat$sample.no
meta_dat <- data.frame(meta_dat[,c(1,3,4)])
head(meta_dat)

factorsfm <- data.frame(meta_dat[c(1:3,7:9),])
head(factorsfm)
factorsfm <- droplevels(factorsfm)

###converting data into a NOISeq object
library(NOISeq)
require(NOISeq)
datafm <- NOISeq::readData(data = countsfmnz, factors = factorsfm) #using readData function from NOISeq
datafm

#check what information is included
str(datafm)
head(assayData(datafm)$exprs) ###give the expression data with rpkm counts
head(pData(datafm))
head(featureData(datafm)@data)

#### Quality control of count data
#need not for already normalized count I think

##Generating data for exploratory plots
#can't do so due to lack of information on biotype detection, sequencing depth and expression 
#quantification, sequencing bias detection and batch effect exploration

##count distribution per sample
countsplo <- dat(datafm, factor = NULL, type="countsbio")
explo.plot(countsplo, toplot = 1, samples =NULL, plottype = "boxplot")

###not give much info as all samples are normalized

#sensitivity plot

explo.plot(countsplo, toplot = 1, samples =NULL, plottype = "barplot")
#showing number of features with low counts for each samples 

###Sequencing bias detection

##not needed as the countsa re already normalized

#length bias plot: he “lengthbias” plot describes the relationship between the feature length 
#and the expression values. For each bin, the 5% trimmed mean of the corresponding expression
#values (CPM if norm=FALSE or values provided if norm=TRUE) is computed and depicted in Y axis.

#lengthbiasfm <- dat(datafm, factor = "Vehicle", norm = TRUE, type = "lengthbias")
#explo.plot(lengthbiasfm, samples = NULL, toplot = "global")
#show(lengthbiasfm)
#can't run as we haven't the feature length

####RNA composition

###to check bias cd plot is used

fmcd <- dat(datafm, type = "cd", norm = TRUE, refColumn = 1)

explo.plot(fmcd)

###pca plot

fmPCA = dat(datafm, type = "PCA")
explo.plot(fmPCA, factor = "gender")

###Quality control report

QCreport(datafm, samples = NULL, factor = "gender", norm = TRUE)
#diagnostic test failed.
#result shows normalization is required to correct for bias....confused!!!


###Normalization
#not needed as normalization is done already.


###Low-count filtering
#Excluding features with low counts improves, in general, differential expression results, 
#no matter the method being used, since noise in the data is reduced.

#NOISeq includes three methods to filter out features with low counts: 1. CPM 2. Wilcoxon test 
#3. Proportion test

##Using CPM
fmfilt <- filtered.data(countsfmnz, factor = factorsfm$gender, norm = TRUE, 
                       depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
#Filtering out low count features...
#19028 features are to be kept for differential expression analysis with filtering method 1
#The “Sensitivity plot” described in previous section can help to take decisions on the CPM 
#threshold to use in methods 1 and 3. 

##Using proportion test
#cant do as we haven't sequence depth info

####Differential expression

#NOISeq-real: using available replicates
fmnoiseq1 <- noiseq(datafm, k = 0.5, norm = "n", factor = "gender", nss = 0, replicates = "technical")
head(fmnoiseq1@results[[1]])

fmnoiseq2 <- noiseq(datafm, k = 0.5, norm = "n", factor = "gender", nss = 5, replicates = "biological")
head(fmnoiseq2@results[[1]])

#NOISeqBIO: recommended when we have biological replicates
fmnoiseqbio <- noiseqbio(datafm, k = 0.5, norm = "n", nclust = 50, factor = "gender", 
                         r = 20, adj = 1.5, plot = TRUE, a0per = 0.9, random.seed = 12345, 
                         filter = 1)

######## How to select the differentially expressed features

fmnoiseq1.deg <- degenes(fmnoiseq1, q = 0.8, M = NULL)

fmnoiseq1.deg1 <- degenes(fmnoiseq1, q = 0.8, M = "up")

fmnoiseq1.deg2 <- degenes(fmnoiseq1, q = 0.8, M = "down")

fmnoiseq2.deg <- degenes(fmnoiseq2, q = 0.8, M = NULL)

fmnoiseq2.deg1 <- degenes(fmnoiseq2, q = 0.8, M = "up")

fmnoiseq2.deg2 <- degenes(fmnoiseq2, q = 0.8, M = "down")

fmnoiseqbio.deg <- degenes(fmnoiseqbio, q = 0.95, M = NULL)

fmnoiseqbio.deg1 <- degenes(fmnoiseqbio, q = 0.95, M = "up")

fmnoiseqbio.deg2 <- degenes(fmnoiseqbio, q = 0.95, M = "down")
