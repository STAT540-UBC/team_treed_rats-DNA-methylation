Methods
-------

**Whole genome bisulphite sequencing analysis:**

+ Aligned and called CpG methylation using Bismark.
+ Methylation calls are smoothed across the genome using BSmooth - a local likelihood estimator conceptually similar to loess smoothing or running average. 
+ Differentially methylated regions between male and female were determined and the nearest gene to the DMR was determined using HOMER.[[4]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) 

**Sanity checking RNASeq data:**

+ Pearson correlated each sample’s normalized read counts to determine if any of the samples appeared to be outliers. 
+ Correlations were visualized using a heatmap.
+ A density plot was made of the samples to ensure there were no unexpected spikes in expression that could have resulted from a technical error and affect the DEGs.

**RNASeq analysis:**

+ SAILFISH was used to estimate isoform abundances from the RNASeq reads and reference sequences.[[5]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md)
+ SAILFISH uses an alignment free algorithm to estimate the abundances without complexity of read mapping (instead using k-mer indexing and counting), making this a very fast and reliable tool. 

**Find differentially expressed genes between male vs. female:**

+ A number of R packages were used to find DEGs between male and female; edgeR, limma, DESeq and NOISeq. 
+ The main approach was the use of glmQLFit in edgeR which addresses two types of dispersion; the gene specific dispersion modelled by a QL parameter, and the other is the global NB parameter over all of the genes. 

**DMR and DEG overlap analysis:**

+ When overlapping the DMRs and DEGs, the distance was limited to 150kbs as this has been described as the distance chromatin loops occur most frequently.[[6]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md)
+ Each gene with higher expression in female, and the associated DMR is less methylated in female, was considered to be epigenetically regulated. 
