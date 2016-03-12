# Results

## Milestone 1 - Data Preprocessing

### Sanity Checking DNA Methylation (Tony)
[1-Methylation_sanity_check.md](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/1-Methylation_sanity_check.md)
* Motivation
  * Performed sanity checks on DNA methylation data to check if anything funny is happening
* Results
  * Discovered that libraries have low coverage, so decided to pool replicates together
  * Noticed that female libraries have less coverage overall, which may influence results
  * Good that Estradiol and Male samples cluster together

### Sanity Checking RNAseq Data (Emma and David) 
[1-RNA_Seq_Sanity_Check.md](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/1-RNA_Seq_Sanity_Check.md)
* Motivation 
 * Same as with methylation data, the RNAseq data must be checked for any obvious anomalies or problems.
* Results 
 * A lot of genes were found to have an RPKM value of 0.
 * Replications in each group were highly correlated with one another. There is one particular sample that looks to be less correlated, although when actually looking at the values 0.97 correlation is good. Based on this, the gene was decided to be kept in. 
