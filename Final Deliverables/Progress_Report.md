# Progress Report
Emma  
2 March 2016  



**Progress Report on Data Analysis**
=====================================

As stated in the project proposal, our goals were: 

* Find epigenetically regulated genes that determine gender outcome (by comparing male and female)
* For each epigenetically regulated gene discovered in (1), determine which are changed by estradiol (testosterone) or DNMT inhibitor.

**Tasks completed so far:** 

* Complete sanity checks on methylation data.
    * It was discovered here that the data had very poor coverage and so the decision was made to combine the replicates to improve the data. 
* Complete sanity checks on RNA Seq data. 
    * The replicates in the RNA data were of better quality than methylation. Groups were highly correlated with one another, so samples have not been merged, nor any removed.
* Find differentially methylated regions (DMRs) between males and females.
    * bssmooth package was used to find DMRs between these two groups. 
* Find differentially expressed genes (DEGs) between males and females. 
    * This was carried out using a multiple t test and fold change, following checking the distributions of rpkm values over male and female samples. T test and P-value adjustment was completed using FDR. 
    
**Research:**

Further research into genes that are involved in gender is also being undertaken so that results found from the analysis can be compared to other literature to assess reliability. 

Work on finding DEGs appears to follow a similar analysis pipeline as ours, which is encouraging! 

**Current Work**

Currently, the work is focused on aligning the DMRs to neaby genes. Once this has been completed, this list will be compared with the differentially expressed genes from the RNA Seq data to generate a list of genes that are epigenetically different between male and females (**Tony and Rashed**). 

These genes will then be annotated using many databases (**Emma L, Emma T, David**).
  
  
**Future Tasks:**
The next stages of the project are to find DMRs between female and female + estradiol and to find differentially expressed genes between female and female + DNMT inhibitor. 

This will be a lot quicker than the comparison between male and female as the analysis pipeline has now been used and refined. 

The final stages will be to find overlap between the two sets to create a 'master' list of 'gender genes', which will then also be annotated to allow cellular pathway descriptions.  
