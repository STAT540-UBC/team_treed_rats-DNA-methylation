This is the private repository of team **TREED** for their group project as part of STAT 540/ BIOF 540/ GSAT 540.

**Epigenetic Determinants of Gender in Rats**
==============================================

This is the repo for group project of **TEAM TREED**

Our team members are: 

Github ID |  Name
---------|------------
[@hui-tony-zk](https://github.com/hui-tony-zk) | Tony Hui (MSc, Genome Science and Tech - Bioinformatics)
[@RashedHUBC](https://github.com/RashedHUBC) |	Md. Rashedul Hoque (MSc, Bio-Statistics)
[@emminic93](https://github.com/emminic93) |	Emma Titmuss (MSc, Genome Science and Tech)
[@eclaks](https://github.com/eclaks) | Emma Laks (MSc, Genome Science and Tech)
[@david-rattray](https://github.com/David-Rattray) |	David Rattray (MSc, Biochemistry)

Our project is based on the paper ["Brain feminization requires active repression of masculinization via DNA methylation, Nugent et al 2015"](http://www.nature.com/neuro/journal/v18/n5/full/nn.3988.html).

Our [project proposal](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/originalproposal.md) gives the outline of the project.

**Abstract**
--------------
Male rats have lower DNA (cytosine-5)-methyltransferase 3A (DNMT3a) activity and DNA methylation than females. Nugent et al demonstrated using DNMT inhibitors or conditional knockouts that female rats display masculinized behavior. We will analyze their RNA-seq data from male and female rats with and without DNMT inhibitor treatment to relate the involvement of gene expression to differentially methylated regions found in females with and without estradiol therapy and male rats using their Whole Genome Bisulfite Sequencing (WGBS) data. We will conduct some exploratory analysis on both data through appropriate graphical visualization. After that to compare the different gene expressions and epigenetic characteristics we will use t-tests for two-groups (explicitly male and female) and ANOVA or analysis of variance (for more than two groups). Also if data permits some linear models will be checked for fitting to express the relationship in a formal way. 

**Motivation**
---------------
The paper by Nugent et al demonstrates that: 
* Male rats have lower DNMT3a activity and DNA methylation than females
* Inhibiting DNMT masculinized neuronal markers and sexual behaviour in female rats, and has no effect on male behaviour. Females with conditional knockout of isoform DNMT3 also display male sexual behaviour  
* Even outside of the restricted period of development where they are sensitive to hormone therapy, DNMT inhibition still masculinizes females 
* Some changes in gene expression are as a result of DNMT inhibition.
  + The authors did not, however, attempt to link these changes to their methylation data which is why we have selected to further evaluate this paper.

**Objectives**
----------
Based on the data already collected from the paper we aim to:

1. Find epigenetically regulated genes that determine gender outcome (by comparing male and female)
2. For each epigenetically regulated gene discovered in (1), determine which are changed by estradiol (testosterone) or DNMT inhibitor.

**The Data**
--------------
The data provided from the [paper](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66203) contains:
* [RNA-seq](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data) - processed data with gene name and RPKM (Reads Per Kilobase of transcript per Million)   
 + Female x3
 + Female + DNMT inhibitor x3
 + Male x3
 + Male + DNMT inhibitor x3
* [Whole Genome Bisulfite Sequencing (WGBS)](http://www.ncbi.nlm.nih.gov/bioproject/?term=275796) - DNA-methylation (~270 million raw reads) -  
 + Male x3
 + Female x3
 + Female + estradiol x3

**Methodology and Division of Labour** 
---------------------------------------

1. Firstly, [align WGBS reads and call methylation with bismark](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/1) (**Tony**)
2. At the same time… 
 + [Find differentially methylated regions (DMR) between male vs female, and find nearest gene for each DMR](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/8) (Tony)
 + [RNA alignment with SAILFISH](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/24) (**Tony**)
 + [RNA-seq data file preparation for sanity check and analysis](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data) (**Rashed**)
 + [RNA sanity checks] (https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/1-RNA_Seq_Sanity_Check.md) -[issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/17) (**Emma T., David**)
 + [Background research on differentially expressed genes between male vs female] (https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/3) (**Rashed**)
 + [Find differentially expressed genes between male vs female](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/25) (**Tony, Rashed**) 
3. Generate a list of gene and region pairs by finding overlap between 2a and 2b. This list of gene/region pairs represent epigenetically regulated genes that are important for gender (**Tony**)
4. At the same time...
 + [Find DMRs between female and female + estradiol, and only consider DMRs that overlap or are close to regions in (3)](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/29).  (**Tony**)
 + [Find differentially expressed genes between female and female + DNMT inhibitor and only consider genes that overlap with (3)](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/28). (**Tony**)
5. Select genes that overlap between 4a and 4b to form list of master list of gender genes that can be artificially altered via epigenetic reprogramming (**Tony**)
6. To conclude, and properly answer the biological question we will describe cellular pathways in the master list in (5).
 + [Background research into rats, gender, epigenetics.](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/16)  (**Emma L.**)
 + [For known genes, describe gene pathways and functions](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/19) (using [GSEA](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/Background_Research/GSEA), [UNIPROT](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/Background_Research/annotation_Uniprot_lit)) (**Emma L**)
7. Administrative tasks and preparation of deliverables:
 + [Project proposal](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/6) - (**Everyone**)
 + Repo readme directories with arrangement (**Emma T., Emma L, Tony, Rashed**)
 + [Progress report](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/20) (**Emma T.**)
 + [Poster](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/23) (**Emma L., Emma T., Tony**)
 + [Bibliography](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) (**Emma L.**)
8. Methodological research for explanation of the various statistical tools and methods used in the project (**Rashed**)
