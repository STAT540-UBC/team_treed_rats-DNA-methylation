This is the private repository of team **TREED** for their group project as part of STAT 540/ BIOF 540/ GSAT 540.

**Epigenetic Determinants of Gender in Rats**
==============================================

This is the repo for group project of **TEAM TREED**

**Table 1. Our team members:**

Github ID |  Name
---------|------------
[@hui-tony-zk](https://github.com/hui-tony-zk) | Tony Hui (MSc, Genome Science and Tech - Bioinformatics)
[@RashedHUBC](https://github.com/RashedHUBC) |	Md. Rashedul Hoque (MSc, Bio-Statistics)
[@emminic93](https://github.com/emminic93) |	Emma Titmuss (MSc, Genome Science and Tech)
[@eclaks](https://github.com/eclaks) | Emma Laks (MSc, Genome Science and Tech)
[@david-rattray](https://github.com/David-Rattray) |	David Rattray (MSc, Biochemistry)

Our project is based on the paper ["Brain feminization requires active repression of masculinization via DNA methylation, Nugent et al 2015"](http://www.nature.com/neuro/journal/v18/n5/full/nn.3988.html).

Our [project proposal](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/originalproposal.md) gives the outline of the project.

**Introduction**
--------------
Genetic sex (XX vs. XY) has been held as the dominant sexual differentiation model, causing differentiation of the gonads, which secrete sex hormones, such as estradiol, to masculinize the brain. However, recent evidence suggests that environmental and epigenetic influences also contribute to sex differences.[[1,2,3]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) Enzymes such as methyltransferases influence the epigenome via the methylation of the genetic code. Male rats have lower DNA (cytosine-5)-methyltransferase 3A (DNMT3a) activity and DNA methylation than females.[[1]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) Nugent et al demonstrated using DNMT inhibitors or conditional knockouts that female rats display masculinized behavior.[[1]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) We used the raw reads generated from Nugent et al to find differentially expressed genes (DEGs) and differentially methylated regions (DMRs) between all conditions to identify gender-related genes that are epigenetically regulated.

<img src="https://raw.githubusercontent.com/STAT540-UBC/team_treed_rats-DNA-methylation/master/Background_Research/rats_figure.png?token=APx5yAMWHcgkQaDqsYeneu6Y2M4TN2-Sks5XEG_CwA%3D%3D" width="450px" height="308px" />

**Figure 1. Rat experiments carried out by Nugent et al to demonstrate behavioural changes in female rats upon treatment with estradiol or zebularine (zeb), a DNMT inhibitor.**

**Motivation**
---------------
The paper by Nugent et al demonstrates that: 
* Male rats have lower DNMT3a activity and DNA methylation than females.
* Inhibiting DNMT masculinized neuronal markers and sexual behaviour in female rats, and has no effect on male behaviour. Females with conditional knockout of isoform DNMT3 also display male sexual behaviour.
* Even outside of the restricted period of development where they are sensitive to hormone therapy, DNMT inhibition still masculinizes females.
* Some changes in gene expression are as a result of DNMT inhibition.
  + The authors did not, however, attempt to link these changes to their methylation data which is why we have selected to further evaluate this paper.

**Objective**
----------
Find overlaps between differentially expressed genes (DEGs) and differentially methylated genes (DMRs) (figure 2A) to reveal potential epigenetically-regulated genes involved in masculinization and feminization of the rat brain by looking for overlapping genes in samples (figure 2B).

![](https://raw.githubusercontent.com/STAT540-UBC/team_treed_rats-DNA-methylation/master/Background_Research/figure_dmr_deg.png?token=APx5yDnzBXvjeTY4m0jLSqEuP0H5-oqYks5XEGj7wA%3D%3D)
**Figure 2. Strategy to (A) identify DEG/DMR overlap and (B) identify genes of interest for gender regulation.**

**The Data**
--------------
 
**Table 2. The data provided from the [paper](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66203) contains [RNA-seq](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data) - processed data with gene name and RPKM (Reads Per Kilobase of transcript per Million) and [Whole Genome Bisulfite Sequencing (WGBS)](http://www.ncbi.nlm.nih.gov/bioproject/?term=275796) - DNA-methylation (~270 million raw reads), both taken from the preoptic area of the rat brain.**

 Sample |  [RNASeq](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data/new_data_Tony_TPM) (Day 2) | [WGBS](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/methylation_data) (Day 4)
--------|-----------------|-------------
male | 3 replicates | 3 replicates merged into 1
male treated with zebularine |	3 replicates | -
female |	3 replicates | 3 replicates merged into 1
female treated with zebularine  | 3 replicates | -
female treated with estradiol |	- | 3 replicates merged into 1


**Analysis and tasks directory** 
---------------------------------------
[Methodology](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/methods.md) with division of labour, and links to analysis and issues:

1. Firstly, align WGBS reads and call methylation with bismark  (**Tony**)
    - [click for issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/1)
    - also had to chase down authors for correct BAM file - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/5)
2. At the same time… 
 + Find differentially methylated regions (DMR) between male vs female, and find nearest gene for each DMR using HOMER.[[4]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) (**Tony**)
     + [results (with code)](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/2-Calling_DMRs.md)
     + [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/8)
 + RNA alignment with SAILFISH.[[5]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) (**Tony**)
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/24)
 + RNA-seq data file preparation for sanity check and analysis. (**Rashed**)
     - [processed data](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data)
 + RNA sanity checks (**Emma T., David**)
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/17)
     - [results](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/1-RNA_Seq_Sanity_Check.md) 
 + Background research on differentially expressed genes methods (**Rashed**)
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/3)
 + Find differentially expressed genes between male vs female (**Tony, Rashed**) 
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/25)
     - [results](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/3-TPM_RNA_Seq_differential__expression.md)
     - [m/f DEG list](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/RNAseq_result/DE_genes/glmQLFit_DE_genes.tsv). 
3. Generate a list of gene and region pairs by finding overlap between 2a and 2b. This list of gene/region pairs represent epigenetically regulated genes that are important for gender  (**Tony**)
    - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/18)
    - [list of DEGs between f/m close to a convergent DMR](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/methylation_data/3-all_gene_associated_DMRs.tsv)
4. At the same time...
 + Find DMRs between female and female + estradiol, and only consider DMRs that overlap or are 150Kb or less distance away from transcriptional start sites of genes in (3)[[6]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md).  (**Tony**)
     - [result](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/2-Calling_DMRs.md)
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/29)
 + Find differentially expressed genes between female and female + DNMT inhibitor and only consider genes that overlap with (3) (**Tony**)
     - [results](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/3.1-DE_genes_femaleVSzeb.md)
     - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/28)
     - [list of f/z-f DEGS](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/RNAseq_result/DE_genes/3.1-femVSfemZeb_allDE_genes.tsv.) 
     - [list of f/fz&m DEGs](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/RNAseq_result/DE_genes/3.1-femVSfemZeb_glmQLFit_DE_genes.tsv)
5. Select genes that overlap between 4a and 4b to form list of master list of gender genes that can be artificially altered via epigenetic reprogramming (**Tony**)
    - [issue](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/28)
    - [list of DEGs between f/m and f/fz close to a convergent DMR](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/methylation_data/3-DE_gene_associated_DMRs.tsv).
6. To conclude, and properly answer the biological question we will describe the biological relevance of genes in the master list in (5).
 + [Background research into rats, gender, epigenetics.](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/16).  (**Emma L.**)
 + [For known genes, describe gene pathways and functions](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/19) (using [GSEA](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/Background_Research/GSEA), [UNIPROT](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/Background_Research/annotation_Uniprot_lit)). (**Emma L**)
7. Administrative tasks and preparation of deliverables:
 + [Project proposal](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/6). (**Everyone**)
 + Repo readme directories with arrangement. (**Emma L., Emma T, Tony, Rashed**)
 + [Progress report](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/20). (**Emma T.**)
 + [Poster](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/issues/23). (**Emma L., Emma T., Tony, Rashed**)
 + [Bibliography](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md). (**Emma L.**)
8. Methodological research for explanation of the various statistical tools and methods used in the project. (**Rashed**)

**[Results](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/README.md)** 
---------------------------------------

RNASeq data was provided in RPKM and it was also aligned to an outdated version 4 (2004) of the rat genome. Only one gene was found using this for DEG analysis, so read counts were regenerated by `SAILFISH` for further analysis. For DEG analysis, methods used in `R` packages `edgeR` (`glmFit` & `glmQLFit` functions), `NOISeq`, `DESeq` and `limma` were all tested. `glmQLFit` (quasi-likelihood negative binomial model)  was selected to maximize sensitivity even at the cost of false discovery rate, as it gave the most DEGs, and better explained the data variation by incorporating two dispersions. The use of `glmQLFit` was further supported based on its overlap with DEGs detected by other methods. 

Following RNAseq DEG analysis, 164 genes were differentially expressed between male and female rats and of these, 53 genes were found to be differentially expressed between females and zeb-treated females following the same masculinizing gene expression pattern.

WGBS sanity checks revealed significantly fewer reads for one female sample files hosted on [NCBI SRA](http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP055171); the first author updated the data on request. The female rats also had a lower number of WGBS reads relative to the male and estradiol-treated females, and the WGBS data had low coverage overall, so the replicates for this dataset were merged following protocol from the literature.[[7]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) Using the merged data, 16308 DMRs were found between female and male rats in total, and after filtering DMRs that are promoters of other genes and those that are >150kb away from the nearest promoter of DEG genes of interest, 267 DMRs remain as putative gender regulating regions. 21465 DMRs were found between female and estradiol, and 100 of these overlapped with the DMRs found in females vs. males. The methylation of 99% of these regions were the same in male and estradiol-treated females. Of these 100 regions, 13 associated with genes that were differentially expressed across females vs. males and zeb-treated females. These are our ten potential epigenetically regulated genes involved in sex determination. 

**Table 3. Putative epigenetically regulated genes involved in sex determination found through DEG/DMR overlap analysis.**

gene | description | DMR annotation | distance to transcriptional start site | influence
-----|-------------|----------------|----------------------------------------|----------
Gipr |	gastric inhibitory polypeptide receptor	| promoter-TSS	| 456	| feminizing
LOC100362027	| ribosomal protein L30-like	| intron	| -116711	| feminizing
Fbxw10	| F-box and WD repeat domain containing 10	| exon	| -110780	| feminizing
Tcf12	| transcription factor 12 |	intergenic	| 120886	| feminizing
Plch1	| phospholipase C, eta 1 | intergenic	| 107433	| feminizing
| | | intergenic	| 106442	| feminizing		
| | | intergenic	| -143383	| feminizing		
Vars |	valyl-tRNA synthetase	| intron |	-93035	| feminizing
Foxj2	| forkhead box J2	| intron	| -54560	| feminizing
Hsp90aa1	| heat shock protein 90aa1	| intergenic	| 103431	| masculinizing
| | | intron |	80776	| masculinizing		
Adcy6	| adenylate cyclase 6	| intron	| -138546	| masculinizing
Tpp2	| tripeptidyl peptidase II	| intron |	-66157	| masculinizing

**Discussion** 
---------------------------------------

Our study presents a method for analysis of WGBS and RNASeq data to integrate DMRs and DEGs, and uses it to identify ten putative genes that determine sexual differentiation and behaviour and can be epigenetically regulated. Two of these genes, TCF12 and ADCY6, are already known to have have biological significance to neurology.[[8,9,10]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) TCF12, a transcriptional regulator involved in the initiation of neuronal differentiation was found to be both hypomethylated and upregulated in females relative to to males and masculinized females, suggesting it could be a epigenetically regulated feminizing gene.[[8]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) ADCY6 was both upregulated and associated with a hypomethylated DMR in males and masculinized females compared to vehicle females, making it a putative epigenetically regulated masculinizing gene. ADCY6 is an adenylate cyclase; this family of proteins catalyzes the formation of the signaling molecule cAMP and is known to be involved in intracellular signalling,[[9]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) and in developmentally and neurophysiologically relevant signaling pathways, including neural plasticity.[[10]](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md) Our findings provide insight into the mechanisms of gender, and may help explain non-heterosexual gender identities in mammals. 
 
**Useful links**
----------------
 + [RNASeq data](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/RNASeq_data/new_data_Tony_TPM)
 + [WGBS data](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/tree/master/methylation_data)
 + [Methods](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/methods.md)
 + [Result outputs](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Data_Analysis/README.md)

**Deliverables** 
---------------------------------------
 + [Project proposal](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/originalproposal.md)
 + [Progress report](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/Progress_Report.md)
 + [Poster](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/poster_treed.pdf)
 + [Bibliography](https://github.com/STAT540-UBC/team_treed_rats-DNA-methylation/blob/master/Final%20Deliverables/bibliography.md)
