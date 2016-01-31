### Workflow for obtaining DNA methylation calls from FASTQ

1. Downloaded FastQ (single-end reads) from SRA using fastq-dump v2.4.2

| Accession       | Sample         |
|------------|------------|
| SRR1810939 | male3      |
| SRR1810938 | male2      |
| SRR1810937 | male1      |
| SRR1810936 | estradiol3 |
| SRR1810931 | female2    |
| SRR1810930 | female1    |
| SRR1810934 | estradiol2 |
| SRR1810933 | estradiol1 |
| SRR1810932 | female3    |

2. Downloaded bismark v 0.14.5 - http://www.bioinformatics.babraham.ac.uk/projects/bismark/
3. Downloaded rattus_norvegicus genome from Ensemble - ftp://ftp.ensembl.org/pub/release-83/fasta/rattus_norvegicus/dna/
  * Rattus_norvegicus.Rnor_6.0.dna_rm.toplevel.fa.gz - Hard-masked sequence
4. Followed instructions for Bismark http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.14.4.pdf
  1. Prepared rat genome with bismark_genome_preparation (default settings)
  2. 
