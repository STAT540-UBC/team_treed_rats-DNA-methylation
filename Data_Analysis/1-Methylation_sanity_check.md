# Processing DNA Methylation Data

# Processing and sanity checking DNA methylation data

Warning, the data files come out to be 2.5 GB! Check that your computer has more than 2.5 GB of RAM, or else your computer will crash!

## Load libraries


```r
require(data.table)
require(foreach)
require(doMC)
require(knitr)
require(dplyr)
require(ggplot2)
require(pheatmap)
# setwd("Data_Analysis")
registerDoMC(9)
```

## Load data


```r
cpg_file_names <- dir("../methylation_data/", "^SRR.*.CG.bed", full.names = T)

file_table <- read.delim("../methylation_data/sample_table.txt", header = FALSE) %>% arrange(V1)

registerDoMC(9)
cpg_files <- foreach(i=seq_along(cpg_file_names), .combine = c) %dopar% {
  tmp <- fread(input = cpg_file_names[i]) %>% 
    setnames(c("chr","pos","end","converted_C","cov")) %>% 
    filter(chr != "chrMT") %>%
    select(-end)
#     mutate(sample = file_table$V2[i]) %>%
#     mutate(sampel = gsub("[1-3]", "", sample)
  # save(tmp, file = paste0("../methylation_data/",file_table$V2[i],".RData"), compress = TRUE)
  list(tmp)
}

names(cpg_files) <- file_table$V2
```

## Merge replicates


```r
samples <- gsub("[1-3]", "", names(cpg_files)) %>% unique()

samples_replicates <- foreach(i=samples, .combine=c) %do% {
  list(paste0(i, 1:3))
} 

registerDoMC(9)
cpg_files_merged <- foreach(i=seq_along(samples_replicates), .combine = c) %do% {
  tmp <- rbindlist(cpg_files[samples_replicates[[i]]]) %>% group_by(chr, pos) %>% summarize(converted_C = sum(converted_C), cov = sum(cov))
  list(tmp)
}

names(cpg_files_merged) <- samples

# for (i in names(cpg_files_merged)) { 
#   tmp <- cpg_files_merged[[i]]
#   save(tmp, file = paste0("../methylation_data/", i, ".RData"), compress = TRUE)
#   }

rm(cpg_files)
```

## Analysis begins

### Number of CpGs

#### Per sample


```r
lapply(cpg_files_merged, nrow) %>% data.frame() %>% t() %>% data.frame() %>% add_rownames() %>% arrange(rowname)
```

```
## Source: local data frame [3 x 2]
## 
##     rowname        .
##       (chr)    (int)
## 1 estradiol 19313755
## 2    female 18172632
## 3      male 20194551
```

#### Per chromosome per sample


```r
cpg_per_chrom_per_sample <- foreach(i=1:length(cpg_files_merged), .combine = rbind) %dopar%{
  cpg_files_merged[[i]] %>% group_by(chr) %>% tally() %>% mutate(sample = names(cpg_files_merged)[i])
} 

cpg_per_chrom_per_sample %>% 
  mutate(facet = gsub("[1-3]", "", sample),
         chr = factor(chr, levels = paste0("chr", c(1:20, "X", "Y")))) %>%
  ggplot(aes(chr, n, color = facet, group = sample)) +
  # facet_wrap(~facet) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](1-Methylation_sanity_check_files/figure-html/cpg_per_chrom_per_sample-1.png)

Well that sucks. Looks like the female libraries have less reads overall. Luckily, the coverage seems even.

### Average coverage per CpG


```r
foreach(i=1:length(cpg_files_merged), .combine = rbind) %dopar%{
  cpg_files_merged[[i]] %>% ungroup() %>% summarize(average_Cov = sum(cov)/n()) %>% mutate(sample = names(cpg_files_merged)[i])
} %>% arrange(sample)
```

```
##    average_Cov    sample
## 1:    5.304188 estradiol
## 2:    3.287158    female
## 3:    7.130281      male
```

### Coverage distribution


```r
cov_distribution <- foreach(i=1:length(cpg_files_merged), .combine = rbind) %dopar%{
  cpg_files_merged[[i]] %>% group_by(cov) %>% tally() %>% mutate(sample = names(cpg_files_merged)[i])
}

cov_distribution %>%
  mutate(facet = gsub("[1-3]", "", sample)) %>%
  ggplot(aes(cov, n, group = sample, color = facet))+
  geom_line() +
  xlim(0,30)
```

```
## Warning: Removed 493 rows containing missing values (geom_path).
```

![](1-Methylation_sanity_check_files/figure-html/cov_distribution-1.png)

Distribution of coverage looks poisson!

### Average methylation per sample


```r
foreach(i=1:length(cpg_files_merged), .combine = rbind) %dopar%{
  cpg_files_merged[[i]] %>% ungroup %>% summarize(mean_meth = mean(converted_C/cov) %>% round(3)) %>% mutate(sample = names(cpg_files_merged)[i])
} %>% arrange(sample)
```

```
##    mean_meth    sample
## 1:     0.780 estradiol
## 2:     0.775    female
## 3:     0.776      male
```

Looks like methylation is generally very similar across samples (xxcept female 1, what's up with that?)

### Distribution of DNA methylation per sample


```r
dist_DNA_meth <- foreach(i=1:length(cpg_files_merged), .combine = rbind) %dopar%{
  tmp <- cpg_files_merged[[i]]
  tmp$meth <- tmp$converted_C/tmp$cov
  tmp <- hist(tmp$meth, plot=FALSE, breaks = seq(0,1,0.05))
  data.table(
    breaks = tmp$breaks[-1],
    counts = tmp$counts 
  ) %>% mutate(sample = names(cpg_files_merged)[i])
}

dist_DNA_meth %>%
  group_by(sample) %>%
  mutate(total = sum(counts),
         freq = counts/total) %>%
  mutate(facet = gsub("[1-3]", "", sample)) %>%
  ggplot(aes(x = breaks-0.025, y = freq, group = sample, color = facet)) +
  geom_line()
```

![](1-Methylation_sanity_check_files/figure-html/dist_DNA_meth-1.png)

Yep, super low coverage = higher chance of 0% or 100%...

### Pairwise correlation of DNA methylation


```r
all_cpgs <- rbindlist(cpg_files_merged) %>% select(chr, pos) %>% unique()

cpg_files_merged_cbind <- foreach(i=1:length(cpg_files_merged), .combine = cbind) %do% {
  file_name <- names(cpg_files_merged)[[i]]
  tmp <- cpg_files_merged[[i]] %>% mutate(meth = converted_C/cov)
  left_join(all_cpgs, tmp) %>%
    select(meth) %>%
    setnames(paste0(file_name, c("_meth")))
}
```

```
## Joining by: c("chr", "pos")
## Joining by: c("chr", "pos")
## Joining by: c("chr", "pos")
```

```r
cpg_files_merged_cbind <- cbind(all_cpgs, cpg_files_merged_cbind)
```


```r
pairwise_cor <- cpg_files_merged_cbind %>%
  select(contains("meth")) %>%
  cor(use = "pairwise.complete.obs")

diag(pairwise_cor) <- 1

hclust <- hclust(as.dist(1-pairwise_cor))

pheatmap(pairwise_cor, cluster_rows = hclust, cluster_cols = hclust, display_numbers = T)
```

![](1-Methylation_sanity_check_files/figure-html/pairwise_cor-1.png)

Correlations look much better now



