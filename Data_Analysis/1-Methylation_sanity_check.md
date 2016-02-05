# Processing DNA Methylation Data

# Processing and sanity checking DNA methylation data

Warning, the data files come out to be 2.5 GB! Check that your computer has more than 2.5 GB of RAM, or else your computer will crash!

## Load libraries


```r
require(data.table)
require(foreach)
require(doMC)
require(knitr)
# opts_knit$set(root.dir = '../..')
require(dplyr)
require(ggplot2)
require(pheatmap)
```

## Load data


```r
cpg_file_names <- dir("../methylation_data/", "*.CG.bed", full.names = T)

registerDoMC(9)
cpg_files <- foreach(x=cpg_file_names, .combine = c) %dopar% {
  tmp <- fread(input = x) %>% 
    setnames(c("chr","pos","end","converted_C","cov")) %>% 
    filter(chr != "chrMT") %>%
    mutate(meth = (converted_C/cov) %>% round(3)) %>%
    select(-end, -converted_C)
  list(tmp)
}

file_table <- read.delim("../methylation_data/sample_table.txt", header = FALSE) %>% arrange(V1)

names(cpg_files) <- file_table$V2
```


```r
all_cpgs <- rbindlist(cpg_files) %>% select(chr, pos) %>% unique()

cpg_files_merged <- foreach(i=1:length(cpg_files), .combine = cbind) %dopar% {
  file_name <- names(cpg_files)[[i]]
  left_join(all_cpgs, cpg_files[[i]]) %>%
    select(meth, cov) %>%
    setnames(paste0(file_name, c("_meth", "_cov")))
}

cpg_files_merged <- cbind(all_cpgs, cpg_files_merged)

# save("cpg_files_merged", file = "methylation_data/methylation_data.RData", compress = TRUE)
```

## Analysis begins

### Number of CpGs

#### Per sample


```r
lapply(cpg_files, nrow) %>% data.frame() %>% t() %>% data.frame() %>% add_rownames()
```

```
## Source: local data frame [9 x 2]
## 
##      rowname        .
##        (chr)    (int)
## 1    female1  8922668
## 2    female2  8702855
## 3    female3 13632672
## 4 estradiol1 13811801
## 5 estradiol2 14178594
## 6 estradiol3 16619149
## 7      male1 18058509
## 8      male2 15173572
## 9      male3 17688957
```

#### Per chromosome per sample


```r
cpg_per_chrom_per_sample <- foreach(i=1:length(cpg_files), .combine = rbind) %dopar%{
  cpg_files[[i]] %>% group_by(chr) %>% tally() %>% mutate(sample = names(cpg_files)[i])
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
foreach(i=1:length(cpg_files), .combine = rbind) %dopar%{
  cpg_files[[i]] %>% summarize(average_Cov = sum(cov)/n()) %>% mutate(sample = names(cpg_files)[i])
}
```

```
##    average_Cov     sample
## 1:    1.409210    female1
## 2:    1.414381    female2
## 3:    1.901558    female3
## 4:    2.013134 estradiol1
## 5:    2.103462 estradiol2
## 6:    2.696568 estradiol3
## 7:    3.191067      male1
## 8:    2.138167      male2
## 9:    3.048415      male3
```

### Coverage distribution


```r
cov_distribution <- foreach(i=1:length(cpg_files), .combine = rbind) %dopar%{
  cpg_files[[i]] %>% group_by(cov) %>% tally() %>% mutate(sample = names(cpg_files)[i])
}

cov_distribution %>%
  mutate(facet = gsub("[1-3]", "", sample)) %>%
  ggplot(aes(cov, n, group = sample, color = facet))+
  geom_line() +
  xlim(0,50)
```

```
## Warning: Removed 226 rows containing missing values (geom_path).
```

![](1-Methylation_sanity_check_files/figure-html/cov_distribution-1.png)

Distribution of coverage looks poisson at least despite low coverage.

### Average methylation per sample


```r
foreach(i=1:length(cpg_files), .combine = rbind) %dopar%{
  cpg_files[[i]] %>% summarize(mean_meth = mean(meth) %>% round(3)) %>% mutate(sample = names(cpg_files)[i])
}
```

```
##    mean_meth     sample
## 1:     0.763    female1
## 2:     0.777    female2
## 3:     0.776    female3
## 4:     0.772 estradiol1
## 5:     0.776 estradiol2
## 6:     0.777 estradiol3
## 7:     0.771      male1
## 8:     0.774      male2
## 9:     0.775      male3
```

Looks like methylation is generally very similar across samples

### Distribution of DNA methylation per sample


```r
dist_DNA_meth <- foreach(i=1:length(cpg_files), .combine = rbind) %dopar%{
  tmp <- hist(cpg_files[[i]]$meth, plot=FALSE, breaks = seq(0,1,0.05))
  data.table(
    breaks = tmp$breaks[-1],
    counts = tmp$counts 
  ) %>% mutate(sample = names(cpg_files)[i])
}

dist_DNA_meth %>%
  group_by(sample) %>%
  mutate(total = sum(counts),
         freq = counts/total) %>%
  ggplot(aes(x = breaks-0.025, y = freq, color = sample)) +
  geom_line()
```

![](1-Methylation_sanity_check_files/figure-html/dist_DNA_meth-1.png)

Yep, super low coverage = higher chance of 0% or 100%...

### Pairwise correlation of DNA methylation


```r
pairwise_cor <- cpg_files_merged %>%
  select(contains("meth")) %>%
  cor(use = "pairwise.complete.obs")

diag(pairwise_cor) <- NA

hclust <- hclust(as.dist(1-pairwise_cor))

pheatmap(pairwise_cor, cluster_rows = hclust, cluster_cols = hclust)
```

![](1-Methylation_sanity_check_files/figure-html/pairwise_cor-1.png)

Holy crap, the correlation values are atrocious. We might need to increase coverage by pooling replicates.

At least male and estradiol are clustering together.

