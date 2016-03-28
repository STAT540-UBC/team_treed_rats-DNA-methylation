# 2-Calling_DMRs.Rmd
Tony  
February 18, 2016  


```r
require(data.table)
require(foreach)
require(doMC)
require(bsseq)
require(ggplot2)
require(tidyr)
require(knitr)
require(dplyr)
require(GenomicRanges)
```

# Merge data


```r
load("../methylation_data/cpg_files_merged.RData")

bssmooth_cpgs <- lapply(cpg_files_merged, function(x) {
  tmp <- x %>% select(chr, pos) %>% filter(pos >0)
}) %>% rbindlist() %>% unique

registerDoMC(length(cpg_files_merged))
bssmooth <- foreach(f=names(cpg_files_merged), .combine = cbind) %dopar% {
  tmp <- cpg_files_merged[[f]] %>% setnames(c("chr", "pos", paste0(f, "_meth"), paste0(f, "_cov"))) %>% data.table()
  left_join(bssmooth_cpgs, tmp, by = c("chr","pos")) %>% select(-chr,-pos)
} %>% cbind(bssmooth_cpgs, .)

bssmooth[is.na(bssmooth)] <- 0

samples <- names(cpg_files_merged)

rm(bssmooth_cpgs)

bssmooth_smooth <- BSseq(chr = bssmooth$chr, pos = bssmooth$pos, 
                         M = bssmooth %>% select(contains("meth")) %>% as.matrix(), 
                         Cov = bssmooth %>% select(contains("cov")) %>% as.matrix(),
                         sampleNames = samples)

chroms <- granges(bssmooth_smooth) %>% seqnames %>% levels

# rm(bssmooth)

bssmooth_smooth <- BSmooth(bssmooth_smooth, verbose = TRUE, parallelBy = "chromosome", mc.cores = length(chroms))

pData(bssmooth_smooth)$col <- c("#7fc97f","#beaed4","#fdc086")

# save(... = bssmooth_smooth, file = "../methylation_data/bssmooth_smooth.RData", compress = T)
```


```r
samples <- pData(bssmooth_smooth) %>% rownames %>% gsub("_meth_cov", "", x = .)

bssmooth_dt <- cbind(
  as.data.frame(granges(bssmooth_smooth)) %>% setnames(colnames(.[1]), "chr"),
  getMeth(bssmooth_smooth) %>% as.data.frame() %>% setnames(paste0(samples, "_meth")),
  getCoverage(bssmooth_smooth) %>% as.data.frame() %>% setnames(paste0(samples, "_cov"))
) %>% data.table() %>% select(-end, -width, -strand)

save(bssmooth_smooth, bssmooth_dt, file = "../methylation_data/bssmooth_smooth.RData", compress = F)
# load("bssmooth_smooth.RData")
```

# Call DMRs for female vs male


```r
load("/projects/epigenomics/users/thui/stat540/methylation_data/bssmooth_smooth.RData")

rn6_genes <- fread("homer/rn6_tss_raw.gtf", skip = 1) %>%
  select(nearest.promoter = V1, name = V7) %>% as.data.frame() %>% tbl_df() 
```


```r
bssmooth_dt_maleVsFemale <- bssmooth_dt %>%
  tbl_df() %>%
  select(-contains("estradiol")) %>%
  # filter(female_cov > 0, male_cov > 0) %>%
  mutate(diff = female_meth - male_meth)
```


```r
top_percent <- 1

# summary(bssmooth_dt_maleVsFemale$diff)
diff_quantile <- quantile(x = bssmooth_dt_maleVsFemale$diff, probs = c(top_percent/100/2, 1-(top_percent/100/2)))
diff_quantile
```

```
##       0.5%      99.5% 
## -0.1876567  0.1935138
```

```r
# ggplot(bssmooth_dt_maleVsFemale, aes(diff)) +
#   geom_density() +
#   geom_vline(xintercept = diff_quantile)
```


```r
binsize <- 500
min_cpg <- 3

bssmooth_dt_maleVsFemale_dCPG <- bssmooth_dt_maleVsFemale %>%
  filter(diff < diff_quantile[1] | diff > diff_quantile[2]) %>%
  group_by(chr) %>%
  mutate(dist = c(binsize+1,diff(start)),
         diff_cumul = diff*lag(diff, n = 1)) %>%
  ungroup() %>% mutate(chr = as.character(chr)) %>%
  mutate(bin = 1+cumsum(dist > binsize | diff_cumul < 0 | is.na(diff_cumul))) %>%
  setnames(colnames(.)[2], "pos")

bssmooth_dt_maleVsFemale_DMR <- bssmooth_dt_maleVsFemale_dCPG %>%
  filter(female_cov > 0) %>%
  group_by(chr, bin) %>%
  summarize(start=min(pos),
            end=max(pos), 
            mean_female = mean(female_meth) %>% round(3), 
            mean_male = mean(male_meth) %>% round(3),
            female_cov = sum(female_cov),
            male_cov = sum(male_cov),
            # remove DMRs within transition zones
            slope_clusters = min(
              abs(female_meth[1] - female_meth[n()]),
              abs(male_meth[1] - male_meth[n()])
              ),
            num_cpg = n()
            ) %>%
  filter(num_cpg >= min_cpg) %>%
  filter(slope_clusters < 0.2) %>%
  # filter(female_cov >= 3) %>%
  select(-slope_clusters, -num_cpg)
```

## Annotate DMRs


```r
# save(bssmooth_dt_maleVsFemale_DMR, file = "methylation_data/maleVSfemaleDMRs.RData")
write.csv(x = bssmooth_dt_maleVsFemale_DMR, file = "maleVSfemaleDMRs.csv", quote = F, row.names = F)

write.table(x = bssmooth_dt_maleVsFemale_DMR %>% select(bin,chr,start,end) %>% mutate(strand = "+"), file = "homer/maleVSfemaleDMRs.bed", quote = F, row.names = F, col.names = F, sep = "\t")

system("cd homer; sh get_DE_tss.sh")
system("cd homer; sh homer_annotate.sh maleVSfemaleDMRs.bed")
```


```r
bssmooth_dt_maleVsFemale_DMR_annotated <- read.table("homer/maleVSfemaleDMRs.bed.annotated", header = T, sep = "\t", quote="")

bssmooth_dt_maleVsFemale_DMR_annotated <- bssmooth_dt_maleVsFemale_DMR_annotated %>%
  select(1:11) %>%
  setnames(c("bin", "chr", "start", "end", "strand", "score", "ratio", "annotation", "detailed.annotation", "dist.to.tss", "nearest.promoter")) %>%
  left_join(., bssmooth_dt_maleVsFemale_DMR) %>%
  tbl_df()
```

```
## Joining by: c("bin", "chr", "start", "end")
```

```
## Warning in left_join_impl(x, y, by$x, by$y): joining character vector and
## factor, coercing into character vector
```

```r
maleVsFemale_DMR_genes <- inner_join(bssmooth_dt_maleVsFemale_DMR_annotated, rn6_genes) %>%
  select(chr, start, end, mean_female, mean_male, annotation, dist.to.tss, nearest.promoter, name, bin)
```

```
## Joining by: "nearest.promoter"
```

```
## Warning in inner_join_impl(x, y, by$x, by$y): joining factor and character
## vector, coercing into character vector
```

### Remove irrelevant DMRs

* Remove DMRs that are promoters of other genes
* Remove DMRs < 50000 away from nearest promoter


```r
maleVsFemale_DMR_genes_filtered <- maleVsFemale_DMR_genes %>%
  mutate(annotation = gsub("\\(.*", "", annotation) %>% gsub(" ", "", .)) %>%
  filter(abs(dist.to.tss) < 50000) %>%
  # filter(slope_clusters < 0.1) %>%
  filter(ifelse(grepl("promoter-TSS", annotation), abs(dist.to.tss) < 3000, TRUE))

ggplot(maleVsFemale_DMR_genes_filtered, aes("DMR location", fill = annotation)) +
  geom_bar(position = position_fill()) +
  ylab("Fraction of all DMRs") +
  xlab("")
```

![](2-Calling_DMRs_files/figure-html/unnamed-chunk-5-1.png)

## Finalize file


```r
de_genes <- read.table("../Data_Analysis/RNAseq_result/DE_genes/glmQLFit_DE_genes.tsv", header = TRUE)

dmr_set <- maleVsFemale_DMR_genes_filtered %>%
  mutate(hypo_in_female = mean_female < mean_male) %>%
  select(chr, start, end, annotation, dist.to.tss, gene = nearest.promoter, name, hypo_in_female) %>%
  left_join(., de_genes) %>%
  mutate(epi_regulation = hypo_in_female == gExp_up_in_female) %>%
  select(-V7) %>%
  # filter(!annotation %in% c("exon", "intron", "intergenic")) %>%
  arrange(name, dist.to.tss)
```

```
## Joining by: "gene"
```

```
## Warning in left_join_impl(x, y, by$x, by$y): joining factor and character
## vector, coercing into character vector
```

```r
table(dmr_set$epi_regulation)
```

```
## 
## FALSE  TRUE 
##   103    79
```

```r
epi_genes <- dmr_set %>%
  group_by(gene, epi_regulation) %>%
  tally() %>%
  spread(key = epi_regulation, value = n, fill = 0) %>%
  setnames(c("gene", "non_epi_regulated", "epi_regulated")) %>%
  mutate(frac_regulated = epi_regulated/(epi_regulated+non_epi_regulated))

epi_gene_list <- epi_genes %>%
  filter(frac_regulated >= 0.5) %>%
  .$gene

qplot(epi_genes$frac_regulated)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](2-Calling_DMRs_files/figure-html/unnamed-chunk-6-1.png)

```r
final_table <- dmr_set %>% 
  filter(gene %in% epi_gene_list) %>%
  filter(epi_regulation)
  
write.table(final_table, file = "2_genes_DMR_associations.tsv", row.names = F, quote = F, sep = "\t")
```

# Now, the same but for female vs female+ estradiol


```r
bssmooth_dt_FemaleVsEstradiol <- bssmooth_dt %>%
  tbl_df() %>%
  select(-starts_with("male")) %>%
  mutate(diff = female_meth - estradiol_meth)
```


```r
top_percent <- 1

# summary(bssmooth_dt_maleVsFemale$diff)
diff_quantile <- quantile(x = bssmooth_dt_FemaleVsEstradiol$diff, probs = c(top_percent/100/2, 1-(top_percent/100/2)))
diff_quantile
```

```
##       0.5%      99.5% 
## -0.1992926  0.1906262
```

```r
# ggplot(bssmooth_dt_maleVsFemale, aes(diff)) +
#   geom_density() +
#   geom_vline(xintercept = diff_quantile)
```


```r
bssmooth_dt_FemaleVsEstradiol_dCPG <- bssmooth_dt_FemaleVsEstradiol %>%
  filter(diff < diff_quantile[1] | diff > diff_quantile[2]) %>%
  group_by(chr) %>%
  mutate(dist = c(binsize+1,diff(start)),
         diff_cumul = diff*lag(diff, n = 1)) %>%
  ungroup() %>% mutate(chr = as.character(chr)) %>%
  mutate(bin = 1+cumsum(dist > binsize | diff_cumul < 0 | is.na(diff_cumul))) %>%
  setnames(colnames(.)[2], "pos")

bssmooth_dt_FemaleVsEstradiol_DMR <- bssmooth_dt_FemaleVsEstradiol_dCPG %>%
  filter(female_cov > 0) %>%
  # filter(bin == 2) %>%
  group_by(chr, bin) %>%
  summarise(start = min(pos), end = max(pos)+2, 
            mean_female = mean(female_meth) %>% round(3), 
            mean_estradiol = mean(estradiol_meth) %>% round(3),
            female_cov = sum(female_cov),
            estradiol_cov = sum(estradiol_cov),
            # remove DMRs within transition zones
            slope_clusters = min(
              abs(female_meth[1] - female_meth[n()]),
              abs(estradiol_meth[1] - estradiol_meth[n()])
              ),
            num_cpg = n()
            ) %>%
  filter(num_cpg >= min_cpg) %>%
  select(-slope_clusters, -num_cpg)
```

## Subset DMRs that are overlapping the epigenetically regulated male + female DMRs


```r
#maleVsFemale_DMR_genes_filtered
#final_table

bssmooth_dt_maleVsFemale_DMR_GRanges <- final_table %>%
  arrange(chr, start) %>%
  makeGRangesFromDataFrame(keep.extra.columns = F)

bssmooth_dt_FemaleVsEstradiol_DMR_GRanges <- bssmooth_dt_FemaleVsEstradiol_DMR %>%
  makeGRangesFromDataFrame(keep.extra.columns = F)

overlaps <- findOverlaps(bssmooth_dt_maleVsFemale_DMR_GRanges, bssmooth_dt_FemaleVsEstradiol_DMR_GRanges, select = "all")

bssmooth_dt_FemaleVsEstradiol_DMR_overlap <- bssmooth_dt_FemaleVsEstradiol_DMR[subjectHits(overlaps),] %>% arrange(chr, start) %>% filter(chr != "chrY") %>%
  mutate(hypo_in_female = mean_female < mean_estradiol)

bsmooth_dt_all_DMR_overlap <- (final_table %>% arrange(chr, start))[queryHits(overlaps),] %>% filter(chr != "chrY")

bsmooth_dt_all_DMR_overlap$mean_estradiol <- bssmooth_dt_FemaleVsEstradiol_DMR_overlap$hypo_in_female
```

Called concordant if methylation is in same direction


```r
bsmooth_dt_all_DMR_overlap_final <- bsmooth_dt_all_DMR_overlap %>%
  select(chr, start, end, annotation, dist.to.tss, gene, name, hypo_in_female, gExp_up_in_female) %>%
  unique()

# bsmooth_dt_all_DMR_overlap_final %>% head %>% kable("markdown")
```

## Subset genes that are DE between female and both male and female+zeb


```r
final_DE_genes <- read.table(file = "../Data_Analysis/RNAseq_result/DE_genes/femVSfemZeb_glmQLFit_DE_genes.tsv", header=TRUE) %>%
  select(gene, name = V7)

rn6_de_genes_track <- fread("homer/rn6_tss_raw.gtf", skip = 1) %>%
  filter(V1 %in% final_DE_genes$gene) %>%
  select(chr = V2, start = V4, end = V5, gene = V1)

rn6_genes_track <- fread("homer/rn6_tss_raw.gtf", skip = 1) %>%
  select(chr = V2, start = V4, end = V5, gene = V1)
```


```r
final_DMR_set <- bsmooth_dt_all_DMR_overlap_final %>%
  filter(abs(dist.to.tss) < 20000) %>%
  inner_join(., final_DE_genes)
```

```
## Joining by: c("gene", "name")
```

```
## Warning in inner_join_impl(x, y, by$x, by$y): joining character vector and
## factor, coercing into character vector

## Warning in inner_join_impl(x, y, by$x, by$y): joining character vector and
## factor, coercing into character vector
```

```r
final_DMR_set %>% kable("markdown")
```



|chr   |     start|       end|annotation | dist.to.tss|gene                 |name   |hypo_in_female |gExp_up_in_female |
|:-----|---------:|---------:|:----------|-----------:|:--------------------|:------|:--------------|:-----------------|
|chr1  |  15779390|  15779992|Intergenic |       -2784|ENSRNOT00000088025.1 |Bclaf1 |TRUE           |TRUE              |
|chr10 |  85051230|  85051912|Intergenic |       -2240|ENSRNOT00000012538.5 |Tbx21  |TRUE           |TRUE              |
|chr10 | 109502106| 109502623|Intergenic |       18482|ENSRNOT00000054976.4 |Actg1  |TRUE           |TRUE              |
|chr15 |  61677824|  61678114|intron     |      -14355|ENSRNOT00000015517.5 |Kbtbd6 |TRUE           |TRUE              |
|chr15 |  61682395|  61683201|intron     |       -9526|ENSRNOT00000015517.5 |Kbtbd6 |TRUE           |TRUE              |
|chr3  | 170569405| 170570106|Intergenic |       19443|ENSRNOT00000006991.5 |Tfap2c |TRUE           |TRUE              |

```r
final_DMR_set_visualize <- final_DMR_set %>%
  left_join(., rn6_de_genes_track, by = "gene") %>%
  mutate(start = pmin(start.x, start.y),
         end = pmax(end.x, end.y))

regions_final_DMR_set <- makeGRangesFromDataFrame(final_DMR_set_visualize, seqnames.field = "chr.x", keep.extra.columns = T)
```


```r
pdf("final_DMR_plots.pdf")
for (i in seq_along(regions_final_DMR_set)) {
  title <- paste("Transcript", regions_final_DMR_set$gene[i], "for gene", regions_final_DMR_set$name[i])
  plotRegion(BSseq = bssmooth_smooth, region = regions_final_DMR_set[i], extend = 10000, addRegions = makeGRangesFromDataFrame(final_DMR_set), main = title,
             annoTrack = list(
               de = makeGRangesFromDataFrame(rn6_de_genes_track),
               genes= makeGRangesFromDataFrame(rn6_genes_track)
             ))
}
dev.off()
```

```
## png 
##   2
```





