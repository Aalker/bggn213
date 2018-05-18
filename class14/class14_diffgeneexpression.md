---
title: "Bioinformatics class 14 differential gene expression"
author: "Amanda Alker"
date: "5/18/2018"
output: 
  html_document: 
    keep_md: yes
---



Make sure we can download DESeq2

```r
#see the detailed installation handout on the website
library(BiocInstaller)
```

```
## Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
```

```
## A new version of Bioconductor is available after installing the most
##   recent version of R; see http://bioconductor.org/install
```

```r
#biocLite("DESeq2")
```


Read the data into R


```r
#stringsAsFactors makes catagorical factors but this program doesnt like them (ie treat as normal vector)
metadata <-read.csv("data/airway_metadata.csv", stringsAsFactors = FALSE)
counts <- read.csv("data/airway_scaledcounts.csv", stringsAsFactors = FALSE)
```


Let's look at them


```r
head(counts)
```

```
##           ensgene SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
## 1 ENSG00000000003        723        486        904        445       1170
## 2 ENSG00000000005          0          0          0          0          0
## 3 ENSG00000000419        467        523        616        371        582
## 4 ENSG00000000457        347        258        364        237        318
## 5 ENSG00000000460         96         81         73         66        118
## 6 ENSG00000000938          0          0          1          0          2
##   SRR1039517 SRR1039520 SRR1039521
## 1       1097        806        604
## 2          0          0          0
## 3        781        417        509
## 4        447        330        324
## 5         94        102         74
## 6          0          0          0
```


```r
#GEO collection is like the ENTREZ or SRA of NCBI
#check to make sure the IDs are distributed properly amongst row and columns 
head(metadata)
```

```
##           id     dex celltype     geo_id
## 1 SRR1039508 control   N61311 GSM1275862
## 2 SRR1039509 treated   N61311 GSM1275863
## 3 SRR1039512 control  N052611 GSM1275866
## 4 SRR1039513 treated  N052611 GSM1275867
## 5 SRR1039516 control  N080611 GSM1275870
## 6 SRR1039517 treated  N080611 GSM1275871
```



```r
colnames(counts)[-1] == metadata$id
```

```
## [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```

## Let's start to mess with differential gene expression (this isnt the real workflow)


```r
View(metadata)
```


```r
control <- metadata[metadata [, "dex"] =="control", ]
control
```

```
##           id     dex celltype     geo_id
## 1 SRR1039508 control   N61311 GSM1275862
## 3 SRR1039512 control  N052611 GSM1275866
## 5 SRR1039516 control  N080611 GSM1275870
## 7 SRR1039520 control  N061011 GSM1275874
```

```r
#check the logic to make sure that the metadata corresponds to the rows
#metadata[ with the inside function] will give the controls only

control.mean <- rowSums( counts[ ,control$id] )/4
#control$id will give us just the counts for the control ids

names(control.mean) <- counts$ensgene
#put the ids back on to make sure we dont get confused
```

How do we make this more robust?


```r
#get rid of the 4, because it is hardcoded
#nrow(control) should be substituted for 4
control.mean <- rowSums( counts[ ,control$id] )/nrow(control)
head(control.mean)
```

```
## [1] 900.75   0.00 520.50 339.75  97.25   0.75
```

Now lets do it again with the treated...


```r
treated <- metadata[metadata [, "dex"] =="treated", ]
control
```

```
##           id     dex celltype     geo_id
## 1 SRR1039508 control   N61311 GSM1275862
## 3 SRR1039512 control  N052611 GSM1275866
## 5 SRR1039516 control  N080611 GSM1275870
## 7 SRR1039520 control  N061011 GSM1275874
```

```r
treated.mean <- rowSums( counts[ ,treated$id] )/nrow(treated)

names(treated.mean) <- counts$ensgene
head(treated.mean)
```

```
## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 
##          658.00            0.00          546.00          316.50 
## ENSG00000000460 ENSG00000000938 
##           78.75            0.00
```


```r
meancounts <- data.frame(control.mean, treated.mean)
head(meancounts)
```

```
##                 control.mean treated.mean
## ENSG00000000003       900.75       658.00
## ENSG00000000005         0.00         0.00
## ENSG00000000419       520.50       546.00
## ENSG00000000457       339.75       316.50
## ENSG00000000460        97.25        78.75
## ENSG00000000938         0.75         0.00
```


```r
colSums(meancounts)
```

```
## control.mean treated.mean 
##     23005324     22196524
```


```r
plot(meancounts$control, meancounts$treated)
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
#notice all of the points are on top of each other. Normally, we fix this by plotting on a log scale
```


```r
plot(meancounts$control, meancounts$treated, log = "xy")
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0
## omitted from logarithmic plot
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0
## omitted from logarithmic plot
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
#make the log2foldchange in a new column
meancounts$log2fc <- log2(meancounts[,"treated.mean"]/meancounts[,"control.mean"])
head(meancounts)
```

```
##                 control.mean treated.mean      log2fc
## ENSG00000000003       900.75       658.00 -0.45303916
## ENSG00000000005         0.00         0.00         NaN
## ENSG00000000419       520.50       546.00  0.06900279
## ENSG00000000457       339.75       316.50 -0.10226805
## ENSG00000000460        97.25        78.75 -0.30441833
## ENSG00000000938         0.75         0.00        -Inf
```

We need to filter out the zeros because they are messing with the log plotting of our data

practice making a matrix to play with 

```r
x <- matrix(1:10, ncol = 2, byrow = TRUE)

x[5,2] <- 0
x
```

```
##      [,1] [,2]
## [1,]    1    2
## [2,]    3    4
## [3,]    5    6
## [4,]    7    8
## [5,]    9    0
```

```r
x == 0
```

```
##       [,1]  [,2]
## [1,] FALSE FALSE
## [2,] FALSE FALSE
## [3,] FALSE FALSE
## [4,] FALSE FALSE
## [5,] FALSE  TRUE
```

```r
which(x==0)
```

```
## [1] 10
```

```r
#want it to return row and column
which(x ==0, arr.ind = TRUE)
```

```
##      row col
## [1,]   5   2
```



```r
zero.vals <- which(meancounts[,1:2]== 0, arr.ind = TRUE)
head(zero.vals)
```

```
##                 row col
## ENSG00000000005   2   1
## ENSG00000004848  65   1
## ENSG00000004948  70   1
## ENSG00000005001  73   1
## ENSG00000006059 121   1
## ENSG00000006071 123   1
```

Now we need to remove the zero count containing genes 


```r
to.rm <- unique(zero.vals[,1])
mycounts <- meancounts[-to.rm,]
head(mycounts)
```

```
##                 control.mean treated.mean      log2fc
## ENSG00000000003       900.75       658.00 -0.45303916
## ENSG00000000419       520.50       546.00  0.06900279
## ENSG00000000457       339.75       316.50 -0.10226805
## ENSG00000000460        97.25        78.75 -0.30441833
## ENSG00000000971      5219.00      6687.50  0.35769358
## ENSG00000001036      2327.00      1785.75 -0.38194109
```

Differentially expressed rule of thumb is logfold change of 2


```r
up.ind <- mycounts$log2fc > 2
down.ind <- mycounts$log2fc < -2
head(up.ind)
```

```
## [1] FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
#this gives logical responses
#use sum to give how many are differentially expressed
```


```r
paste("Up:", sum(up.ind))
```

```
## [1] "Up: 250"
```

```r
paste("Down:", sum(down.ind))
```

```
## [1] "Down: 367"
```

```r
#use the paste function to have the names of your choice printed out
```


```r
anno <- read.csv("data/annotables_grch38.csv")
head(anno)
```

```
##           ensgene entrez   symbol chr     start       end strand
## 1 ENSG00000000003   7105   TSPAN6   X 100627109 100639991     -1
## 2 ENSG00000000005  64102     TNMD   X 100584802 100599885      1
## 3 ENSG00000000419   8813     DPM1  20  50934867  50958555     -1
## 4 ENSG00000000457  57147    SCYL3   1 169849631 169894267     -1
## 5 ENSG00000000460  55732 C1orf112   1 169662007 169854080      1
## 6 ENSG00000000938   2268      FGR   1  27612064  27635277     -1
##          biotype
## 1 protein_coding
## 2 protein_coding
## 3 protein_coding
## 4 protein_coding
## 5 protein_coding
## 6 protein_coding
##                                                                                                  description
## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
## 2                                                            tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
## 3 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
## 4                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
## 5                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
## 6                          FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]
```

Use the merge function

```r
results <- merge(mycounts, anno, by.x = "row.names", by.y = "ensgene")
head(results)
```

```
##         Row.names control.mean treated.mean      log2fc entrez   symbol
## 1 ENSG00000000003       900.75       658.00 -0.45303916   7105   TSPAN6
## 2 ENSG00000000419       520.50       546.00  0.06900279   8813     DPM1
## 3 ENSG00000000457       339.75       316.50 -0.10226805  57147    SCYL3
## 4 ENSG00000000460        97.25        78.75 -0.30441833  55732 C1orf112
## 5 ENSG00000000971      5219.00      6687.50  0.35769358   3075      CFH
## 6 ENSG00000001036      2327.00      1785.75 -0.38194109   2519    FUCA2
##   chr     start       end strand        biotype
## 1   X 100627109 100639991     -1 protein_coding
## 2  20  50934867  50958555     -1 protein_coding
## 3   1 169849631 169894267     -1 protein_coding
## 4   1 169662007 169854080      1 protein_coding
## 5   1 196651878 196747504      1 protein_coding
## 6   6 143494811 143511690     -1 protein_coding
##                                                                                                  description
## 1                                                          tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
## 2 dolichyl-phosphate mannosyltransferase polypeptide 1, catalytic subunit [Source:HGNC Symbol;Acc:HGNC:3005]
## 3                                               SCY1-like, kinase-like 3 [Source:HGNC Symbol;Acc:HGNC:19285]
## 4                                    chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
## 5                                                     complement factor H [Source:HGNC Symbol;Acc:HGNC:4883]
## 6                                          fucosidase, alpha-L- 2, plasma [Source:HGNC Symbol;Acc:HGNC:4008]
```

```r
#works with a table of annotations that you want to go with
```


```r
library("AnnotationDbi")
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```r
#biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")
```

```
## 
```


```r
columns(org.Hs.eg.db)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
## [25] "UNIGENE"      "UNIPROT"
```


```r
mycounts$symbol <- mapIds(org.Hs.eg.db,
                keys=row.names(mycounts),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
#mapIds add individual columns to our results table 
head(mycounts)
```

```
##                 control.mean treated.mean      log2fc   symbol
## ENSG00000000003       900.75       658.00 -0.45303916   TSPAN6
## ENSG00000000419       520.50       546.00  0.06900279     DPM1
## ENSG00000000457       339.75       316.50 -0.10226805    SCYL3
## ENSG00000000460        97.25        78.75 -0.30441833 C1orf112
## ENSG00000000971      5219.00      6687.50  0.35769358      CFH
## ENSG00000001036      2327.00      1785.75 -0.38194109    FUCA2
```




```r
library("DESeq2")
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:base':
## 
##     apply
```

```r
citation("DESeq2")
```

```
## 
##   Love, M.I., Huber, W., Anders, S. Moderated estimation of fold
##   change and dispersion for RNA-seq data with DESeq2 Genome
##   Biology 15(12):550 (2014)
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
##     author = {Michael I. Love and Wolfgang Huber and Simon Anders},
##     year = {2014},
##     journal = {Genome Biology},
##     doi = {10.1186/s13059-014-0550-8},
##     volume = {15},
##     issue = {12},
##     pages = {550},
##   }
```


```r
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
```

```
## converting counts to integer mode
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
## design formula are characters, converting to factors
```

```r
dds
```

```
## class: DESeqDataSet 
## dim: 38694 8 
## metadata(1): version
## assays(1): counts
## rownames(38694): ENSG00000000003 ENSG00000000005 ...
##   ENSG00000283120 ENSG00000283123
## rowData names(0):
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(4): id dex celltype geo_id
```

To get the results


```r
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```


```r
res <- results(dds)
res
```

```
## log2 fold change (MLE): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 38694 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE       stat     pvalue
##                 <numeric>      <numeric> <numeric>  <numeric>  <numeric>
## ENSG00000000003 747.19420    -0.35070283 0.1682342 -2.0846111 0.03710462
## ENSG00000000005   0.00000             NA        NA         NA         NA
## ENSG00000000419 520.13416     0.20610652 0.1010134  2.0403876 0.04131173
## ENSG00000000457 322.66484     0.02452714 0.1451103  0.1690242 0.86577762
## ENSG00000000460  87.68263    -0.14714409 0.2569657 -0.5726216 0.56690095
## ...                   ...            ...       ...        ...        ...
## ENSG00000283115  0.000000             NA        NA         NA         NA
## ENSG00000283116  0.000000             NA        NA         NA         NA
## ENSG00000283119  0.000000             NA        NA         NA         NA
## ENSG00000283120  0.974916     -0.6682308  1.694063 -0.3944544  0.6932456
## ENSG00000283123  0.000000             NA        NA         NA         NA
##                      padj
##                 <numeric>
## ENSG00000000003 0.1630257
## ENSG00000000005        NA
## ENSG00000000419 0.1757326
## ENSG00000000457 0.9616577
## ENSG00000000460 0.8157061
## ...                   ...
## ENSG00000283115        NA
## ENSG00000283116        NA
## ENSG00000283119        NA
## ENSG00000283120        NA
## ENSG00000283123        NA
```


```r
summary(res)
```

```
## 
## out of 25258 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 1564, 6.2% 
## LFC < 0 (down)   : 1188, 4.7% 
## outliers [1]     : 142, 0.56% 
## low counts [2]   : 9971, 39% 
## (mean count < 10)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Order the results by p-value


```r
resOrdered <- res[order(res$pvalue), ]
head(resOrdered)
```

```
## log2 fold change (MLE): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 6 rows and 6 columns
##                   baseMean log2FoldChange      lfcSE      stat
##                  <numeric>      <numeric>  <numeric> <numeric>
## ENSG00000152583   954.7709       4.368359 0.23713648  18.42129
## ENSG00000179094   743.2527       2.863888 0.17555825  16.31304
## ENSG00000116584  2277.9135      -1.034700 0.06505273 -15.90556
## ENSG00000189221  2383.7537       3.341544 0.21241508  15.73120
## ENSG00000120129  3440.7038       2.965211 0.20370277  14.55656
## ENSG00000148175 13493.9204       1.427168 0.10036663  14.21955
##                       pvalue         padj
##                    <numeric>    <numeric>
## ENSG00000152583 8.867079e-76 1.342919e-71
## ENSG00000179094 7.972621e-60 6.037267e-56
## ENSG00000116584 5.798513e-57 2.927283e-53
## ENSG00000189221 9.244206e-56 3.500088e-52
## ENSG00000120129 5.306416e-48 1.607313e-44
## ENSG00000148175 6.929711e-46 1.749175e-42
```

```r
#need to rearrange the whole dataframe
```


```r
res05 <- results(dds, alpha=0.05)
summary(res05)
```

```
## 
## out of 25258 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 1237, 4.9% 
## LFC < 0 (down)   : 933, 3.7% 
## outliers [1]     : 142, 0.56% 
## low counts [2]   : 9033, 36% 
## (mean count < 6)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
#change the default cutoff value from 0.1 to 0.5 (p-value)
```


```r
resSig05 <- subset(as.data.frame(res), padj < 0.05)
nrow(resSig05)
```

```
## [1] 2182
```

```r
#use the subset() function to access a smaller part of the dataset
```


```r
resSig01 <- subset(as.data.frame(res), padj < 0.01)
nrow(resSig01)
```

```
## [1] 1437
```

Add annotation to res01 results dataframe


```r
resSig01$symbol <- mapIds(org.Hs.eg.db,
                keys=row.names(resSig01),
                column = "SYMBOL", 
                keytype= "ENSEMBL", 
                multiVals = "first")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

Arrange by adjusted p-value


```r
ord <- order( resSig01$padj )
#View(res01[ord,])
head(resSig01[ord,])
```

```
##                   baseMean log2FoldChange      lfcSE      stat
## ENSG00000152583   954.7709       4.368359 0.23713648  18.42129
## ENSG00000179094   743.2527       2.863888 0.17555825  16.31304
## ENSG00000116584  2277.9135      -1.034700 0.06505273 -15.90556
## ENSG00000189221  2383.7537       3.341544 0.21241508  15.73120
## ENSG00000120129  3440.7038       2.965211 0.20370277  14.55656
## ENSG00000148175 13493.9204       1.427168 0.10036663  14.21955
##                       pvalue         padj  symbol
## ENSG00000152583 8.867079e-76 1.342919e-71 SPARCL1
## ENSG00000179094 7.972621e-60 6.037267e-56    PER1
## ENSG00000116584 5.798513e-57 2.927283e-53 ARHGEF2
## ENSG00000189221 9.244206e-56 3.500088e-52    MAOA
## ENSG00000120129 5.306416e-48 1.607313e-44   DUSP1
## ENSG00000148175 6.929711e-46 1.749175e-42    STOM
```

Write this to a .csv file


```r
write.csv(resSig01[ord,], "signif01_results.csv")
```

##Data Visualization

Get the gene ID for the "CRISPLD2" gene

```r
index <- grep("CRISPLD2", resSig01$symbol)
resSig01 [index, ]
```

```
##                 baseMean log2FoldChange     lfcSE     stat       pvalue
## ENSG00000103196 3096.159       2.626034 0.2674705 9.818031 9.416441e-23
##                         padj   symbol
## ENSG00000103196 3.395524e-20 CRISPLD2
```


```r
rownames(resSig01[index, ])
```

```
## [1] "ENSG00000103196"
```

Now lets plot the counts


```r
plotCounts(dds, gene="ENSG00000103196", intgroup="dex")
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-39-1.png)<!-- -->


```r
#Return the data instead of the counts
d <- plotCounts(dds, gene = "ENSG00000103196", intgroup = "dex", returnData = TRUE)
head(d)
```

```
##                count     dex
## SRR1039508  774.5002 control
## SRR1039509 6258.7915 treated
## SRR1039512 1100.2741 control
## SRR1039513 6093.0324 treated
## SRR1039516  736.9483 control
## SRR1039517 2742.1908 treated
```

Now from this object returned, we can plot a boxplot

```r
boxplot(count ~ dex, data = d)
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-41-1.png)<!-- -->


```r
library(ggplot2)
ggplot(d, aes(dex, count)) + geom_boxplot(aes(fill=dex)) + scale_y_log10() + ggtitle("CRISPLD2")
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

## MA and Volcano plots 

Add a column called sig to our rull "res" that calls true if padj <0.05 and FALSE if not

```r
res$sig <- res$padj<0.05
table(res$sig)
```

```
## 
## FALSE  TRUE 
## 12963  2182
```


```r
sum(is.na(res$sig))
```

```
## [1] 23549
```

MA plot shows average expression on the x-axis and log fold change on the y-axis 
Volcano plot shows the log fold change on the X-axis and the -log10 of the p-value on the Y-axis (the more significant the p-value, the larger the -log10 value will be)

In built MA- plot

```r
plotMA(res, ylim = c(-2,2))
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-45-1.png)<!-- -->

Remove the noise associated with the log2 fold changes from low gene counts


```r
resLFC <- lfcShrink(dds, coef = 2)
resLFC
```

```
## log2 fold change (MAP): dex treated vs control 
## Wald test p-value: dex treated vs control 
## DataFrame with 38694 rows and 6 columns
##                  baseMean log2FoldChange      lfcSE       stat     pvalue
##                 <numeric>      <numeric>  <numeric>  <numeric>  <numeric>
## ENSG00000000003 747.19420    -0.31838595 0.15271739 -2.0846111 0.03710462
## ENSG00000000005   0.00000             NA         NA         NA         NA
## ENSG00000000419 520.13416     0.19883048 0.09744556  2.0403876 0.04131173
## ENSG00000000457 322.66484     0.02280238 0.13491699  0.1690242 0.86577762
## ENSG00000000460  87.68263    -0.11887370 0.20772938 -0.5726216 0.56690095
## ...                   ...            ...        ...        ...        ...
## ENSG00000283115  0.000000             NA         NA         NA         NA
## ENSG00000283116  0.000000             NA         NA         NA         NA
## ENSG00000283119  0.000000             NA         NA         NA         NA
## ENSG00000283120  0.974916    -0.05944174  0.1514839 -0.3944544  0.6932456
## ENSG00000283123  0.000000             NA         NA         NA         NA
##                      padj
##                 <numeric>
## ENSG00000000003 0.1630257
## ENSG00000000005        NA
## ENSG00000000419 0.1757326
## ENSG00000000457 0.9616577
## ENSG00000000460 0.8157061
## ...                   ...
## ENSG00000283115        NA
## ENSG00000283116        NA
## ENSG00000283119        NA
## ENSG00000283120        NA
## ENSG00000283123        NA
```

```r
plotMA(resLFC, ylim = c(-2,2))
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-46-1.png)<!-- -->

Make a volcano plot

```r
ggplot(as.data.frame(res), aes(log2FoldChange, -1*log10(pvalue), col= sig )) + geom_point() + ggtitle("Volcano plot")
```

```
## Warning: Removed 13578 rows containing missing values (geom_point).
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

Log Transformation of results is needed for heatmaps, PCA, or clustering


```r
#variance stabilizing transformation (VST)
vsdata <- vst(dds, blind=FALSE)
```

##PCA 


```r
plotPCA(vsdata, intgroup = "dex")
```

![](class14_diffgeneexpression_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

Now show the session information


```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.3
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggplot2_2.2.1              DESeq2_1.18.1             
##  [3] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
##  [5] matrixStats_0.53.1         GenomicRanges_1.30.3      
##  [7] GenomeInfoDb_1.14.0        org.Hs.eg.db_3.5.0        
##  [9] AnnotationDbi_1.40.0       IRanges_2.12.0            
## [11] S4Vectors_0.16.0           Biobase_2.38.0            
## [13] BiocGenerics_0.24.0        BiocInstaller_1.28.0      
## 
## loaded via a namespace (and not attached):
##  [1] locfit_1.5-9.1         Rcpp_0.12.16           lattice_0.20-35       
##  [4] rprojroot_1.3-2        digest_0.6.15          plyr_1.8.4            
##  [7] backports_1.1.2        acepack_1.4.1          RSQLite_2.1.0         
## [10] evaluate_0.10.1        pillar_1.2.2           zlibbioc_1.24.0       
## [13] rlang_0.2.0            lazyeval_0.2.1         annotate_1.56.2       
## [16] rstudioapi_0.7         data.table_1.10.4-3    blob_1.1.1            
## [19] rpart_4.1-13           Matrix_1.2-14          checkmate_1.8.5       
## [22] rmarkdown_1.9          labeling_0.3           splines_3.4.4         
## [25] BiocParallel_1.12.0    geneplotter_1.56.0     stringr_1.3.0         
## [28] foreign_0.8-69         htmlwidgets_1.2        RCurl_1.95-4.10       
## [31] bit_1.1-12             munsell_0.4.3          compiler_3.4.4        
## [34] pkgconfig_2.0.1        base64enc_0.1-3        htmltools_0.3.6       
## [37] nnet_7.3-12            tibble_1.4.2           gridExtra_2.3         
## [40] htmlTable_1.11.2       GenomeInfoDbData_1.0.0 Hmisc_4.1-1           
## [43] XML_3.98-1.11          bitops_1.0-6           grid_3.4.4            
## [46] xtable_1.8-2           gtable_0.2.0           DBI_0.8               
## [49] magrittr_1.5           scales_0.5.0           stringi_1.1.7         
## [52] XVector_0.18.0         genefilter_1.60.0      latticeExtra_0.6-28   
## [55] Formula_1.2-3          RColorBrewer_1.1-2     tools_3.4.4           
## [58] bit64_0.9-7            survival_2.42-3        yaml_2.1.18           
## [61] colorspace_1.3-2       cluster_2.0.7-1        memoise_1.1.0         
## [64] knitr_1.20
```

