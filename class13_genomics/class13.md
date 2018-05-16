---
title: "Bioinformatics Class 13 Genomics"
author: "Amanda Alker"
date: "5/16/2018"
output: 
  html_document: 
    keep_md: yes
---



## 1000 Genomes Project data

Read 1000 genome data for MXL dataset 


```r
genotype <- read.csv("SampleGenotypes-Homo_sapiens_Variation_Sample_rs12936231.csv")
```

Looks Good!


```r
table(genotype[, 2]) / nrow(genotype) * 100
```

```
## 
##     C|C     C|G     G|C     G|G 
## 34.3750 32.8125 18.7500 14.0625
```

```r
#Access the second column
#Put these characterizations into a table 
#/row genotype creates the distribution
```

## Base quality scores from fastqsanger



```r
#install.packages("gtools")
```


```r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred 
```

```
##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```

## RNA-seq analysis 

Assessing genetic differences on a population scale. We want to find whether there is any association of the 4 asthma associated SNPs on ORMDL3 page 


```r
geno <- read.table("rs8067378_ENSG00000172057.6.txt")
```

What does this data look like?


```r
#numeric summary
summary(geno)
```

```
##      sample     geno          exp        
##  HG00096:  1   A/A:108   Min.   : 6.675  
##  HG00097:  1   A/G:233   1st Qu.:20.004  
##  HG00099:  1   G/G:121   Median :25.116  
##  HG00100:  1             Mean   :25.640  
##  HG00101:  1             3rd Qu.:30.779  
##  HG00102:  1             Max.   :51.518  
##  (Other):456
```

But is it different based on genotype? 


```r
summary(geno$exp[geno$geno == "A/A"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   11.40   27.02   31.25   31.82   35.92   51.52
```

```r
#Gives TRUE/FALSE of G/G
summary(geno$exp[geno$geno == "A/G"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   7.075  20.626  25.065  25.397  30.552  48.034
```

```r
summary(geno$exp[geno$geno == "G/G"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   6.675  16.903  20.074  20.594  24.457  33.956
```

ALWAYS PLOT THE DATA (boxplot)


```r
boxplot(exp ~ geno ,data = geno, notch = TRUE)
```

![](class13_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
#formula arguments needed where y~grp
#y is numeric vector of data value (expression)
#grp split into groups (genotypes )
#notch= TRUE do they overlap?
```


```r
library("ggplot2")
ggplot(geno, aes(geno, exp)) + geom_boxplot()
```

![](class13_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
#ggplot (df, aes) + geom_boxplot()
# + geom_jitter plots the points as well
#alpha = 0.2 creates transparency of boxes and 
```


```r
ggplot(geno, aes(exp, fill = geno)) + geom_density (alpha = 0.2)
```

![](class13_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

