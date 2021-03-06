---
title: "Bioinformatics Class7"
author: "Amanda Alker"
date: "4/25/2018"
output: 
  html_document: 
    keep_md: yes
---



# Functions again

We can source any file of R code with the 'source()' function

```r
source("http://tinyurl.com/rescale-R")
```

Lets make sure things are here


```r
ls()
```

```
##  [1] "both_na"         "both_na2"        "both_na3"       
##  [4] "df1"             "df2"             "df3"            
##  [7] "gene_intersect"  "gene_intersect2" "gene_intersect3"
## [10] "gene_intersect4" "rescale"         "rescale2"
```

check our 'rescale()' function is working


```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```


```r
rescale( c(1:10, "string") )
```

Let's check if 'rescale2()' does any better

```r
rescale2( c(1:10, "string") )
```

## Function for finding missing values in two datasets

Write a 'both_na() function to do this

```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

provides where the NAs are located 

```r
which(is.na(x))
```

```
## [1] 3 5
```

provides how many NAs there are 

```r
sum( is.na(x))
```

```
## [1] 2
```

This will print these vectors seperately

```r
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
is.na(y)
```

```
## [1]  TRUE FALSE  TRUE FALSE FALSE
```

we can use & to combine the vectors

```r
sum(is.na(x) & is.na(y))
```

```
## [1] 1
```

now we can make our function

```r
both_na <-function(x, y) {
  sum( is.na(x) & is.na(y))
}
```

test it

```r
both_na(x, y)
```

```
## [1] 1
```


```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both_na(x, y1)
```

```
## [1] 2
```

```r
both_na(x, y2)
```

```
## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
## shorter object length
```

```
## [1] 3
```


```r
#both_na2(x, y2)
```


```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

ans <-both_na3(x, y)
```

```
## Found 1 NA's at position(s):3
```


```r
ans$which
```

```
## [1] 3
```

##One more example- actually useful

```r
x <- df1$IDs
y <- df2$IDs

x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```

We can try the 'intersect()' function and the '%in%' functions


```r
intersect(x, y)
```

```
## [1] "gene2" "gene3"
```

```r
x %in% y
```

```
## [1] FALSE  TRUE  TRUE
```

```r
y %in% x
```

```
## [1]  TRUE FALSE  TRUE FALSE
```

we can use the logical output in '%in% to get at our matching data..


```r
x[x %in% y]
```

```
## [1] "gene2" "gene3"
```

```r
y[y %in% x]
```

```
## [1] "gene2" "gene3"
```


Lets put these together as columns of a matrix


```r
cbind(x[x %in% y], y[y %in% x])
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```


```r
cbind( c("hello", "Help"), c("Please", "Me"))
```

```
##      [,1]    [,2]    
## [1,] "hello" "Please"
## [2,] "Help"  "Me"
```

```r
rbind( c("hello", "Help"), c("Please", "Me"))
```

```
##      [,1]     [,2]  
## [1,] "hello"  "Help"
## [2,] "Please" "Me"
```

Now we can make our first function for this

```r
gene_intersect <- function(x, y) {
cbind(x[x %in% y], y[y %in% x])
}
```

Test it on x and y

```r
gene_intersect(x, y)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

Lets change to take input data frames rather than vectors

```r
gene_intersect2(df1, df2) 
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```

Looks good, this is our skateboard!

Add some flexibility for col name to match by... 

```r
gene_intersect3(df1, df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```

Make it look nicer- assign variables to the outputs- gene_intersect 4

Lets use the merge()function for this...

```r
merge(df1, df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```


