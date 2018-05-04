
# Can you improve this analysis code?
library(bio3d)
#call library of bio3d to look for accessions

plotprot <- function(x){
  #provide and input and generate an output (name, inputs, body)
  #name the function 'plotprot' to access unique accession
 
   s <- read.pdb(x)
  #read the pdb for the unique accession
  #assign variable 's' to consolidate repetitive commands

  sc <- trim.pdb(s, chain="A", elety="CA")
  #assign the variable 'sc' to trim the pdb fetch to the specific chain/carbons of interest
  
  plotb3(sc$atom$b, sse=s1.chainA, typ="l", ylab="Bfactor")
  #plot the trimmed pdb on a line graph with the y axis label 'Bfactor'
}

plotprot("1AKE")
#generate the output with the unique accession of interest
  
  
  
  
  
#lines(to overwrite) 
  #return(x) optional- can do when necessary- for the values 

 


s1 <- read.pdb(name) 

s1.chainA <- trim.pdb(s1, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```