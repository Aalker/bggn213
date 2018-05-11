Bioinformatics Class 12
================
Amanda Alker
5/11/2018

Call the package

``` r
library(bio3d)
```

Setup HIV-Pr for docking study
==============================

Get the protein

``` r
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

Read this file in and trim out the protein and small molecule ligand from everything else

``` r
hiv <- read.pdb(file.name)
#./ means that it lives locally
```

``` r
ligand <- trim.pdb(hiv, "ligand")
ligand
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "ligand")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

Lets do the same for the protein Extract protein

``` r
protein <- trim.pdb(hiv, "protein")
protein
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "protein")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

Write them up

``` r
write.pdb(ligand, "1hsg_ligand.pdb")
write.pdb(protein, "1hsg_protein.pdb")
```

Convert the docking results from pdbqt file to .pdb in 'Bio3d'
--------------------------------------------------------------

``` r
res <- read.pdb("all.pdbqt", multi = TRUE)

write.pdb(res, "results.pdb")
```

Compare our results to Merk Co Crystal structure

``` r
ori <- read.pdb("ligand.pdbqt")
```

``` r
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

> **Q6.** RMSD based on non hydrogen atoms

``` r
#Get protein (PDB, "protein")
inds <- atom.select(ori, "noh")

rmsd(ori$xyz[,inds$xyz], res$xyz[ , inds$xyz])
```

    ##  [1]  0.458 11.021 10.374  4.301 10.891  3.717  5.764  3.791  5.498 10.759
    ## [11]  4.224  6.308 10.889  8.776
