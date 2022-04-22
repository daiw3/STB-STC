# STB-STC
Code for Super-taxon in Human Microbiome are Identified to Be Associated with Colorectal Cancer

## Core Code
code/SVB-SVC.R is the code for performing the method (STB or STC) on a group of OTUs.
code/utilies includes all necessary codes to support the running of SVB-SVC.R

## Input Data Preparation
SVB-SVC.R requires the OTUs data to be partitioned into blocks based on hiearchical information (Genes, Family, Order, Class) beforehead and save each block as a file.
Please see data/family_STB_blk1.RData for an example to see the format of paritioned OTUs based on family information .

## Excute the code
Run.R take inputs of the block index (--blk), hiearchical level (--level) and method (--STB or STC) to generate results.

## Output file
The output results are saved in a list in the format .rda file.
- super.var.form: gives the cut-off resutls
- super.var.dis: gives the marignal assocation results, including effect sizes, p-values
- otus: gives the selected OTUs
- X: the formed Super-taxon across individuals
- Z: depth importance scores of OTUs used for ranking
- geno.transed: candidate cut-offs used for selecting best cut-off

