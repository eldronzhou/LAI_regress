# LAI_regress
Efficient Tractor based Whole Genome Local Ancestry Regression 

## Installation

To install, you need to first download the repo:

```
https://github.com/eldronzhou/LAI_regress.git
```

You can then compile LAI_regress by running `make`. Make sure you have openblas and lapack installed.

## Running whole genome regression

LAI_regress takes phased whole genome genotype and local ancestry files as the input, and perform regression for local ancestry resolved association. It scales efficiently on All of Us, capable of analyzing 9 million variants for 75,000 African Americans within hours and the cost of fewer than $2 per trait. 

```bash
./lai_regress -vcf path_to_phased.vcf.gz -msp path_to_lai.msp -pheno path_to_pheno.txt -covar path_to_covar.txt -out path_to_out.txt
```
