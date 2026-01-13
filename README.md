# LAI_regress
Efficient Tractor based Whole Genome Local Ancestry Regression 

## Installation

To install, you need to first download the repo:

```
https://github.com/eldronzhou/LAI_regress.git
```

You can then compile SDPR by running `make`. Make sure you have openblas and lapack installed.

## Running whole genome regression

```bash
./lai_regress -vcf path_to_phased.vcf.gz -msp path_to_lai.msp -pheno path_to_pheno.txt -covar path_to_covar.txt -out path_to_out.txt
```
