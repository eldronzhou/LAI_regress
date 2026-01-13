# LAI_regress
Efficient Tractor based Whole Genome Local Ancestry Regression 

## Installation

To install, you need to first download the repo:

```
https://github.com/eldronzhou/LAI_regress.git
```

You can then compile LAI_regress by running `make`. Make sure you have openblas and lapack installed.

## Running whole genome regression

LAI_regress takes whole genome phased genotype and local ancestry files as the input, and perform local ancestry based association between genotype and phenotype. It scales efficiently on All of Us, capable of analyzing 9 million common variants for 70,000 African Americans within hours, at the cost of fewer than $2 per trait. 

```bash
./lai_regress -vcf path_to_phased.vcf.gz -msp path_to_lai.msp -pheno path_to_pheno.txt -covar path_to_covar.txt -out path_to_out.txt
```

Below are the required options.

- vcf (required): Path to the phased genotype file in the gzipped vcf format. Not that currently missing genotypes are not supported. If you use Eagle2 for phasing, there shouldn't be missing genotype in the default setting.
- msp (required): Path to path to the whole genome solved local ancestry files.
- pheno (required): path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included.
- covar (required): path to the covariate file. Covariates will be reading from the first column. There is no header for the covariate file.
- out (required): Path to the output file.
