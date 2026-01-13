# LAI_regress
Efficient Tractor based Whole Genome Local Ancestry Regression 

## Installation

To install, you need to first download the repo:

```
https://github.com/eldronzhou/LAI_regress.git
```

You can then compile LAI_regress by running `make`. Make sure you have openblas and lapack installed.

## Running whole genome regression

LAI_regress takes whole genome phased genotype and local ancestry files as the input, and perform local ancestry based association between genotype and phenotype. It scales efficiently on All of Us, capable of analyzing 9 million common variants for 70,000 African Americans within hours, at the cost of under $2 per trait. 

```bash
./lai_regress -vcf path_to_phased.vcf.gz -msp path_to_lai.msp -pheno path_to_pheno.txt -covar path_to_covar.txt -out path_to_out.txt
```

Below are the required options.

- vcf (required): Path to the phased genotype file in the gzipped vcf format. Not that currently missing genotypes are not supported. If you use Eagle2 for phasing, there shouldn't be missing genotype in the default setting.
- msp (required): Path to the whole genome solved local ancestry files.
- pheno (required): path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included.
- covar (required): path to the covariate file. Covariates will be reading from the first column. There is no header for the covariate file.
- out (required): Path to the output file.

The output file has the following format:

```
id      chr     pos     ref     alt     beta1   se1     beta2   se2     p1      p2      se_diff p_diff
rs576565926     10      46755   A       AC      0.955794        1.41235 0.361792        0.156033        0.498571        0.0204145       1.42103 0.675941
rs572697851     10      46868   AAAG    A       -0.314565       0.815897        0.108156        0.106708        0.699835        0.310795        0.824136        0.608003
...
```
where beta1, se1, p1 are effect sizes, standard error, p value for ancestry 1 based genotype; beta2, se2, p2 are effect sizes, standard error, p value for ancestry 2 based genotype; se_diff and p_diff are standard error and p value of the difference of two ancestral effects.
