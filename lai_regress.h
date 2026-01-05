#include <iostream>
#include <fstream>
#include <vector>

#ifndef LAI_REGRESS_H
#define LAI_REGRESS_H

typedef struct {
    double *pheno;
    double *W;
    double *Qty;
    double *QtX;
    double *XtX;
    double *Xty;
    double *tau;
    double *work;
    std::vector<size_t> ind_idx;
    size_t n_snp;
    size_t n_ind;
    size_t n_cov = 0;
    int lwork;
    std::vector<std::string> chr;
    std::vector<std::string> id;
    std::vector<std::string> pos;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
    std::vector<double> beta1;
    std::vector<double> beta2;
    std::vector<double> se1;
    std::vector<double> se2;
    std::vector<double> p1;
    std::vector<double> p2;
    std::vector<double> se_diff;
    std::vector<double> p_diff;
} Dat;

#endif

void get_size_vcf(const std::string &pheno_path, const std::string &geno_path, Dat *dat);

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, int anc, Dat *dat, const std::string &out);
