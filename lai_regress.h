#include <iostream>
#include <fstream>
#include <vector>

#ifndef LAI_REGRESS_H
#define LAI_REGRESS_H

typedef struct {
    double *pheno;
    double *W;
    // for linear
    double *Qty;
    double *tau;
    double *work;
    // for logistic
    double *eta;
    double *mu;
    double *z;
    double *X;
    double *weight;
    double *XtWX;
    double *XtWz;
    double *beta;
    // other
    std::vector<size_t> ind_idx;
    size_t n_snp;
    size_t n_ind;
    size_t n_cov = 0;
    int lwork;
    int n_anc1 = 0;
    int n_anc2 = 0;
    std::vector<std::string> chr;
    std::vector<std::string> id;
    std::vector<std::string> pos;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
    std::vector<int> mac1;
    std::vector<int> mac2;
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

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, int anc, Dat *dat, const std::string &out, int glm);
