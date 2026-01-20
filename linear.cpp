#include "linear.h"
#include "lai_regress.h"
#include "cblas.h"
#include <limits>
#include <cmath>
#include "t_distribution.hpp"

using std::cout; using std::endl; using std::ifstream;
using std::string; using std::getline;

extern "C" {
    void dgeqrf_(int* m, int* n, double* A, int* lda, double* tau, double* work, int* lwork, int* info);
    void dormqr_(char* side, char* trans, int* m, int* n, int* k, double* A, int* lda,
                 double* tau, double* C, int* ldc, double* work, int* lwork, int* info);
}

void process_cov(Dat *dat) {
    int n = dat->n_ind;
    int p = dat->n_cov;

    int info;
    int lwork = -1;
    double wkopt;

    // W = QR
    double* tau = (double*) malloc(p * sizeof(double));
    dgeqrf_(&n, &p, dat->W, &n, tau, &wkopt, &lwork, &info);

    lwork = (int)wkopt;
    double* work = (double*) calloc(lwork, sizeof(double));

    dgeqrf_(&n, &p, dat->W, &n, tau, work, &lwork, &info);

    // Qt y
    double *Qty = (double*) calloc(dat->n_ind, sizeof(double));
    for (size_t i=0; i<n; i++) {
        Qty[i] = dat->pheno[i];
    }

    char side = 'L';   // Apply from the left
    char trans = 'T';  // Q^T
    int nrhs = 1;
    dormqr_(&side, &trans, &n, &nrhs, &p, dat->W, &n, tau, Qty, &n, work, &lwork, &info);

    // Q Qt y
    trans = 'N';
    for (size_t i=p; i<n; i++) {
	Qty[i] = 0;
    }
    dormqr_(&side, &trans, &n, &nrhs, &p, dat->W, &n, tau, Qty, &n, work, &lwork, &info);
    
    for (size_t i=0; i<n; i++) {
        Qty[i] = dat->pheno[i] - Qty[i];
    }

    dat->Qty = Qty;
    dat->work = work;
    dat->lwork = lwork;
    dat->tau = tau;
}

void linear(Dat *dat, double *X, double *Xr) {
    int n = dat->n_ind;
    int p = dat->n_cov;
    int lwork = dat->lwork;
    int info;

    // QtX overwrite X
    char side = 'L';   // Apply from the left
    char trans = 'T';  // Q^T
    int nrhs = 2;
    dormqr_(&side, &trans, &n, &nrhs, &p, dat->W, &n, dat->tau, X, &n, dat->work, &lwork, &info);
    
    // Q QtX overwrite QtX
    trans = 'N';
    for (size_t i=p; i<n; i++) {
        X[i] = 0;
	X[n+i] = 0;
    }
    dormqr_(&side, &trans, &n, &nrhs, &p, dat->W, &n, dat->tau, X, &n, dat->work, &lwork, &info);

    // Xr = X - Q Qt X
    double tmp = 0;
    for (size_t i=0; i<2*n; i++) {
	tmp = Xr[i];
        Xr[i] -= X[i];
	X[i] = tmp;
    }

    // Xr^t Xr \beta = Xr^t yr
    double XtX[4] = {0,0,0,0};
    double Xty[2] = {0,0}; 
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 2, n,
            1.0, Xr, n, 0.0, XtX, 2);
    cblas_dgemv(CblasColMajor, CblasTrans, n, 2,
            1.0, Xr, n, dat->Qty, 1, 0.0, Xty, 1);

    // Solve for beta1 and beta2
    double a = XtX[0]; double b = XtX[2]; // Note column major and upper
    double c = b; double d = XtX[3];
    double e = Xty[0]; double f = Xty[1];
    
    double det = a*d - b*c;

    // check mac
    double mac1 = 0, mac2 = 0;
    for (size_t i=0; i<n; i++) {
	mac1 += X[i];
    }
    for (size_t i=0; i<n; i++) {
        mac2 += X[n+i];
    }
    
    if (mac1 > n / 2) {
        mac1 = n - mac1;
    }

      
    if (mac1 > n / 2) {
        mac2 = n - mac2;
    }

    double beta1 = 0, beta2 = 0;
    if (mac1 < 20) {
	beta2 = f / d;
    }
    else if (mac2 < 20) {
	beta1 = e / a;
    }
    else {
        beta1 = (d*e - b*f) / det;
	beta2 = (-c*e + a*f) / det;
    }

    double sigma2 = 0, r = 0;
    for (size_t i=0; i<n; i++) {
	r = dat->Qty[i] - X[i] * beta1 - X[i+n] * beta2;
	sigma2 += r * r;
    }

    double se1 = 0, se2 = 0, nu = 0;
    if (mac1 < 20) {
        sigma2 /= n - p - 1;
	se2 = sqrt(sigma2 / d);
	se1 = std::numeric_limits<double>::quiet_NaN();
	beta1 = std::numeric_limits<double>::quiet_NaN();
	nu = n - p - 1;
    }
    else if (mac2 < 20) {
	sigma2 /= n - p - 1;
	se1 = sqrt(sigma2 / a);
	se2 = std::numeric_limits<double>::quiet_NaN();
	beta2 = std::numeric_limits<double>::quiet_NaN();
	nu = n - p - 1;
    }
    else {
        sigma2 /= n - p - 2;
	se1 = sqrt(sigma2 * d / det);
	se2 = sqrt(sigma2 * a / det);
	nu = n - p - 2;
    }

    double t1 = beta1 / se1;
    double t2 = beta2 / se2;
    double p1 = t_distribution::two_tailed_p(t1, nu);
    double p2 = t_distribution::two_tailed_p(t2, nu);
    double diff = beta1 - beta2;
    double se_diff = sqrt(sigma2 * (a + d + b + c) / det);
    double z_diff = diff / se_diff;
    double p_diff = t_distribution::chi2_pvalue(z_diff*z_diff, 1);

    dat->beta1.push_back(beta1);
    dat->se1.push_back(se1);
    dat->beta2.push_back(beta2);
    dat->se2.push_back(se2); 
    dat->p1.push_back(p1);
    dat->p2.push_back(p2);
    dat->se_diff.push_back(se_diff); 
    dat->p_diff.push_back(p_diff);
}


