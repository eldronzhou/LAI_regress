#include "lai_regress.h"
#include "linear.h"
#include "logistic.h"
#include <assert.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <zlib.h>
#include <cstring>
#include "cblas.h"
#include "math.h"

using std::cout; using std::endl; using std::ifstream;
using std::string; using std::getline;


// get n_snp and n_ind
void get_size_vcf(const string &pheno_path, const string &geno_path, Dat *dat) {
    size_t n_ind = 0;
    size_t n_invalid = 0;
    size_t n_snp = 0;

    string line;
    ifstream infile1(pheno_path.c_str());
    string id;
    string y;
    size_t i = 0;
    while (infile1 >> id >> id >> y) {
	try {
	    std::stod(y); 
	    dat->ind_idx.push_back(i);
	    n_ind++;
	}
	catch (std::invalid_argument&) {
	    n_invalid++;
	}
	i++;
    }
    dat->n_ind = n_ind;
    cout << "Warning: " + std::to_string(n_invalid) + \
	" individuals with invalid phenotypes." << endl;

    gzFile infile2 = gzopen(geno_path.c_str(), "rb");

    char buffer[4096];

    while (gzgets(infile2, buffer, 4096)) {
	// Remove the newline character if it exists
        size_t length = strlen(buffer);
	line.append(buffer);
        if (length > 0 && buffer[length - 1] != '\n') {
	    continue;
        }

	if (line.find("##") == 0) {
	    line.clear();
        }
        else if (line.find("#") == 0) {
	    line.clear();
        }
        else {
            line.clear();
	    n_snp++;
        }
    }

    gzclose(infile2);

    /*ifstream infile2(geno_path.c_str());
    while(getline(infile2, line)) {
	if (line.find("##") == 0) {
	    continue;
	}
	else if (line.find("#") == 0) {
	    continue;
	}
	else {
	    n_snp++;	
	}
    }*/

    dat->n_snp = n_snp;
    cout << "In total " + std::to_string(n_snp) + " SNPs and " \
	+ std::to_string(n_ind) + " individuals to be readed." << endl;
}

void read_pheno(const std::string &pheno_path, Dat *dat) {
    ifstream infile(pheno_path.c_str());

    cout << "Reading phenotype file from: " + pheno_path + "." << endl;

    string id;
    string y;
    double *pheno = (double *) malloc(dat->n_ind*sizeof(double));

    size_t i = 0, idx = 0;
    while (infile >> id >> id >> y) {
	if (i == dat->ind_idx[idx]) {
	    pheno[idx] = stod(y);
	    idx++;
	}
	i++;
    }
    dat->pheno = pheno;

    cout << "Readed phenotype from " + std::to_string(idx) + " individuals." << endl;
}

void read_cov(const std::string &cov_path,  Dat *dat) {
    ifstream infile(cov_path.c_str());

    cout << "Reading covariate file from: " + cov_path + "." << endl;

    size_t n_cov = 0, i = 0, idx = 0;

    string line, token;

    while (getline(infile, line)) {
	n_cov = 0;
	std::istringstream iss(line);
	if (i == 0) {
	    while (getline(iss, token, '\t')) {
		n_cov++;
	    }
	    cout << "Reading " + std::to_string(n_cov) + " covariates." << endl;
	    dat->n_cov = n_cov;
	    dat->W = (double *) malloc(dat->n_ind*n_cov*sizeof(double));
	    n_cov = 0;
	    iss.clear();
	    iss.str(line);
	}
	if (i == dat->ind_idx[idx]) {
	    while (getline(iss, token, '\t')) {
		dat->W[n_cov*dat->n_ind+idx] = stod(token); // n by p but lapack uses column major
		n_cov++;
	    }
	    idx++;
	}
	i++;
    }
    infile.close();
}

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, int anc, Dat *dat, const std::string &out_path, int glm) {
    
    string line1, line2;
    string token1, token2;

    ifstream mspfile(msp_path.c_str());
    cout << "Reading RFmix msp file from: " + msp_path + "." << endl;
    
    //ifstream infile(vcf_path.c_str());
    //cout << "Reading VCF file from: " + vcf_path + "." << endl;
    gzFile infile = gzopen(vcf_path.c_str(), "rb");
    char buffer[4096];

    // skip first two lines of msp file
    getline(mspfile, line1);
    getline(mspfile, line1);

    // skip the header of vcf file
    int n_ind = 0;
    while (gzgets(infile, buffer, 4096)) {
        line2.append(buffer);
	size_t length = strlen(buffer);
        if (length > 0 && buffer[length - 1] != '\n') {
	    continue;
        }

	if (line2.find("##") == 0) {
	    line2.clear();
	    continue;
	}
	else if (line2.find("#") == 0) {
	    int idx_2 = 0;
	    std::istringstream iss2(line2);
	    while (getline(iss2, token2, '\t')) {
		idx_2++;
	    }
	    n_ind = idx_2 - 9;
	    line2.clear();
	}
	else {
	    break;
	}
    }

    int *hap_lanc = (int *) malloc(2*n_ind*sizeof(int));
    unsigned spos = 0, epos = 0, pos = 0;
    int chr_vcf = 1, chr_msp = 1;
    
    // read the pos from the first line of vcf file
    std::istringstream iss2(line2);
    int idx2 = 0;
    for (; idx2<9; idx2++) {
	getline(iss2, token2, '\t');
	if (idx2 == 0) {
	   if (token2.find("chr") == 0) {
               chr_vcf = std::stoi(token2.substr(3));
	   }
	   else { 
	       chr_vcf = std::stoi(token2);
	   }
	}
	if (idx2 == 1) {
	    pos = std::stoul(token2);
	}
    }

    size_t idx_snp = 0;
    double *X = (double *) calloc(2*dat->n_ind, sizeof(double));
    double *Xr = (double *) calloc(2*dat->n_ind, sizeof(double));

    while (getline(mspfile, line1)) {
	// read msp file
	std::istringstream iss1(line1);
	for (int idx1=0; idx1<2*n_ind+6; idx1++) {
	    getline(iss1, token1, '\t'); 
	    if (idx1 == 0) {
	        if (token1.find("chr") == 0) {
	            chr_msp = std::stoi(token1.substr(3));
		}
		else {
		   chr_msp = std::stoi(token1);
		}
	    }
	    else if (idx1 == 1) {
		spos = std::stoul(token1);
	    }
	    else if (idx1 == 2) {
		epos = std::stoul(token1);
	    }
	    else if (idx1 >= 6) {
		hap_lanc[idx1-6] = std::stoi(token1);
		if (hap_lanc[idx1-6] != 0 && hap_lanc[idx1-6] != 1 && hap_lanc[idx1-6] != 2) {
		    cout << "MSP field must be either 0, 1 or 2." << endl;
		    return;
		}
	    }
	}

	if ((chr_vcf != chr_msp) || (pos < spos) || (pos > epos)) {
	    cout << "Inconsistent starting position: chr_vcf: " + std::to_string(chr_vcf) + \
		" chr_msp: " + std::to_string(chr_msp) + " pos: " + \
		std::to_string(pos) + " spos: " + std::to_string(spos) + " epos: " + std::to_string(spos) << endl;
	    exit(EXIT_FAILURE);
	} 
	     
	// read vcf file
	while ((chr_vcf == chr_msp && pos >= spos && pos <= epos) || idx_snp == dat->n_snp-1) {
	    // reset the stream
	    std::istringstream iss2(line2);

	    size_t k = 0;
	    for (idx2=0; idx2<n_ind+9; idx2++) {
		getline(iss2, token2, '\t');

		if (idx2 == 0) {
		    dat->chr.push_back(token2);
		}

		if (idx2 == 1) {
		    dat->pos.push_back(token2);
		}

		if (idx2 == 2) {
		    dat->id.push_back(token2);
		}

		if (idx2 == 3) {
		    dat->ref.push_back(token2);
		}

		if (idx2 == 4) {
		    dat->alt.push_back(token2);
		}

		if (idx2 >= 9) {

		    // individuals kept in analysis
		    if (idx2-9 != dat->ind_idx[k]) {
			continue;
		    }

		    // check phasing
		    if (token2[1] != '|') {
			cout << "Genotype must be phased." << endl;
			return;
		    }
		    
		    // read the genotype
		    if (token2[0] == '.' || token2[3] == '.') {
			cout << "Missing genotype not supported yet." << endl;
			return;
		    }
		    if (hap_lanc[2*(idx2-9)] == anc) {
			X[k] += std::stod(&token2[0]);
		    }
		    else {
			X[dat->n_ind+k] += std::stod(&token2[0]);
		    }

		    if (hap_lanc[2*(idx2-9)+1] == anc) {
			X[k] += std::stod(&token2[2]);
		    }
		    else  {
			X[dat->n_ind+k] += std::stod(&token2[2]);
		    }
		    k++;
		}
	    }

	    if (idx_snp == 123429) {
	       std::ofstream out("./X.txt");
	       for (size_t i=0; i<2*dat->n_ind; i++) {
		       out << X[i] << " " << endl;
	       }
	       cout << "Done!" << endl;
	       out.close();
	       /*memcpy(Xr, X, 2*dat->n_ind*sizeof(double));
	       linear(dat, X, Xr);
	       return;*/
	    }

	    //regress here
	    if (glm) {
	        logistic(dat, X);
	    }
	    else {
		memcpy(Xr, X, 2*dat->n_ind*sizeof(double));
	        linear(dat, X, Xr);
	    }

	    if (idx_snp % 10000 == 0) {
               cout << idx_snp << endl;
	    }

	    // read the next line of vcf and update pos
	    assert(idx2 == n_ind+9);
	  
	    memset(X, 0, 2*dat->n_ind*sizeof(double));
	    line2.clear();

	    if (gzgets(infile, buffer, 4096)) {
		size_t length = strlen(buffer);
		line2.append(buffer);
		while (length > 0 && buffer[length - 1] != '\n') {
		    gzgets(infile, buffer, 4096);
		    length = strlen(buffer);
		    line2.append(buffer);
		}

		iss2.clear();
		iss2.str(line2);
		for (idx2=0; idx2<9; idx2++) {
		    getline(iss2, token2, '\t');
		    if (idx2 == 0) {
			if (token2.find("chr") == 0) {
			   chr_vcf = std::stoi(token2.substr(3));
			}
			else {
		            chr_vcf = std::stoi(token2);
			}
		    }
		    if (idx2 == 1) {
			pos = std::stoul(token2);
		    }
		}
		idx_snp++;
	    }
	    else {
		break;
	    }
	}	
    }
    
    cout << "Analyzed " << std::to_string(idx_snp+1) << \
	" SNPs from " << std::to_string(dat->n_ind) << " individuals." << endl;
    free(hap_lanc);
    gzclose(infile);

    free(X);
    free(Xr);

    // write to the output
    std::ofstream out(out_path);
    out << "id" << "\t" << "chr" << "\t" << "pos" << "\t" << "ref" << "\t" << "alt" << "\t" \
            << "mac1" << "\t" << "mac2" << "\t" << "beta1" << "\t" << "se1" << "\t" << "beta2" << "\t" << "se2" << "\t" \
            << "p1" << "\t" << "p2" << "\t" << "se_diff" << "\t" << "p_diff" << endl;
    for (size_t i=0; i<dat->n_snp; i++) {    
        out << dat->id[i] << "\t" << dat->chr[i] << "\t" << dat->pos[i] << "\t" << dat->ref[i] << "\t" \
        << dat->alt[i] << "\t" << dat->mac1[i] << "\t" << dat->mac2[i] << "\t" << dat->beta1[i] << "\t" << \
	dat->se1[i] << "\t" << dat->beta2[i] \
	<< "\t" << dat->se2[i] << "\t" << dat->p1[i] << "\t" << dat->p2[i] << "\t" << \
	dat->se_diff[i] << "\t" << dat->p_diff[i] << endl;
    }
    out.close();
}

int main(int argc, char *argv[]) {

    std::string vcf_path, pheno_path, covar_path, msp_path, out_path;

    int i = 1, anc = 0, glm = 0;
    while (i < argc) {
        if (strcmp(argv[i], "-vcf") == 0) {
            vcf_path = argv[i+1];
            i += 2;
        }
	else if (strcmp(argv[i], "-pheno") == 0) {
	    pheno_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-covar") == 0) {
	    covar_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-anc") == 0) {
	   anc = std::stoi(argv[i+1]);
	   i += 2;
	}
        else if (strcmp(argv[i], "-msp") == 0) {
            msp_path = argv[i+1];
            i += 2;
        }
        else if (strcmp(argv[i], "-glm") == 0) {
	    glm = 1;
	    i += 1;
	}
	else if (strcmp(argv[i], "-out") == 0) {
            out_path = argv[i+1];
            i += 2;
        }
        else {
            cout << "Invalid option: " << argv[i] << endl;
            return 0;
        }
    }

    Dat dat;

    get_size_vcf(pheno_path.c_str(), vcf_path.c_str(), &dat);
    
    read_pheno(pheno_path.c_str(), &dat);

    read_cov(covar_path.c_str(), &dat);

    if (glm) {
	fit_cov(&dat);
    }
    else {
        process_cov(&dat);
    }

    read_lanc(vcf_path.c_str(), msp_path.c_str(), anc, &dat, out_path.c_str(), glm);

    free(dat.pheno);
    free(dat.W);
    
    if (!glm) {
        free(dat.Qty);
        free(dat.tau);
        free(dat.work);
        free(dat.QtX);
        free(dat.XtX);
        free(dat.Xty);
        free(dat.tau);
        free(dat.work);
    }
    else {
	free(dat.eta);
	free(dat.mu);
	free(dat.weight);
	free(dat.z);
	free(dat.X);
	free(dat.XtWX);
	free(dat.XtWz);
	free(dat.beta);
    }
    return 0;
}
