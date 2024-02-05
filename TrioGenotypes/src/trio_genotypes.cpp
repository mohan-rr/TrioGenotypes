#include <string>
#include <vector>
#include <fstream>
#include <fstream>
#include <Rcpp.h>

void open_bed_file(std::string &bed_filename, std::ifstream& bed_file, int n_ind, int n_snp);

Rcpp::List count_genotypes(std::string bed_filename, std::vector<int> trio_c, std::vector<int> trio_f, std::vector<int> trio_m, std::vector<int> trio_pheno, int n_ind, int n_snp){
	std::ifstream bed_file;
	open_bed_file(bed_filename, bed_file, n_ind, n_snp);

	// snp-major .bed files store data contiguously with 4 individual genotypes per byte. the last byte is padded.
	// the order of individuals are in the same order as in the .fam file
	int nbytes = n_ind/4 + (n_ind % 4 != 0);
	int remain = n_ind % 4;

    //genotype combinations (M,F,C)
	Rcpp::IntegerMatrix R_case_trio_counts (15, n_snp);
	Rcpp::IntegerMatrix R_ctrl_trio_counts (15, n_snp);

	char buffer[nbytes];
	char read_bit;
	int mother_geno, father_geno, child_geno;
	unsigned int genotype;
	std::vector<char> geno_vec;
	geno_vec.reserve(n_ind);

	for(int h = 0; h < n_snp; h++){
		bed_file.read(buffer, nbytes);
		for(int i = 0; i < nbytes - 1; i++){
			read_bit = buffer[i];
			for (int j = 0; j < 4; j++){
				genotype = read_bit & (char)3;
				read_bit = read_bit >> 2;
				geno_vec[i*4 + j] = genotype;
			}
		}
		read_bit = buffer[nbytes-1];
		for (int j = 0; j < remain; j++){
			genotype = read_bit & (char)3;
			read_bit = read_bit >> 2;
			geno_vec[(nbytes-1)*4 + j] = genotype;
		}
		
		Rcpp::Rcout << "trio_c.size() " << trio_c.size() << "\n" ; 
		//father, mother, child. 4 genotype codes from bed file : 2 homo, hetero, missing
		// 0 - A1 homozygote , 1 - heterozygote , 2 - missing , 3 - A2 homozygote
		// we will use A1 as the risk allele
		int ctrl_geno_count[4][4][4] = {0};
		int case_geno_count[4][4][4] = {0}; 

		for (size_t i = 0; i < trio_c.size(); i++)
		{
			mother_geno = geno_vec[trio_m[i]];
			father_geno = geno_vec[trio_f[i]];
			child_geno  = geno_vec[trio_c[i]];
			if (trio_pheno[i] == 1) ctrl_geno_count[mother_geno][father_geno][child_geno] ++ ;
			if (trio_pheno[i] == 2) case_geno_count[mother_geno][father_geno][child_geno] ++ ;
		}
		//list the 15 combinations 
		// don't like typing them all out but it's branchless, 
		// could have calculated the index but the missing index is 2 which throws it off
		R_case_trio_counts(h, 0)  = case_geno_count[0][0][0];
		R_case_trio_counts(h, 1)  = case_geno_count[0][1][0];
		R_case_trio_counts(h, 2)  = case_geno_count[0][1][1];
		R_case_trio_counts(h, 3)  = case_geno_count[1][0][0];
		R_case_trio_counts(h, 4)  = case_geno_count[1][0][1];
		R_case_trio_counts(h, 5)  = case_geno_count[0][3][1];
		R_case_trio_counts(h, 6)  = case_geno_count[3][0][1];
		R_case_trio_counts(h, 7)  = case_geno_count[1][1][0];
		R_case_trio_counts(h, 8)  = case_geno_count[1][1][1];
		R_case_trio_counts(h, 9)  = case_geno_count[1][1][3];
		R_case_trio_counts(h, 10) = case_geno_count[1][3][1];
		R_case_trio_counts(h, 11) = case_geno_count[1][3][3];
		R_case_trio_counts(h, 12) = case_geno_count[3][1][1];
		R_case_trio_counts(h, 13) = case_geno_count[3][1][3];
		R_case_trio_counts(h, 14) = case_geno_count[3][3][3];

		R_ctrl_trio_counts(h, 0)  = ctrl_geno_count[0][0][0]; 
		R_ctrl_trio_counts(h, 1)  = ctrl_geno_count[0][1][0];
		R_ctrl_trio_counts(h, 2)  = ctrl_geno_count[0][1][1]; 
		R_ctrl_trio_counts(h, 3)  = ctrl_geno_count[1][0][0]; 
		R_ctrl_trio_counts(h, 4)  = ctrl_geno_count[1][0][1];
		R_ctrl_trio_counts(h, 5)  = ctrl_geno_count[0][3][1]; 
		R_ctrl_trio_counts(h, 6)  = ctrl_geno_count[3][0][1]; 
		R_ctrl_trio_counts(h, 7)  = ctrl_geno_count[1][1][0];
		R_ctrl_trio_counts(h, 8)  = ctrl_geno_count[1][1][1]; 		
		R_ctrl_trio_counts(h, 9)  = ctrl_geno_count[1][1][3]; 
		R_ctrl_trio_counts(h, 10) = ctrl_geno_count[1][3][1];
		R_ctrl_trio_counts(h, 11) = ctrl_geno_count[1][3][3]; 
		R_ctrl_trio_counts(h, 12) = ctrl_geno_count[3][1][1]; 
		R_ctrl_trio_counts(h, 13) = ctrl_geno_count[3][1][3];
		R_ctrl_trio_counts(h, 14) = ctrl_geno_count[3][3][3]; 
	}
	return Rcpp::List::create( 
	  Rcpp::Named("case")  = R_case_trio_counts, 
	  Rcpp::Named("ctrl")  = R_ctrl_trio_counts 
	) ;
}


void open_bed_file(std::string &bed_filename, std::ifstream& bed_file, int n_ind, int n_snp){
		
	bed_file.open(bed_filename.c_str(), std::ifstream::binary);
	if(!bed_file.is_open()) Rcpp::stop("\nCannot open .bed file");
	
	char buffer[3];
	bed_file.read(buffer, 3);
	int magicNo1 = buffer[0];
	int magicNo2 = buffer[1];
	int mode = buffer[2];

	if(magicNo1 != 108 || magicNo2 != 27){
	Rcpp::stop("Detected an old version .bed file, please use a new(er) version of plink and the --make-bed command");
	}
	if(mode != 1)
	{
	Rcpp::stop("currently only snp-major mode .bed files are supported. please use plink --make-bed to output a snp-major .bed file");		
	}

	//the bed file should be ciel(nind/4)*n_snp bytes + 3 bytes long
	// add in check to see that filesize is correct. 
	// ifstream::tellg 
}