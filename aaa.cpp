#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <cstdint>

struct trios{
	std::vector<unsigned int> ind, mother, father;
}

struct snpinfo{
	std::vector<std::string> snpID, a1, a2;
}

read_bim(std::string &bim_filename, snpinfo& snpinfo_){
	std::ifstream bim_file;
	bim_file.open(bim_filename.c_str());
	if(!bim_file.is_open())	Rcpp::stop("can't open .bim file");
	
	std::string chr, snpID, cM, BP, a1, a2;
	while(!bim_file.eof()){
		bim_file >> chr >> snpID >> cM >> BP >> a1 >> a2;
		snpinfo_.snpID.push_back(snpID);
		snpinfo_.a1.push_back(a1);
		snpinfo_.a2.push_back(a2);
	}

	Rcpp::message()
}

count_genotypes(const std::vector<unsigned char> & geno_vec, trios & case_trios){

}

main(std::string bed_filename, std::string bim_filename, std::string fam_filename, nind, nsnp){
	read_fam()
	
	read_bim()
	read_bed()
	uint nbytes = nind / 4;
	uint remain = nind % 4;

	for(int h = 0; h < nSNP; h++){
		for(int i = 0; i < nbytes; i++){
			bed_file.read(buffer, 1);
			read_bit = buffer[0];
			for (int j = 0; j < 4; j++){
				genotype = read_bit & (char)3;
				read_bit = read_bit >> 2;
				geno_vec.push_back(genotype);
			}
		}
		bed_file.read(buffer,1);
		read_bit = buffer[0];
		for (int i = 0; i < remain; i++){
			genotype = read_bit & (char)3;
			read_bit = read_bit >> 2;
			geno_vec.push_back(genotype);
		}
		for i in case_trios.ind.length();
		{
			mother_geno = geno_vec[case_trios.mother[i]];
			father_geno = geno_vec[case_trios.father[i]];
			child_geno  = geno_vec[case_trios.ind[i]];
			geno_count[father_geno][mother_geno][child_geno] ++ ;
		}
	}
}


inline void read_bed(std::string &bed_filename, std::ifstream& bed_file, const uint nSNP, const uint nind){
		
	char buffer[3];
	unsigned char one = 1;
	unsigned char read_bit, allele1, allele2;
	uint snp_count = 0;
	uint bitCount;
	uint indCount = 0;
	bed_file.open(bed_filename.c_str(), ios::binary);
	if(!readBinaryGenotypeData.is_open()) Rcpp::stop("\nCannot open .bed file");
	
	bed_file.read(buffer, 3);
	uint magicNo1 = buffer[0];
	uint magicNo2 = buffer[1];
	uint mode = buffer[2];

	if(magicNo1 != 108 || magicNo2 != 27){
	Rcpp::stop("Detected an old version .bed file, please use a new(er) version of plink and the --make-bed command");
	}
	if(mode != 1)
	{
	Rcpp::stop("currently only snp-major mode .bed files are supported. please use plink --make-bed to output a snp-major .bed file");		
	}
}


inline void read_fam(std::string& fam_filename, unsigned int &nind, trios &case_trios, trios &ctrl_trios){

	std::ifstream fam_file;
	fam_file.open(fam_filename.c_str());
	if(!fam_file.is_open())	Rcpp::stop("can't open .fam file");

	std::string fam_ID, IID, FID, MID, sex, pheno;
	std::vector<std::string> v_famID, v_FID, v_MID, v_pheno;
	std::map<std::string, unsigned int> IID_ord;
	v_famID.reserve(nind); v_FID.reserve(nind); v_MID.reserve(nind); v_pheno.reserve(nind);
	uint i = 0;
	do{
		fam_file >> fam_ID >> IID >> FID >> MID >> sex >> pheno;
		IID_ord.emplace(fam_ID + IID, i);
		v_famID.push_back(fam_ID);
		v_FID.push_back(FID);
		v_MID.push_back(MID);
		v_pheno.push_back(pheno);
		i++;
	}while(!fam_file.eof());

	std::vector<uint> trio_ind, trio_m, trio_f;
	trio_ind.reserve(nind/2); trio_m.reseve(nind/2); trio_f.reserve(nind/2);
	map<std::string, uint>::iterator search_ID;  //search for the ordinal ID of father/mother

	for(uint j = 0; j <= i ; j++ ){
		if (FID[j] != "0" && MID[j] != "0"){
			trio_ind.push_back(j);
			search_ID = IID_ord.find(v_famID[j] + v_FID[j]);
			if (search_ID == IID_ord.end()) {
				Rcpp::stop("father ID not not in IID");
			}
			trio_f.push_back(search_ID->second);
			search_ID = IID_ord.find(v_famID[j] + v_MID[j]);
			if (search_ID == IID_ord.end()) {
				Rcpp::stop("mother ID not not in IID");
			}
			trio_m.push_back(search_ID->second)	;	
		}
	}
}
