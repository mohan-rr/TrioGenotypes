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

// struct trio{
// 	int ind, mother, father;
// 	trio(int ind, int mother, int father) 
// 		: ind(ind), mother(mother), father(father) {}
// };

// class read_bim
// {
// private: 
// 	std::vector<std::string> v_snpID;
// 	std::vector<std::string> allele1;
// 	std::vector<std::string> allele2;
// 	std::string chr, snpID, cM, BP, a1, a2;
// public:
// 	read_bim(std::string bim_filename, int nsnp);
// 	read_line();
// }
// class read_fam
// {

// }
read_bim(std::ifstream &bim_filename){
	bim_file >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >>alleleName1 >> alleleName2;
}

count_genotypes(){

}
read_single_snp()
trio_genotype(std::vector<int8_t> & geno_vec, )
read_bed(std::string &bed_filename, std::ifstream& bed_file, const uint& nSNP){
		
	char buffer[3];
	unsigned char one = 0x01;
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
	Rcpp::stop("currently only snp-major mode .bed files are supported. please use plink --make-bed to output a snp-major .bed file")		
	}

	bitCount = 9;
	while(snp_count < nSNP){
		while(indCount < nind){	
		if(bitCount == 9)
		{
			bed_file.read(buffer, 1);
			if(bed_file.eof()) Rcpp::stop("Error: .bed file is incomplete!\n");
			read_bit = buffer[0];
			bitCount = 1;
		}

		allele1  = read_bit & one; //read the least significant bit				
		read_bit = read_bit >> 1; //shift bits to the right
		allele2  = read_bit & one; //read the new least significant bit				
		read_bit = read_bit >> 1; //shift bits to the right for next time

		bitCount += 2;	

		if(allele1 == 1 && allele2 == 1)
		{
			allele1 = 0; allele2 = 0;
		}
		else if(allele1 == 0 && allele2 == 0)
		{
			allele1 = 1; allele2 = 1;
		};

		}
		snp_count++;
	}
}
//read bim file to put chr, snpID, alleles into the returned dataframe

//main loop will go through the bed file. it's in SNP major mode ... do we want to keep in snp major orientation?
//or read into subject vectors
//the counting loop will need to go through each trio ... but we can do the loop through snps actually ... 
//loop through SNPs, then for each SNP use the indices of the trios to do the counting
// then only need to keep one vector for the SNPs and overwrite for each marker
// so read each line of bed file into the marker vector. can do int8. -1 for missing, then 0/1/2 for rest
//there must be a way to do this branchless with just a binary operation, can just do elif ladder for now 
// ---> later can do a bitwise mask of M | F | C after shifting and it will give the 15 unique combinatinos
//      get last 2 bits of each, then shift them so you have MFC, then just take that number and it will 
//      for now we do ifelse ladder
// but also need an if statement for missing
void read_fam(std::string& fam_filename, uint nind, std::vector<uint>& trio_ind, std::vector<uint>& trio_m, std::vector<uint>& trio_f
				std::vector<int>& trio_ph){

	std::ifstream fam_file;
	fam_file.open(fam_filename.c_str());
	if(!fam_file.is_open())	Rcpp::stop("can't open .fam file");

	std::string fam_ID, IID, FID, MID, sex, pheno;
	std::vector<std::string> v_famID, v_FID, v_MID, v_pheno;
	std::map<std::string, uint> IID_ord;
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
