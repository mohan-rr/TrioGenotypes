#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

struct trio{
	int ind, mother, father;
	trio(int ind, int mother, int father) 
		: ind(ind), mother(mother), father(father) {}
};

class read_bim
{
private: 
	std::vector<std::string> v_snpID;
	std::vector<std::string> allele1;
	std::vector<std::string> allele2;
	std::string chr, snpID, cM, BP, a1, a2;
public:
	read_bim(std::string bim_filename, int nsnp)
	read_line()
}

class read_fam
{

}
read_bim(std::ifstream &bim_file){
	bim_file >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >>alleleName1 >> alleleName2;
}

read_fam(std::ifstream &fam_file, std::string& fam_ID, std::string& IID, std::string& FID, std::string& MID, std::string& sex, std::string& pheno){

IID_ord std::map<std::string, uint>;
uint i = 0;
do{
	i++;
	fam_file >> fam_ID >> IID >> FID >> MID >> sex >> pheno;
	IID_ord.emplace(fam_ID + IID, i);
}while(!fam_file.eof());

std::vector<uint> trio_ind, trio_m, trio_f;
trio_ind.reserve(nind/2); trio_m.reseve(nind/2); trio_f.reserve(nind/2);
for(uint j = 0; j < i ; j++ ){

}
//return the 3 trio vectors, and the pheno vector




//read line
// while not end of file
// i++ 
ind_order.insert(fam_ID+IID, i)