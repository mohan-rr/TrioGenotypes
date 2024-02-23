#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// count_genotypes
Rcpp::List count_genotypes(std::string bed_filename, std::vector<int> trio_c, std::vector<int> trio_f, std::vector<int> trio_m, std::vector<int> trio_pheno, int n_ind, int n_snp);
RcppExport SEXP _TrioGenotypes_count_genotypes(SEXP bed_filenameSEXP, SEXP trio_cSEXP, SEXP trio_fSEXP, SEXP trio_mSEXP, SEXP trio_phenoSEXP, SEXP n_indSEXP, SEXP n_snpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bed_filename(bed_filenameSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type trio_c(trio_cSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type trio_f(trio_fSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type trio_m(trio_mSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type trio_pheno(trio_phenoSEXP);
    Rcpp::traits::input_parameter< int >::type n_ind(n_indSEXP);
    Rcpp::traits::input_parameter< int >::type n_snp(n_snpSEXP);
    rcpp_result_gen = Rcpp::wrap(count_genotypes(bed_filename, trio_c, trio_f, trio_m, trio_pheno, n_ind, n_snp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TrioGenotypes_count_genotypes", (DL_FUNC) &_TrioGenotypes_count_genotypes, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_TrioGenotypes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
