
trio_genotypes <- function(filepath){
    fam_file <- normalizePath(paste0(filepath,".fam"), mustWork = TRUE)
    bim_file <- normalizePath(paste0(filepath,".bim"), mustWork = TRUE)
    bed_file <- normalizePath(paste0(filepath,".bed"), mustWork = TRUE)

    cat(".fam file: \n", fam_file)
    cat(".bim file: \n", bim_file)
    cat(".bed file: \n", bed_file)

    #columns : 1 - Family ID, 2 - ind ID, 3 - father, 4 - mother, 5 -sex, 6 - phenotype (1/2)
    fam_file <- fread(fam_file, sep="\t", colClasses = list(character=1:4), header=FALSE, data.table=FALSE)

    trios <- which((fam_file[,3] != "0") & (fam_file[,4] != "0"))
    check_fathers_present <- fam_file[trios,4] %in% fam_file[,1]
    check_mothers_present <- fam_file[trios,3] %in% fam_file[,1]

    if ( sum(! check_fathers_present) > 0){
        print("some father ID's not present in individual ID's: ")
        print(fam_file[! check_fathers_present,])
        stop()
    }
    if ( sum(! check_mothers_present) > 0){
        print("some mother ID's not present in individual ID's: ")
        print(fam_file[! check_mothers_present])
        stop()
    }

    full_ID <- paste0(fam_file[,1] , "_", fam_file[,2])
    trio_c  <- trios
    trio_f  <- which(full_ID %in% paste0(fam_file[trios,1],"_",fam_file[trios,3]) )
    trio_m  <- which(full_ID %in% paste0(fam_file[trios,1],"_",fam_file[trios,4]) )
    trio_pheno <- fam_file[trios,6]

    cat("number of individuals in .fam file: ", length(full_ID))
    cat("number of case trios: ", sum(trio_pheno == 2))
    cat("number of control trios: ", sum(trio_pheno == 1))

    snp_info <- fread(bim_file, drop = c(3,4), sep="\t")
    names(snp_info) <- c("chr", "snp", "A1", "A2")

    cat("number of SNP's: ", nrow(snp_info))
    
    n_case_trio  = 0L
    n_ctrl_trio  = 0L
    n_individual = 0L

    trio_counts <- .Call(count_genotypes, bed_file, trio_c, trio_f, trio_m, length(full_ID), nrow(snp_info), package="geno_count")
    colnames(trio_counts$case) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    colnames(trio_counts$ctrl) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    return(list(snp_info = snp_info, case_trios = trio_counts$case, ctrl_trios = trio_counts$ctrl))
}