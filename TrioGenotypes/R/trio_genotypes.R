
trio_genotypes <- function(filepath, missing_as_ctrl = TRUE){
    fam_file <- normalizePath(paste0(filepath,".fam"), mustWork = TRUE)
    bim_file <- normalizePath(paste0(filepath,".bim"), mustWork = TRUE)
    bed_file <- normalizePath(paste0(filepath,".bed"), mustWork = TRUE)

    cat("\n fam file: \n", fam_file)
    cat("\n bim file: \n", bim_file)
    cat("\n bed file: \n", bed_file)
    cat("\n ------------- \n ")
    
    fam_file <- check_fam_file(fam_file, missing_as_ctrl)
    snp_info <- fread(bim_file, drop = c(3,4), sep="\t")
    names(snp_info) <- c("chr", "snp", "A1", "A2")
    cat("\n number of SNP's: ", nrow(snp_info))

    # trio_counts <- count_genotypes( toString(bed_file), trio_c, trio_f, trio_m, trio_pheno, length(full_ID), nrow(snp_info))
    trio_counts <- .Call(`_TrioGenotypes_count_genotypes`, toString(bed_file), 
                        fam_file$trio_c, fam_file$trio_f, fam_file$trio_m, fam_file$trio_pheno, 
                        fam_file$n_ind, nrow(snp_info))

    # res <- microbenchmark(count_genotypes(bed_file, trio_c, trio_f, trio_m, trio_pheno, length(full_ID), nrow(snp_info)))
    rownames(trio_counts$case) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    rownames(trio_counts$ctrl) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    return(invisible(list(snp_info = snp_info, case_trios = trio_counts$case, ctrl_trios = trio_counts$ctrl)))
}


check_fam_file <- function(fam_file, missing_as_ctrl) { 

    fam_file <- fread(fam_file, colClasses = list(character=1:4), header=FALSE, data.table=FALSE)
    #columns : 1 - Family ID, 2 - ind ID, 3 - father, 4 - mother, 5 -sex, 6 - phenotype (1/2)
    trio_c  <- which((fam_file[,3] != "0") & (fam_file[,4] != "0"))
    trio_f  <- paste0(fam_file[trio_c,1], "_", fam_file[trio_c, 3])
    trio_m  <- paste0(fam_file[trio_c,1], "_", fam_file[trio_c, 4])

    #incase there are duplicate ind. ID (w/ different family ID), so have to combine famID + IID
    full_ID <- paste0(fam_file[,1] , "_", fam_file[,2])
    missing <- which(!(trio_f %in% full_ID) | !(trio_m %in% full_ID))

    trio_c  <- trio_c[-missing]
    trio_f  <- match(trio_f, full_ID)[-missing]
    trio_m  <- match(trio_m, full_ID)[-missing]
    trio_pheno <- fam_file[trio_c,6]

    cat("\n trio children have the following phenotypes: \n")
    print(unique(trio_pheno))   
    cat("\n '1' is ctrl, '2' is case. Any other phenotype value considered missing data.")

    if (missing_as_ctrl){
        cat("\n treating missing phenotype trio as control")
        missing_pheno <- (trio_pheno != 1) & (trio_pheno != 2)
        cat("\n ", sum(missing_pheno), "missing phenotype trios treated as control")
        trio_pheno[missing_pheno] <- 1
    }

    cat("\n ------------- \n ")
    cat("\n number of individuals in .fam file: ", length(full_ID))
    cat("\n number of case trios: ", sum(trio_pheno == 2))
    cat("\n number of control trios: ", sum(trio_pheno == 1))

    return(list(trio_c = trio_c, trio_f = trio_f, trio_m = trio_m, trio_pheno = trio_pheno, n_ind = length(full_ID)))
}