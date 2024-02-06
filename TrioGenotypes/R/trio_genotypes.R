
trio_genotypes <- function(filepath){
    fam_file <- normalizePath(paste0(filepath,".fam"), mustWork = TRUE)
    bim_file <- normalizePath(paste0(filepath,".bim"), mustWork = TRUE)
    bed_file <- normalizePath(paste0(filepath,".bed"), mustWork = TRUE)

    cat("\n fam file: \n", fam_file)
    cat("\n bim file: \n", bim_file)
    cat("\n bed file: \n", bed_file)

    #columns : 1 - Family ID, 2 - ind ID, 3 - father, 4 - mother, 5 -sex, 6 - phenotype (1/2)
    fam_file <- fread(fam_file, sep="\t", colClasses = list(character=1:4), header=FALSE, data.table=FALSE)

    trios <- which((fam_file[,3] != "0") & (fam_file[,4] != "0"))
    check_fathers_present <- fam_file[trios,3] %in% fam_file[,2]
    check_mothers_present <- fam_file[trios,4] %in% fam_file[,2]
    #check the genealogy is complete... should always be when generated from PLINK, but nonetheless...
    if ( sum(! check_fathers_present) > 0){
        print(fam_file[! check_fathers_present,])
        stop("some father ID's not present in individual ID's")
    }
    if ( sum(! check_mothers_present) > 0){
        print(fam_file[! check_mothers_present])
        stop("some mother ID's not present in individual ID's")
    }

    trio_c  <- trios
    trio_f  <- paste0(fam_file[trios,1], "_", fam_file[trios, 3])
    trio_m  <- paste0(fam_file[trios,1], "_", fam_file[trios, 4])
    #incase there are duplicate ind. ID (w/ different family ID), so have to combine famID + IID
    full_ID <- paste0(fam_file[,1] , "_", fam_file[,2])
    trio_f  <- match(trio_f, full_ID)
    trio_m  <- match(trio_m, full_ID)
    trio_pheno <- fam_file[trios,6]

    cat("\n number of individuals in .fam file: ", length(full_ID))
    cat("\n number of case trios: ", sum(trio_pheno == 2))
    cat("\n number of control trios: ", sum(trio_pheno == 1))

    snp_info <- fread(bim_file, drop = c(3,4), sep="\t")
    names(snp_info) <- c("chr", "snp", "A1", "A2")
    cat("\n number of SNP's: ", nrow(snp_info))

    # trio_counts <- count_genotypes( toString(bed_file), trio_c, trio_f, trio_m, trio_pheno, length(full_ID), nrow(snp_info))
    trio_counts <- .Call(`_TrioGenotypes_count_genotypes`, toString(bed_file), trio_c, trio_f, trio_m, trio_pheno, length(full_ID), nrow(snp_info))
    # res <- microbenchmark(count_genotypes(bed_file, trio_c, trio_f, trio_m, trio_pheno, length(full_ID), nrow(snp_info)))
    rownames(trio_counts$case) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    rownames(trio_counts$ctrl) <- c("222", "212", "211", "122", "121", "201", "021", "112", "111", "110", "101", "100", "011", "010", "000")
    return(invisible(list(snp_info = snp_info, case_trios = trio_counts$case, ctrl_trios = trio_counts$ctrl)))
}
