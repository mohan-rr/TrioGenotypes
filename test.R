library(data.table)
fam_file <- fread("./g1k_CEU_filtered.fam", sep="\t", colClasses = list(character=1:4), header=FALSE, data.table=FALSE)
trios <- which((fam_file[,3] != "0") & (fam_file[,4] != "0"))
full_ID <- paste0(fam_file[,1] , "_", fam_file[,2])
trio_c  <- trios
trio_f  <- paste0(fam_file[trios,1], "_", fam_file[trios, 3])
trio_m  <- paste0(fam_file[trios,1], "_", fam_file[trios, 4])
trio_f  <- match(trio_f, full_ID)
trio_m  <- match(trio_m, full_ID)

trio_pheno <- fam_file[trios,6]
ped <- fread("./test.ped", select = 1:100)

genotype <- paste0(unlist(ped[,7]), unlist(ped[,8]))
paste(genotype[trio_m], genotype[trio_f], genotype[trio_c])
paste(trio_m, trio_f, trio_c)


library(TrioGenotypes)
trio_genotypes("./g1k_CEU_filtered")
