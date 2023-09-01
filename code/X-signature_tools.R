library(signature.tools.lib)            # https://github.com/Nik-Zainal-Group/signature.tools.lib
library(dplyr)


mut_mat <- read.csv("data/data_for_MutationalPatterns.dat", sep = "\t", row.names = 1, header = TRUE)
cosmic3 <- as.matrix(read.table("../input/COSMIC_v3_SBS_GRCh38_for_MutationalPatterns.txt", sep = ",", header = TRUE, check.names = FALSE))
res <- Fit(catalogues = mut_mat, signatures = cosmic3, useBootstrap = TRUE, nboot = 10, nparallel = 1)


res_to_save <- t(res$exposures)
res_to_save <- res_to_save[!(row.names(res_to_save) == "unassigned"),]
res_to_save <- as.data.frame(res_to_save) %>% mutate_if(is.numeric, round, digits = 2)
write.csv(res_to_save, file = "signature_results/signature_tools-contribution.dat")
