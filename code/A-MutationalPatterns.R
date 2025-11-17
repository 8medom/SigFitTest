library(MutationalPatterns)             # https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html
library(dplyr)

cosmic3 <- as.matrix(read.table("../input/COSMIC_v3.4_SBS_GRCh37_for_MutationalPatterns.txt", sep = ",", header = TRUE, check.names = FALSE))
mut_mat <- read.csv("data/data_for_MutationalPatterns.dat", sep = "\t", row.names = 1, header = TRUE)
res <- fit_to_signatures(mut_mat + 0.0001, cosmic3)
res_to_save <- as.data.frame(res$contribution) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_to_save, file = "signature_results/MutationalPatterns-contribution.dat")
