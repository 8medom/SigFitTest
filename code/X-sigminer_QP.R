library(sigminer)   # https://github.com/ShixiangWang/sigminer
library(dplyr)


cosmic3 <- as.matrix(read.table("../input/COSMIC_v3_SBS_GRCh38_for_MutationalPatterns.txt", sep = ",", header = TRUE, check.names = FALSE))
mut_mat <- as.matrix(read.csv("data/data_for_MutationalPatterns.dat", sep = "\t", row.names = 1, header = TRUE))
res <- sig_fit(mut_mat, cosmic3, method = "QP")
res_to_save <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_to_save, file = "signature_results/sigminer_QP-contribution.dat")
