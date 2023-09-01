library(SigsPack)       # https://bioconductor.org/packages/release/bioc/html/SigsPack.html
library(dplyr)


cosmic3 <- as.matrix(read.table("../input/COSMIC_v3_SBS_GRCh38_for_MutationalPatterns.txt", sep = ",", header = TRUE, check.names = FALSE))
mut_mat <- as.matrix(read.csv("data/data_for_MutationalPatterns.dat", sep = "\t", row.names = 1, header = TRUE))


res <- signature_exposure(mut_mat, P = cosmic3)
res_to_save <- as.data.frame(res$exposures) %>% mutate_if(is.numeric, round, digits = 4)
old_names <- colnames(res_to_save)
new_names <- c()
for (i in 1:length(old_names)){
  new_names <- append(new_names, sprintf("S%d", i - 1))
  }
colnames(res_to_save) <- new_names
write.csv(res_to_save, file = "signature_results/SigsPack-contribution.dat")
