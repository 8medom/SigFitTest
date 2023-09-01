library(YAPSA)      # http://bioconductor.org/packages/3.12/bioc/vignettes/YAPSA/inst/doc/YAPSA.html
library(dplyr)


cosmic3 <- read.csv("../input/COSMIC_v3_SBS_GRCh38_for_YAPSA.txt", sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
mut_mat <- read.csv("data/data_for_YAPSA.dat", sep = "\t", row.names = 1, header = TRUE)


res <- LCD(mut_mat,cosmic3)
res_to_save <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_to_save, file = "signature_results/YAPSA-contribution.dat")
