library(deconstructSigs)    # https://github.com/raerose01/deconstructSigs
library(coop)
library(dplyr)


cosmic3 <- read.table("../input/COSMIC_v3_SBS_GRCh38-RELEVANT.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
cosmic3 <- as.data.frame(t(cosmic3))
mut_matrix <- read.csv("data/data_for_deconstructSigs.dat", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
mut_matrix <- as.data.frame(t(mut_matrix)) # to get the same format as randomly.generated.tumors
identical(colnames(cosmic3), colnames(mut_matrix))


df <- data.frame()
for (sample in rownames(mut_matrix)){
  cat('Processing sample', sample, '\n')
  res = whichSignatures(tumor.ref = mut_matrix, sample.id = sample, signatures.ref = cosmic3,
                        contexts.needed = TRUE, signature.cutoff = 0)
  vals <- as.data.frame(c(res$weights[1,]))
  rownames(vals) <- sample
  df <- rbind(df, vals)
  }
df <- df %>% mutate_if(is.numeric, round, digits = 4)
write.csv(t(df), file = "signature_results/deconstructSigs-contribution.dat")
