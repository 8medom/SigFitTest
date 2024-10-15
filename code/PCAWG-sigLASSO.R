library(siglasso)       # https://github.com/gersteinlab/siglasso, https://github.gersteinlab.org/siglasso/
library(dplyr)


# prepare input data
cosmic3 <- read.table("../input/COSMIC_v3_SBS_GRCh38_for_sigLASSO.txt", row.names = 1, header = TRUE, check.names = FALSE, sep = ",")
spectrum <- data.matrix(read.table("../real mutational catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat", row.names = 1, header = TRUE, check.names = FALSE))
rownames(spectrum) <- sub("\\]", "", sub("\\[", "", rownames(spectrum)))
spectrum <- spectrum[rownames(cosmic3), ]


# fit signatures
res <- siglasso(spectrum, signature = data.matrix(cosmic3), plot = F)
res_round <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_round, file = "../WGS_PCAWG-146_samples-sigLASSO-contribution.dat")
