library(siglasso)       # https://github.com/gersteinlab/siglasso, https://github.gersteinlab.org/siglasso/
library(dplyr)


cosmic3 <- read.table("../input/COSMIC_v3_SBS_GRCh38_for_sigLASSO.txt", row.names = 1, header = TRUE, check.names = FALSE, sep = ",")
spectrum <- data.matrix(read.table("data/data_for_sigLASSO_spectrum.dat", row.names = 1, header = TRUE, check.names = FALSE))
res <- siglasso(spectrum, signature = data.matrix(cosmic3))
res_round <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_round, file = "signature_results/sigLASSO-contribution.dat")
