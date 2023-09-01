library(dplyr)
library(sigfit)     # https://github.com/kgori/sigfit


cosmic3 <- read.table("../input/COSMIC_v3_SBS_GRCh38_for_sigfit-RELEVANT.txt", row.names = 1, header = TRUE, check.names = FALSE, sep = ",")
context_file <- read.table("data/data_for_sigfit.dat", header = TRUE, check.names = FALSE)
counts <- build_catalogues(context_file)
dim(counts)


fit <- fit_signatures(counts = counts, signatures = cosmic3, chains = 1, seed = 0)    # iter = 6000, warmup = 2000
exposures <- retrieve_pars(fit, par = "exposures", hpd_prob = 0.90)
x <- as.data.frame(exposures$mean) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(t(x), file = "signature_results/sigfit-contribution.dat")
y <- as.data.frame(exposures$lower_90) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(t(y), file = "signature_results/sigfit-contribution_lower90.dat")
