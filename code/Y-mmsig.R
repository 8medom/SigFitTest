library(mmsig)          # https://github.com/evenrus/mmsig
library(dplyr)


cosmic3 <- read.table("../input/COSMIC_v3_SBS_GRCh38_for_mmsig-RELEVANT.txt", sep = ",", header = TRUE, check.names = FALSE)
mut_mat <- read.csv("data/data_for_MutationalPatterns.dat", sep = "\t", row.names = 1, header = TRUE)


set.seed(0)
sig_out <- mm_fit_signatures(muts.input=mut_mat,
                             sig.input=cosmic3,
                             input.format = "classes", # other option: classes
                             sample.sigt.profs = NULL,  # NULL = use all signatures provided in the reference. Optionally provide list with signatures to consider for each sample.
                             strandbias = FALSE,
                             bootstrap = FALSE,
                             refcheck=FALSE,    # check that input mutational catalog (if vcf-format) is aligned to hg19
                             cos_sim_threshold = 0.01,  # cosine similarity threshold below which signatures are removed from the final profile
                             dbg=FALSE)


df <- sig_out$estimate
df <- df[, !names(df) %in% c("mutations")]
res_to_save <- as.data.frame(df) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(t(res_to_save), file = "signature_results/mmsig-contribution.dat")
