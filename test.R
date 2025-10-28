library(dplyr)
library(survival)
library(survminer)
library(preprocessCore)

source("/home/kyriakoc/Documents/code/originNet/src/classification.R")
# source all scripts in R/
r_files <- list.files(path = "../LGPclassifiers/R", pattern = "*.R", full.names = TRUE)
sapply(r_files, source)

rna.tpm <- readRDS("/home/kyriakoc/Documents/code/LGP_Build/build/rna.tpm.rds")
# remove TRANSCEND sample that results to no TPM values after gene collapsing
rna.tpm$TRANSCEND <- rna.tpm$TRANSCEND[, which(!colnames(rna.tpm$TRANSCEND) %in% c("017001-0021-0298_TRANSCEND"))]

rna.sf <- readRDS("/home/kyriakoc/Documents/code/LGP_Build/build/rna.sf.counts.rds")
rna.tpm.quant <- readRDS("/home/kyriakoc/Documents/code/LGP_Build/build/rna.tpm.hugo.quant.rds")
rna.pctRanks <- readRDS("/home/kyriakoc/Documents/code/LGP_Build/build/rna.pctRanks.rds")

rna.tpm.hugo <- lapply(rna.tpm, function(x) {
  collapseToGenes(x)
})
coo.reddy <- lapply(seq_along(rna.tpm), function(i) {
  print(names(rna.tpm)[i])
  computeCOO(query = rna.tpm[[i]], useReference = FALSE)
})
coo.reddy <- purrr::reduce(coo.reddy, dplyr::full_join) %>%
  dplyr::rename(reddy = ssREFERENCE, reddy_score = ssREFERENCE_score) %>%
  mutate(reddy = ordered(reddy, levels = c("GCB", "Unclassified", "ABC"))) %>%
  dplyr::rename(ID_Cohort = ID)

coo.ssREFERENCE <- lapply(seq_along(rna.tpm), function(i) {
  print(names(rna.tpm)[i])
  computeCOO(query = rna.tpm[[i]])
})
coo.ssREFERENCE <- purrr::reduce(coo.ssREFERENCE, dplyr::full_join) %>%
  mutate(ssREFERENCE = ordered(ssREFERENCE, levels = c("GCB", "Unclassified", "ABC"))) %>%
  dplyr::rename(ID_Cohort = ID)

rna.tpm.quant <- rna.tpm.quant[!names(rna.tpm.quant) %in% c("LENZ", "REMODLB")]
coo.originNet <- lapply(seq_along(rna.tpm.quant), function(i) {
  print(names(rna.tpm.quant)[i])
  computeOriginNet(rna.tpm.quant[[i]])
})
coo.originNet <- purrr::reduce(coo.originNet, dplyr::full_join) %>%
  dplyr::rename(ID_Cohort = ID)

coo.roche <- lapply(seq_along(rna.sf), function(i) {
  print(names(rna.sf)[i])
  computeRoche(rna.sf[[i]])
})
coo.roche <- purrr::reduce(coo.roche, dplyr::full_join) %>%
  dplyr::rename(ID_Cohort = ID)

coo.all <- coo.reddy %>%
  left_join(coo.ssREFERENCE) %>%
  left_join(coo.originNet) %>%
  left_join(coo.roche) %>%
  relocate(ID_Cohort) %>%
  as_tibble()

coo.all.new <- coo.reddy %>%
  left_join(coo.ssREFERENCE) %>%
  as_tibble()

coo.all <- readRDS(paste0("/home/kyriakoc/Documents/code/LGP_Build/build/coo.all.rds"))
# pass to coo.all new reddy and ssREFERENCE calls
coo.all <- coo.all %>%
  left_join(coo.all.new, by = "ID_Cohort", suffix = c(".old", ".new"))

# compare old vs new calls
table(coo.all$reddy.old, coo.all$reddy.new)
table(coo.all$ssREFERENCE.old, coo.all$ssREFERENCE.new)

saveRDS(coo.all, paste0(derived_path, "coo.all.rds"))
