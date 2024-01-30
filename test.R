# locate the project dir
path <- system("find ~/Documents -name LGP_Build", intern = TRUE)
setwd(paste0(path))
rce_path <- paste(readLines("rce_path.txt"), collapse = " ")
data_path <- paste0(rce_path, "/raw_data/")
derived_path <- paste0(path, "/derived/")
output_path <- paste0(path, "/output/")

library(dplyr)
library(survival)
library(survminer)
library(preprocessCore)

makeRPM <- function(query) {
  # Normalize to reads per million
  query.tpm <- apply(query, 2, function(x) {
    1000000 * x / sum(x, na.rm = TRUE)
  })
}
id2geneName <- readRDS("/stash/results/dev/ugidosgm/DLBCL/id2geneName.rds")
id2geneName <- id2geneName[order(id2geneName$gene_id), ]

# GET THE SAME RESULTS USING THE PKG
rna_counts <- readRDS("../LGP_Build/derived/rna.counts.rds")
goya.genes.tpm <- rna_counts$GOYA

# execute stepwise
goya.tpm <- collapseToGenes(goya.genes.tpm, id2geneName = LGPclassifiers::geneName.map, featureType = "gene_id")
goya.scaled <- scaleTPM(goya.tpm, ref.mean = LGPclassifiers::robust.ref.mean, ref.sd = LGPclassifiers::robust.ref.sd)
goya.coo <- runReddyCOO(goya.scaled)
# run the wrapper computeCOO
goya.coo2 <- computeCOO(query = goya.genes.tpm, useReference = TRUE)
# all equal
table(goya.coo$ssREFERENCE, goya.coo2$ssREFERENCE)

meta_merged <- readRDS("../LGP_Build/derived/meta.merged.rds")
meta_merged %>%
  filter(Cohort == "GOYA") %>%
  select(ID, Nanostring_COO) %>%
  full_join(goya.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# REMODLB
remodlb.logtpm  <- rna_counts$REMODLB
remodlb.tpm <- makeRPM(remodlb.logtpm)
remodlb.coo <- computeCOO(query = remodlb.tpm, featureType = "gene", useReference = FALSE)

meta_merged %>%
  filter(Cohort == "REMODLB") %>%
  select(ID, Nanostring_COO) %>%
  full_join(remodlb.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# NDMER from unscalling Andy's scaled entry
# This does not work
# ndmer.tpm  <- (all_cohorts$NDMER$scaled + LGPclassifiers::robust.ref.mean) * LGPclassifiers::robust.ref.sd
# ndmer.tpm <- makeRPM(ndmer.tpm)
# ndmer.coo <- computeCOO(query = ndmer.tpm, featureType = "gene")

meta_merged %>%
  filter(Cohort == "MER", Stage == "Baseline") %>%
  select(ID, Nanostring_COO) %>%
  full_join(ndmer.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# NDMER from doing TPM in gene counts
ndmer.counts  <- rna_counts$NDMER
# Perform TPM normalization
dge <- edgeR::DGEList(ndmer.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
ndmer.tpm <- edgeR::cpm(dge.tpm)
ndmer.tpm <- makeRPM(ndmer.tpm)
ndmer.coo <- computeCOO(query = ndmer.tpm, featureType = "gene_id")
table(ndmer.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "MER", Stage == "Baseline") %>%
  full_join(ndmer.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# rrmer from doing TPM in gene counts
rrmer.counts  <- rna_counts$RRMER
# Perform TPM normalization
dge <- edgeR::DGEList(rrmer.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
rrmer.tpm <- edgeR::cpm(dge.tpm)
rrmer.tpm <- makeRPM(rrmer.tpm)
rrmer.coo <- computeCOO(query = rrmer.tpm, featureType = "gene_id")
table(rrmer.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "MER", Stage == "Relapse") %>%
  full_join(rrmer.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# ROBUST from doing TPM in gene counts
robust.counts  <- rna_counts$ROBUST
# Perform TPM normalization
dge <- edgeR::DGEList(robust.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
robust.tpm <- edgeR::cpm(dge.tpm)
robust.tpm <- makeRPM(robust.tpm)
robust.coo <- computeCOO(query = robust.tpm, featureType = "gene_id")
table(robust.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "ROBUST") %>%
  full_join(robust.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()


# COMMERCIAL from doing TPM in gene counts
commercial.counts  <- rna_counts$COMMERCIAL
# Perform TPM normalization
dge <- edgeR::DGEList(commercial.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
commercial.tpm <- edgeR::cpm(dge.tpm)
commercial.tpm <- makeRPM(commercial.tpm)
commercial.coo <- computeCOO(query = commercial.tpm, featureType = "gene_id")
table(commercial.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "COMMERCIAL") %>%
  full_join(commercial.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# TRANSCEND from doing TPM in gene counts
transcend.counts  <- rna_counts$TRANSCEND
# Perform TPM normalization
dge <- edgeR::DGEList(transcend.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
transcend.tpm <- edgeR::cpm(dge.tpm)
transcend.tpm <- makeRPM(transcend.tpm)
transcend.coo <- computeCOO(query = transcend.tpm, featureType = "gene_id")
table(transcend.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "TRANSCEND") %>%
  full_join(transcend.coo, by = "ID") %>%
  select(Hans_COO, ssREFERENCE) %>%
  table()

# CC99282 from doing TPM in gene counts
cc282.counts  <- rna_counts$CC99282
# Perform TPM normalization
dge <- edgeR::DGEList(cc282.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
cc282.tpm <- edgeR::cpm(dge.tpm)
cc282.tpm <- makeRPM(cc282.tpm)
cc282.coo <- computeCOO(query = cc282.tpm, featureType = "gene_id")
table(cc282.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "CC99282") %>%
  full_join(cc282.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

# CC122 from doing TPM in gene counts
cc122.counts  <- rna_counts$CC122
# Perform TPM normalization
dge <- edgeR::DGEList(cc122.counts)
dge.tpm <- edgeR::calcNormFactors(dge, method = "TMM")
cc122.tpm <- edgeR::cpm(dge.tpm)
cc122.tpm <- makeRPM(cc122.tpm)
cc122.coo <- computeCOO(query = cc122.tpm, featureType = "gene_id")
table(cc122.coo$ssREFERENCE)

meta_merged %>%
  filter(Cohort == "CC122") %>%
  full_join(cc122.coo, by = "ID") %>%
  select(Nanostring_COO, ssREFERENCE) %>%
  table()

robust.tpm <- readRDS(paste0(data_path, "ndData/ROBUST/robustTPMs_all1098.RDS"))
tpm.samples <- colnames(robust.tpm)

robust.samples <- meta_merged %>%
  filter(Cohort == "ROBUST") %>%
  select(ID) %>%
  unlist()

tpm.samples[tpm.samples %in% c("X7061001", "X706.1001")]
colnames(robust.tpm)[!colnames(robust.tpm) %in% robust.samples]
grep("X7061001", colnames(robust.tpm))

robust.manifest <- read.table(paste0(data_path, "ndData/ROBUST/robust_clinical_screening_v2.txt"), sep = "\t", header = TRUE)

mer_meta <- meta_merged %>%
  filter(Cohort == "MER" | Cohort == "nonMER")

saveRDS(mer_meta, paste0(derived_path, "mer.meta.rds"))
