#' @title Collapse transcript_ids or gene_ids to gene name level
#' @description
#' Collapse transcript_ids or gene_ids to gene name level
#' @param tpm.mat matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param bm data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param featureType either "gene_id", "transcript_id" (default "gene_id")
#' @importFrom dplyr .data
#' @export
collapseToGenes <- function(tpm.mat, bm = NULL, featureType = "gene_id") {

  if (is.null(bm)) {
    stop("A data frame mapping Ensembl gene/transcript ids to gene names is required")
  }

  sums <- round(apply(tpm.mat, 2, sum), digits = 0)
  if (any(sums < 1000000 * 0.9999 | sums > 1000000 * 1.0001)) {
    stop("Input tpm.mat is not a TPM matrix! Check column sums (some do not add up to 1 million)")
  }
  rownames(tpm.mat) <- gsub(x = rownames(tpm.mat), pattern = "\\.\\d+", replacement = "")

  if (!featureType %in% c("gene_id", "transcript_id")) {
    stop('Not valid featureType ("gene_id" or "transcript_id")')
  }

  # remove unnecesary transcripts/gene_ids to speed up
  tpm.mat <- tpm.mat[rowSums(tpm.mat) > 0, ]

  # get the total expression of all transcripts/gene_ids for each gene
  tpm.mat <- as.data.frame(tpm.mat) %>%
    tibble::rownames_to_column("ID") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(gene_name = bm[match(.data$ID, bm[, featureType]), "gene_name"]) %>%
    dplyr::select(-.data$ID) %>%
    dplyr::group_by(.data$gene_name) %>%
    dplyr::summarise_all(mean, na.rm = TRUE) %>%
    dplyr::filter(!is.na(.data$gene_name)) %>%
    tibble::column_to_rownames("gene_name") %>%
    as.matrix()

  # Normalize again to reads per million
  tpm.mat <- apply(tpm.mat, 2, function(x) {
    10^6 * x / sum(x, na.rm = TRUE)
  })

  return(tpm.mat)
}

