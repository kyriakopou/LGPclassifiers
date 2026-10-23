#' Helper function to calculate gene mean
#' @noRd
gene.mean <- function(X, mat, bm, featureType) {
  if (featureType == "gene_id") {
    apply(mat[rownames(mat) %in% bm$gene_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, mean, na.rm = T)
  } else {
    apply(mat[rownames(mat) %in% bm$transcript_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, mean, na.rm = T)
  }
}

gene.sum <- function(X, mat, bm, featureType) {
  if (featureType == "gene_id") {
    apply(mat[rownames(mat) %in% bm$gene_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, sum, na.rm = T)
  } else {
    apply(mat[rownames(mat) %in% bm$transcript_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, sum, na.rm = T)
  }
}


#' Helper function to calculate geometric mean
#' @noRd
geometric.mean <- function(x) {
  v <- na.omit(x)
  exp(mean(log(v[v > 0])))
}

#' Helper function to calculate Hotelling's T squared
#' @noRd
computeTsquared <- function(pc, newdata, lim, lambda = NULL) {
  scores <- newdata %*% pc[, lim]
  if (is.null(lambda)) {
    lambda <- apply(scores, 2, var)  # Eigenvalues
  }
  T_squared <- rowSums((scores**2 %*% diag(1 / lambda)))
  return(list(T_squared = T_squared, lambda = lambda))
}

#' Helper function to project new data onto PCA model of reference dataset
#' @noRd
getPCAprojection <- function(newdata, pc, ncomp, eigenvalues) {
  T_squared <- computeTsquared(pc, newdata, ncomp, lambda = eigenvalues[ncomp])
  scores <- newdata %*% pc
  return(list(T2 = T_squared$T_squared, scores = scores))
}

#' Helper function to compute confidence probability of COO calls
#' @noRd
confidenceCOO <- function(query.quant, df) {
  probConfIntData <- LGPclassifiers::probConfIntData
  # load("data/probConfIntData.rda")



  # Project on PCA to get OOD
  new_data <- scale(t(query.quant)[, rownames(probConfIntData$pcaModelRobust$PCA)],
    center = probConfIntData$robust.mean[rownames(probConfIntData$pcaModelRobust$PCA)],
    scale = F)
  project <- getPCAprojection(newdata = new_data, pc = probConfIntData$pcaModelRobust$PCA,
    ncomp = probConfIntData$pcaModelRobust$ncomp,
    eigenvalues = probConfIntData$pcaModelRobust$eigenvalues)

  df <- df %>%
    mutate(T2.on = project$T2) %>%
    # scale wrt ROBUST T2 min, max values
    mutate(T2.on_scaled = scale_values(T2.on,
      min = min(probConfIntData$pcaModelRobust$T2),
      max = max(probConfIntData$pcaModelRobust$T2))) %>%
    # compute originNet score distances (default values for min/max)
    mutate(originNet_score_scaled = scale_values((originNet_score)))

  df$coo_dist_scaled <- sapply(df$originNet_score_scaled, function(x) {
    a <- scale_values(1907.8)
    b <- scale_values(2433.5)
    min(abs(x - a), abs(x - b))
  })

  df <- df %>%
    # mutate(logit_odds.on = probConfIntData$probConfCoef[1] + probConfIntData$probConfCoef[2] * T2.on_scaled + probConfIntData$probConfCoef[3] * originNet_score_scaled) %>%
    mutate(logit_odds.on = probConfIntData$probConfCoef[1] + probConfIntData$probConfCoef[2] * originNet_score_scaled) %>%
    mutate(prob.ABC =  sigmoid(logit_odds.on)) %>%
    mutate(prob.GCB = 1 - prob.ABC) %>%
    mutate(conf.coo = pmax(prob.GCB, prob.ABC))

  return(df)
}

#' Helper function to compute confidence probability of DZSig calls
#' @noRd
confidenceDZ <- function(query.quant, df) {
  probConfIntData.dz <- LGPclassifiers::probConfIntData.dz
  # load("data/probConfIntData.dz.rda")

  # Project on PCA to get OOD
  new_data <- scale(t(query.quant)[, rownames(probConfIntData.dz$pcaModelRobust$PCA)],
    center = probConfIntData.dz$robust.mean[rownames(probConfIntData.dz$pcaModelRobust$PCA)],
    scale = F)
  project <- getPCAprojection(newdata = new_data, pc = probConfIntData.dz$pcaModelRobust$PCA,
    ncomp = probConfIntData.dz$pcaModelRobust$ncomp,
    eigenvalues = probConfIntData.dz$pcaModelRobust$eigenvalues)

  df <- df %>%
    mutate(T2.dz = project$T2) %>%
    # scale wrt ROBUST T2 min, max values
    mutate(T2.dz_scaled = scale_values(T2.dz,
      min = min(probConfIntData.dz$pcaModelRobust$T2),
      max = max(probConfIntData.dz$pcaModelRobust$T2))) %>%
    # compute originNet score distances (default values for min/max)
    mutate(DZ_score_scaled = scale_values(DZSig_score, min = -45, max = 25))

  df <- df %>%
    # mutate(logit_odds.dz = probConfIntData.dz$probConfCoef[1] + probConfIntData.dz$probConfCoef[2] * T2.dz_scaled + probConfIntData.dz$probConfCoef[3] * DZ_score_scaled) %>%
    mutate(logit_odds.dz = probConfIntData.dz$probConfCoef[1] + probConfIntData.dz$probConfCoef[2] * DZ_score_scaled) %>%
    mutate(prob.dz.pos =  sigmoid(logit_odds.dz)) %>%
    mutate(prob.dz.neg = 1 - prob.dz.pos) %>%
    mutate(conf.GCB.dz = pmax(prob.dz.pos, prob.dz.neg))

  return(df)
}

#' Helper function to scale a vector of values between 0-1
#' @noRd
scale_values <- function(x, min = -1181.849, max = 4326.485) {
  x <- ifelse(x < min, min, ifelse(x > max, max, x))
  (x - min) / (max - min)
}

#' Helper function to compute sigmoid distribution for probabilities
#' @noRd
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}
