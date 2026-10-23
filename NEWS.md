# LGPclassifiers 2.0.0

- Breaking: `dz.prob` and `dz.neg` conditional probabilities removed; only the marginal probability is returned
- originNet DZ trained on DLBCL90 scores
- DHit signature classifier added (`computeDHit`, `runDHitSig`), with `DHITstatus` returning 3 classes, and new `dhit.model` and `dhit_parameters` datasets
- New functions `getClassFromDLBCL90_DZ` and `getOriginNetConfidence`
- New datasets `probConfIntData` and `probConfIntData.dz`
- DZsig model revised (now 50 genes) and COO model retrained
- Public release: the TME26 classifier (`computeTME`, `tmeClassifier`, `modelM3.coefficients`) and the Roche COO model (`computeRoche`, `gneSeqCOO`, `ref.goya`) are not included in this version

# LGPclassifiers 1.0.1

- Bugfix for missing dependency

# LGPclassifiers 1.0.0

- First release including originNet single-sample classifier

# LGPclassifiers 0.1.0

- First release including TME26 classifier

# LGPclassifiers 0.0.0.9000

- First development version
