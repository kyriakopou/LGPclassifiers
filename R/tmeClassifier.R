#' @title Apply TME26 classifier
#' @description
#' Normalize tpm sample counts based on housekeeping gene levels
#' @param exprs Input standardized expression data 
#' @export
#' 
tmeClassifier = function(exprs)
{
  # Read coefficients 
  # modelM3.coefficientsTable=data.frame(read.table(paste(models_path,paste="LinearModel_M3_Coefficients.txt",sep=""),sep="\t"))
  modelM3.coefficientsTable[,1] = as.character(modelM3.coefficientsTable[,1])
  modelM3.coefficientsTable[-1,1] = toupper(as.character(modelM3.coefficientsTable[-1,1]))
  
  # modelM3.threshold=read.table(paste(models_path,paste="LinearModel_M4_Threshold.txt",sep=""),sep="\t")
  modelM3.threshold=modelM3.threshold[1,1]
  modelM3.genes=as.character(modelM3.coefficientsTable[-1,1])
  myModel=modelM3.coefficientsTable[,2]
  names(myModel)=modelM3.coefficientsTable[,1]
  intercept <- myModel[1]
  
  myModel = myModel[which(names(myModel) %in% rownames(exprs))]
  scaledData = exprs
  scores.M3=apply(scaledData[toupper(names(myModel)),], 2, function(x,y){sum(x*y)}, myModel)+intercept
  calls.M3 <- scores.M3>modelM3.threshold
  data_pred <- data.frame(ID=colnames(scaledData),"RNAseq_M3_Scores"=scores.M3,"RNAseq_M3_Calls"=calls.M3)
  
  inter = length(intersect(modelM3.coefficientsTable[,1], rownames(exprs)))
  missing = modelM3.genes[-which(modelM3.genes %in% rownames(exprs))]
  if (inter<25)
    print(paste("Warning: ",25-inter," of 25 TME genes missing: ",missing,sep=""))
  
  return(data_pred)
}
