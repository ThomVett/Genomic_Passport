  # Script to fin HLA Alleles for a given patient, european ancestry classifier
library(HIBAG)

args = commandArgs(trailingOnly=TRUE)

patientName=args[1]
folder = args[2]

if(identical(match("model.list",ls()),NA_integer_)) model.list <- get(load("European-HLA4-hg19.RData"))

#patientName <- "10437_chr6"
#setwd(folder)

patient <- c(paste(patientName,".bed",sep=""),
               paste(patientName,".fam",sep=""),
               paste(patientName,".bim",sep=""))
  
yourgeno <- hlaBED2Geno(bed.fn= patient[1], fam.fn=patient[2],
                          bim.fn=patient[3],assembly="hg19")
  
HLA_Loc = c("A","B","C","DRB1","DQA1","DQB1","DPB1")
  
Result = matrix(NA,nrow=4,ncol=7)
for (allele in 1:7){
  hla.id <- HLA_Loc[allele]
  model <- hlaModelFromObj(model.list[[hla.id]])
  
  cbind(frequency = model$hla.freq)
  pred.guess <- predict(model, yourgeno, type="response+prob",match.type="Position")
  Locus = c(pred.guess$locus,pred.guess$value$allele1[1],
            pred.guess$value$allele2[1],pred.guess$value$pro[1])
  Result[,allele] = Locus
}
resultsSTR = paste(patientName,"_results.csv",sep="")
write.csv(Result,file=resultsSTR)