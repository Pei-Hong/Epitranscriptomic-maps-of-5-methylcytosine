#!/bin/Rcscipt
args <- commandArgs()
# print(args)
inputFile = read.table(args[6], header = T, sep = "\t")

# inputFile[, c("OR", "p")] = as.data.frame(t(as.data.frame(apply(inputFile, 1, function(x){
# 		t = fisher.test(matrix(as.numeric(x[12:15]), ncol = 2))
# 		if((x[12] == 0 & x[13] == 0) | (x[12] == 0 & x[14] == 0)){
# 		  return(c(1,1))
# 		}else{
# 		  p = t$p.value
# 		  OR = t$estimate
# 		  return(c(OR, p))
# 		}
# 	}))))
inputFile[, c("OR", "p")] = as.data.frame(t(as.data.frame(apply(inputFile, 1, function(x){
  if((as.numeric(x[12]) == 0 & as.numeric(x[13]) == 0) | (as.numeric(x[12]) == 0 & as.numeric(x[14]) == 0)){
    # print(as.numeric(x[12]))
    return(c(1,1))
  }else{
    t = fisher.test(matrix(as.numeric(x[12:15]), ncol = 2))
    p = t$p.value
    OR = t$estimate
    return(c(OR, p))}
}))))
inputFile$p.adj = p.adjust(inputFile$p,method = "fdr")
outfile = paste0(strsplit(args[6], ".temp.txt")[[1]][1], ".txt", quote = "")
# print(outfile)
write.table(inputFile, file = outfile, col.names = T, row.names = F, quote = F, sep= "\t")

