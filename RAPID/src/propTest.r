#!/bin/Rcscipt
args <- commandArgs()
# print(args)
inputFile = read.table(args[6], header = T, sep = "\t")

inputFile[, c("OR", "p", "ES")] = as.data.frame(t(as.data.frame(apply(inputFile, 1, function(x){
  if((as.numeric(x[12]) == 0 & as.numeric(x[13]) == 0) | (as.numeric(x[12]) == 0 & as.numeric(x[14]) == 0)){
    # print(as.numeric(x[12]))
    return(c(1,1,1))
  }else{
    act_in = as.numeric(x[12])
    act_out = as.numeric(x[13])
    n = act_in + act_out
    p = (act_in + as.numeric(x[14]))/(n + as.numeric(x[14]) + as.numeric(x[15]))
    t = prop.test(x = act_in, n = n, p = p)
    p.val = t$p.value
    OR = (act_in/act_out)/(as.numeric(x[14])/as.numeric(x[15]))
    ES = abs(2 * (asin(sqrt(act_in/n)) - asin(sqrt(p))))
    return(c(OR, p.val, ES))}
}))))
inputFile$p.adj = p.adjust(inputFile$p,method = "fdr")
inputFile$EnrichValue = log2(inputFile$OR) * (-1 * log10(inputFile$p.adj))
outfile = paste0(strsplit(args[6], ".temp.txt")[[1]][1], ".txt", quote = "")
# print(outfile)
write.table(inputFile, file = outfile, col.names = T, row.names = F, quote = F, sep= "\t")

