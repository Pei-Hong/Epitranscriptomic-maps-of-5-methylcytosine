#!/bin/Rcscipt
args <- commandArgs()
# print(args)
inputFile = read.table(args[6], header = T, sep = "\t")

inputFile[, c("OR", "p", "ES")] = as.data.frame(t(as.data.frame(apply(inputFile, 1, function(x){
  if((as.numeric(x[2]) == 0 & as.numeric(x[3]) == 0) | (as.numeric(x[2]) == 0 & as.numeric(x[4]) == 0)){
    # print(as.numeric(x[12]))
    return(c(1,1,1))
  }else{
    act_in = as.numeric(x[2])
    act_out = as.numeric(x[3])
    n = act_in + act_out
    p = (act_in + as.numeric(x[4]))/(n + as.numeric(x[4]) + as.numeric(x[5]))
    t = prop.test(x = act_in, n = n, p = p)
    p.val = t$p.value
    OR = (act_in/act_out)/(as.numeric(x[4])/as.numeric(x[5]))
    ES = abs(2 * (asin(sqrt(act_in/n)) - asin(sqrt(p))))
    return(c(OR, p.val, ES))}
}))))
inputFile$p.adj = p.adjust(inputFile$p,method = "fdr")
outfile = paste0(strsplit(args[6], ".temp.txt")[[1]][1], ".txt", quote = "")
# print(outfile)
write.table(inputFile, file = outfile, col.names = T, row.names = F, quote = F, sep= "\t")

