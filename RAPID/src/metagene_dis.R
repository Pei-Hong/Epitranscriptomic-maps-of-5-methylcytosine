#!/bin/Rscript
library(ggplot2)

args <- commandArgs()
# print(args)
inputFile = read.table(args[6], header = T, sep = "\t")

inputFile = inputFile[, c(12:14)]
names(inputFile) = c("region", "len", "dis")

mean_3utr = mean(inputFile[inputFile$region == '3utr', "len"])
mean_5utr = mean(inputFile[inputFile$region == '5utr', "len"])
mean_cds = mean(inputFile[inputFile$region == 'cds', "len"])

nor_3utr = 3 * (mean_3utr)/(mean_3utr + mean_5utr + mean_cds)
nor_5utr = 3 * (mean_5utr)/(mean_3utr + mean_5utr + mean_cds)
nor_cds = 3 * (mean_cds)/(mean_3utr + mean_5utr + mean_cds)

inputFile$dis = as.numeric(as.character(inputFile$dis))
inputFile = na.omit(inputFile)
inputFile$dis_per = inputFile$dis / inputFile$len

# print(which(inputFile$dis_per >1))
inputFile$dis_nor = ifelse(inputFile$region == "3utr", nor_5utr + nor_cds + nor_3utr* inputFile$dis_per,
						ifelse(inputFile$region == "5utr", nor_5utr * inputFile$dis_per, nor_5utr + nor_cds * inputFile$dis_per) )
out_file = paste0(strsplit(args[6], "[.]dis")[[1]][1], ".pdf")
# print(args)
print(out_file)
pdf(file = out_file, height = 3, width = 4)
ggplot(inputFile, aes(dis_nor)) + geom_density(size =1.5, bw = 0.15, color = "#5B9BD5") + theme_classic() + geom_vline(xintercept = c(nor_5utr, c(nor_cds + nor_5utr)), linetype = 4, color = "grey") 
dev.off()



