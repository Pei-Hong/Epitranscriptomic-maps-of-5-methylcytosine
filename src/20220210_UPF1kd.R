library(ggplot2)
library(dplyr)
library(DESeq2)
library(pheatmap)
library(reshape2)
library(VennDiagram)
library(easyGgplot2)

######## DEG analysis ###########
## read table
data.file = read.table("/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/50_NSUN2kd_UPF1kd/21_featurecount/combine_trim.counts",
                       header = T, stringsAsFactors = F, sep = "\t")
head(data.file)
data.matrix = data.file[, c("X9.LDH9543_L2.bam", "X10.LDH9544_L2.bam",  "X11.LDH9545_L2.bam",  "X12.LDH9546_L2.bam",  
                         "A.LDJ8187_L2.bam",  "B.LDJ8188_L2.bam",  "C.LDJ8189_L2.bam",  "D.LDJ8190_L2.bam",
                         "E.LDJ8191_L3.bam",  "F.LDJ8192_L3.bam",  "G.LDJ8193_L3.bam",  "H.LDJ8194_L3.bam")]
names(data.matrix)
treatment = rep(c("wt_scr", "wt_si", "ko_scr", "ko_si"), 3)
tmp = data.frame(file= c("9-LDH9543_L2.bam", "10-LDH9544_L2.bam",  "11-LDH9545_L2.bam",  "12-LDH9546_L2.bam",  
                         "A-LDJ8187_L2.bam",  "B-LDJ8188_L2.bam",  "C-LDJ8189_L2.bam",  "D-LDJ8190_L2.bam",
                         "E-LDJ8191_L3.bam",  "F-LDJ8192_L3.bam",  "G-LDJ8193_L3.bam",  "H-LDJ8194_L3.bam"),
                 sample= rep(c("N2wt_UPF1scr", "N2wt_UPF1si", "N2ko_UPF1scr", "N2ko_UPF1si"), 3))
rownames(data.matrix) = data.file$Geneid
annotation_col = data.frame(treat = factor(treatment, levels = c("wt_scr", "wt_si", "ko_scr", "ko_si")))
rownames(annotation_col) = names(data.matrix)

## bulid DEseq object and get size factors
sample1 = round(data.matrix)
rownames(sample1) = rownames(data.matrix) 
sample1 = sample1[rowSums(sample1[, c(1,5,9)] >=50) >=1 & 
                    rowSums(sample1[, c(2,6,10)] >=50) >=1 &
                    rowSums(sample1[, c(3,7,11)] >=50) >=1 &
                    rowSums(sample1[, c(4,8,12)] >=50) >=1, ]
condition.1 <- factor(treatment, levels = c("wt_scr", "wt_si", "ko_scr", "ko_si"))
coldata.1 <- data.frame(row.names = colnames(sample1), condition.1)
dds.1 <- DESeqDataSetFromMatrix(countData=sample1, colData=coldata.1, design=~condition.1)
dds.sizefactor <- estimateSizeFactors(dds.1)
sizeFactors <- sizeFactors(dds.sizefactor)
sizeFactors(dds.1) <- sizeFactors

sample1$gene = rownames(sample1)
deg_marker = c("COLGALT1", "FKBP10", "TMED9", "PEA15", "ITM2B")
HK_genes = c(deg_marker, "NSUN2", "UPF1", "GAPDH", "ACTB")
sample1$label = factor(ifelse(sample1$gene %in% HK_genes, sample1$gene, ""),
                       levels = c(HK_genes, ""))
sample1 = rbind(sample1[sample1$label == "", ], sample1[sample1$label != "", ])
sample1 %>%
    ggplot(aes(X9.LDH9543_L2.bam, X10.LDH9544_L2.bam, color = label)) + geom_point(alpha = 1, size = 0.8) + theme_classic() + 
      scale_x_log10() + scale_y_log10() + geom_abline(intercept = 0, slope = 1) + geom_text(aes(label = label)) +
  scale_color_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu"))(length(HK_genes)), "grey"))
sample1 %>% mutate(X9.LDH9543_L2.bam = X9.LDH9543_L2.bam/sizeFactors(dds.1)["X9.LDH9543_L2.bam"], 
                   X10.LDH9544_L2.bam = X10.LDH9544_L2.bam/sizeFactors(dds.1)["X10.LDH9544_L2.bam"]) %>%
    ggplot(aes(X9.LDH9543_L2.bam, X10.LDH9544_L2.bam, color = label)) + geom_point(alpha = 0.9, size = 0.8) + theme_classic() + 
      scale_x_log10() + scale_y_log10() + geom_abline(intercept = 0, slope = 1)  + #geom_text(aes(label = label)) + 
  scale_color_manual(values = c(colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu"))(length(HK_genes)), "grey"))
## check for replices by size factor ##
pheatmap(log(t(t(sample1[, 1:12])/sizeFactors(dds.1)+1)), scale = "row", show_rownames = F, 
         annotation_col = annotation_col)
test.marker = as.data.frame(t(t(sample1[c("NSUN2", "UPF1"), 1:12])/sizeFactors(dds.1)))
test.marker= test.marker %>% mutate(gene = rownames(test.marker)) %>% 
  melt(id = "gene") %>% mutate(treat = rep(condition.1, each = 2), rep = rep(c(1,2,3), each = 8))
ggplot(test.marker, aes(rep, value, fill = treat)) + 
  geom_bar(stat = "identity", position = "dodge") + theme_classic() + facet_grid(gene ~ .)
sample1.cor = cor(log(t(t(sample1[, 1:12])/sizeFactors(dds.1)+1)))
pheatmap(sample1.cor, show_colnames = F, annotation_col = annotation_col)

## estimate DEG 
dds.1 <- estimateDispersions(dds.1)
dds.1 <- nbinomWaldTest(dds.1)
res.wt <- results(dds.1, contrast = c('condition.1', 'wt_si', 'wt_scr'))
res.wt.data = as.data.frame(res.wt)
res.ko <- results(dds.1, contrast = c('condition.1', 'ko_si', 'ko_scr'))
res.ko.data = as.data.frame(res.ko)

res.N2ko <- results(dds.1, contrast = c('condition.1', 'ko_scr', 'wt_scr')) 
res.N2ko.data = as.data.frame(res.N2ko)
dim(as.data.frame(results(dds.1, contrast = c('condition.1', 'wt_si', 'wt_scr'))))
res.N2koUPF1si <- results(dds.1, contrast = c('condition.1', 'ko_si', 'wt_si'))
res.N2koUPF1si.data = as.data.frame(res.N2koUPF1si)

###############

####### Association analysis of UPF1 and NSUN2 #############
# res.wt.data     UPF1 ko/wt in NSUN2 wt cell line
res.wt.data$label = factor(ifelse(res.wt.data$log2FoldChange >= log2(1.2) & res.wt.data$padj <= 0.05, "up",
                                  ifelse(res.wt.data$log2FoldChange <= -log2(1.2) & res.wt.data$padj <= 0.05, "down", "unchange")),
                           levels = c("up", "down", "unchange"))
ggplot(res.wt.data, aes(log2(baseMean+1), log2FoldChange, color = label)) + geom_point(alpha = 0.4, size = 0.5) + theme_classic() +
  scale_color_manual(values = c("red", "blue", "grey")) + #geom_hline(yintercept = c(-log2(1.2), log2(1.2)))
  geom_hline(yintercept = c(0)) + ylim(c(-5,5)) + facet_grid()
# res.ko.data     UPF1 ko/wt in NSUN2 ko cell line
res.ko.data$label = factor(ifelse(res.ko.data$log2FoldChange >= log2(1.2) & res.ko.data$padj <= 0.05, "up",
                                  ifelse(res.ko.data$log2FoldChange <= -log2(1.2) & res.ko.data$padj <= 0.05, "down", "unchange")),
                           levels = c("up", "down", "unchange"))
ggplot(res.ko.data, aes(log2(baseMean+1), log2FoldChange, color = label)) + geom_point(alpha = 0.4, size = 0.5) + theme_classic() +
  scale_color_manual(values = c("red", "blue", "grey")) + #geom_hline(yintercept = c(-log2(1.2), log2(1.2)))
  geom_hline(yintercept = c(0)) + ylim(c(-5,5))
# res.N2ko.data     NSUN2 ko/wt in UPF1 wt cell line
res.N2ko.data$label = factor(ifelse(res.N2ko.data$log2FoldChange >= log2(1.2) & res.N2ko.data$padj <= 0.05, "up",
                                    ifelse(res.N2ko.data$log2FoldChange <= -log2(1.2) & res.N2ko.data$padj <= 0.05, "down", "unchange")),
                             levels = c("up", "down", "unchange"))
table(res.N2ko.data$label)
ggplot(res.N2ko.data, aes(log2(baseMean+1), log2FoldChange, color = label)) + geom_point(alpha = 0.4, size = 0.5) + theme_classic() +
  scale_color_manual(values = c("red", "blue", "grey")) + #geom_hline(yintercept = c(-log2(1.2), log2(1.2)))
  geom_hline(yintercept = c(0)) + ylim(c(-5,5))
# res.N2koUPF1si.data
res.N2koUPF1si.data$label = factor(ifelse(res.N2koUPF1si.data$log2FoldChange >= log2(1.2) & res.N2koUPF1si.data$padj <= 0.05, "up",
                                          ifelse(res.N2koUPF1si.data$log2FoldChange <= -log2(1.2) & res.N2koUPF1si.data$padj <= 0.05, "down", "unchange")),
                                   levels = c("up", "down", "unchange"))
table(res.N2koUPF1si.data$label)
ggplot(res.N2koUPF1si.data, aes(log2(baseMean+1), log2FoldChange, color = label)) + geom_point(alpha = 0.4, size = 0.5) + theme_classic() +
  scale_color_manual(values = c("red", "blue", "grey")) + #geom_hline(yintercept = c(-log2(1.2), log2(1.2)))
  geom_hline(yintercept = c(0)) + ylim(c(-5,5))

##### res.combine  ####
res.combine = cbind(res.wt.data[, c("log2FoldChange", "padj", "label")], res.ko.data[, c("log2FoldChange", "padj", "label")])
names(res.combine) = c("log2FoldChange.wt", "padj.wt", "label.wt", "log2FoldChange.ko", "padj.ko", "label.ko")
ggplot(res.combine, aes(log2FoldChange.wt, log2FoldChange.ko, color = label.wt)) + geom_point(size = 0.5, alpha = 0.5) + theme_classic() +
  geom_abline(intercept = 0, slope = 1, color = "black") + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_hline(yintercept = c(log2(1.2), -log2(1.2)), linetype = "dotdash") + 
  geom_vline(xintercept = c(log2(1.2), -log2(1.2)), linetype = "dotdash") + #geom_smooth() + #geom_density2d() +
  ylim(-5,5) + xlim(-5,5) +
  scale_color_manual(values = c("red", "blue", "grey"))

res.N2ko.data[c("TRIB3", "FADS2", "FKBP10", "NIBAN2", "OSMR", "SLC39A13", "FXR2", "NIBAN2", "ATG13"), ]
data.matrix[c("TRIB3", "FADS2", "FKBP10", "NIBAN2", "OSMR"), ]
##### m5C  #######
m5C.file = read.table("/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/1_20RC3C_30RC5C/3_RBP/1_HeLa/HeLa5.mRNA",
                      header = F, stringsAsFactors = F, sep = "\t")
# m5C.file = m5C.file[, c( "V1","V2","V3","V4", "V5", "V6","V8","V9")]
names(m5C.file) = c("chr", "start", "end", "sites", "value", "strand", "region", "gene")
m5C.file.mRNA = unique(m5C.file$gene) #1592 #1669

res.combine$if.m5C = "N"
res.combine[intersect(rownames(res.combine), m5C.file.mRNA), "if.m5C"] = "Y"
res.combine[is.na(res.combine)] 

## N2 dependent m5C #######
NSUN2_m5C.file = read.table("/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/Huang_Yang.N2depend.mRNA",
                            header = F, stringsAsFactors = F)
# NSUN2_m5C.file = NSUN2_m5C.file[, c( "V1","V2","V3","V4", "V5", "V6","V8","V9")]
names(NSUN2_m5C.file) = c("chr", "start", "end", "sites", "value", "strand", "region", "gene", "ID")
length(unique(NSUN2_m5C.file$gene))
dim(NSUN2_m5C.file) ## 1273
res.combine$if.N2m5C = "N"
res.combine[intersect(rownames(res.combine), unique(NSUN2_m5C.file$gene)), "if.N2m5C"] = "Y"
res.combine[is.na(res.combine)]

## UPF1 binding ########
UPF1.file3 = read.table("/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/bsRNA/70_protein/RBP3/20_sample/UPF1_binding/UPF1_more1_fpkm.txt",
                        header = T, stringsAsFactors = F, sep = "\t")
res.combine$UPF1.HeLa2 = 0
res.combine[UPF1.file3$Geneid, "UPF1.HeLa2"] = UPF1.file3$enrich
all.equal(res.combine[UPF1.file3$Geneid, "UPF1.HeLa2"], UPF1.file3$enrich)
length(intersect(UPF1.file3$Geneid, rownames(res.combine)))
res.combine[is.na(res.combine)]
res.combine <- na.omit(res.combine)
res.combine$UPF1.HeLa2.trim = ifelse(res.combine$UPF1.HeLa2 <=8, res.combine$UPF1.HeLa2, 8)
# res.combine$UPF1.HeLa2.trim = ifelse(res.combine$UPF1.HeLa2 < 4, 0,
#                                     ifelse(res.combine$UPF1.HeLa2 < 8, res.combine$UPF1.HeLa2, 8))
# res.combine = rbind(res.combine[res.combine$UPF1.HeLa2.trim <2, ],
#                     res.combine[res.combine$UPF1.HeLa2.trim >=2, ])
res.combine$UPF1.type = factor(ifelse(res.combine$UPF1.HeLa2 >=4, 4,
                                      ifelse(res.combine$UPF1.HeLa2 >=2, 2,
                                             ifelse(res.combine$UPF1.HeLa2 >1, 1, 0))), levels = c(0,1,2,4))
res.combine$UPF1.type2 = factor(ifelse(res.combine$UPF1.HeLa2 >=4, 4,
                                      ifelse(res.combine$UPF1.HeLa2 >=2, 2, 0)), levels = c(0,2,4))
table(res.combine$UPF1.type)
table(res.combine$UPF1.type2)
res.combine[is.na(res.combine)]
a = res.combine[res.combine$UPF1.type == 2, ]

## UPF1 overlapped NSUN2 dependenr m5C #######
UPF1_overlap.file.N2 = read.table("/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/UPF1_HeLa_more1_fpkm.overlapm5C.N2.txt",
                                  header = F, stringsAsFactors = F, sep = "\t")

dim(UPF1_overlap.file.N2) ## 104
res.combine$ifoverlap.N2 = "N"
length(intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)))

res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)), "ifoverlap.N2"] = "Y"
res.combine[is.na(res.combine)] 
table(res.combine$ifoverlap.N2)

res.combine$gene = rownames(res.combine)

median(res.combine[res.combine$UPF1.type == 0, "log2FoldChange.wt"])
median(res.combine[res.combine$UPF1.type == 1, "log2FoldChange.wt"])
median(res.combine[res.combine$UPF1.type == 2, "log2FoldChange.wt"])
median(res.combine[res.combine$UPF1.type == 4, "log2FoldChange.wt"])

median(res.combine[res.combine$UPF1.type == 0, "log2FoldChange.ko"])
median(res.combine[res.combine$UPF1.type == 1, "log2FoldChange.ko"])
median(res.combine[res.combine$UPF1.type == 2, "log2FoldChange.ko"])
median(res.combine[res.combine$UPF1.type == 4, "log2FoldChange.ko"])

t.test(res.combine[res.combine$UPF1.type == 0, "log2FoldChange.wt"],
       res.combine[res.combine$UPF1.type == 2, "log2FoldChange.wt"])

##### ########## boxplot for Log2 FoldChange of UPF1 kd ####### #######
# 1
# dim(UPF1_overlap.file2) # 75
length(intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)))
res.combine[res.combine$ifoverlap.N2 == "Y",
            c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type")] %>% 
  mutate(gene = rownames(res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),])) %>% 
  # select(gene, log2FoldChange.wt, log2FoldChange.ko, UPF1.type) %>%
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ variable) +
  scale_color_manual(values = c("black","green", "blue"))
res.combine.sub2 = res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),]
table(res.combine.sub2$ifoverlap.N2)

# 1 m5C not bound UPF1; UPF1 with 2 enrich not around; UPF1 with 4 enrich not around
res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V8)), rownames(res.combine)),
            c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type", "gene")] %>% filter(UPF1.type %in% c(1,2,4)) %>%
  # select(gene, log2FoldChange.wt, log2FoldChange.ko, UPF1.type) %>%
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ variable) +
  scale_color_manual(values = c("grey", "green", "blue"))
res.combine.sub2.1 = res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V8)), rownames(res.combine)),] %>%
  filter(UPF1.type %in% c(1,2,4))
# dim(res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V9)), rownames(res.combine)),])
# dim(res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file$V9)), rownames(res.combine)),])
# median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type2 == 0, "log2FoldChange.wt"])
# median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 0, "log2FoldChange.wt"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 1, "log2FoldChange.wt"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 2, "log2FoldChange.wt"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 4, "log2FoldChange.wt"])

table(res.combine.sub2.1$UPF1.type)
# median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type2 == 0, "log2FoldChange.ko"])
# median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 0, "log2FoldChange.ko"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 1, "log2FoldChange.ko"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 2, "log2FoldChange.ko"])
median(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == 4, "log2FoldChange.ko"])

t.test(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "1", "log2FoldChange.wt"],
       res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "2", "log2FoldChange.wt"])
t.test(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "1", "log2FoldChange.wt"],
       res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "4", "log2FoldChange.wt"])

t.test(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "1", "log2FoldChange.wt"],
       res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "1", "log2FoldChange.ko"], paired = T)
t.test(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "2", "log2FoldChange.wt"],
       res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "2", "log2FoldChange.ko"], paired = T)
t.test(res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "4", "log2FoldChange.wt"],
       res.combine.sub2.1[res.combine.sub2.1$UPF1.type == "4", "log2FoldChange.ko"], paired = T)

# 1 m5C not bound UPF1; UPF1 with 2 enrich around; UPF1 with 4 enrich around
# intersect(res.combine.sub2.1[res.combine.sub2.1$UPF1.HeLa == "0", "gene"], res.combine.sub2[, "gene"])
# 
# rbind(res.combine.sub2.1[res.combine.sub2.1$UPF1.HeLa == "0", c("UPF1.type", "log2FoldChange.wt", "log2FoldChange.ko", "gene")],
#       res.combine.sub2[, c("UPF1.type", "log2FoldChange.wt", "log2FoldChange.ko", "gene")]) %>% filter(UPF1.type %in% c(1,2,4)) %>%
#   melt(id = c("gene", "UPF1.type")) %>% 
#   ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
#   theme_classic() + facet_grid(. ~ variable) +
#   scale_color_manual(values = c("grey","green", "blue"))
length(intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)))
res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),
            c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type", "gene")] %>% filter(UPF1.type %in% c(1,2,4)) %>%
  # select(gene, log2FoldChange.wt, log2FoldChange.ko, UPF1.type) %>%
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ variable) +
  scale_color_manual(values = c("grey","#999933","#339933", "#006699"))

a = res.combine[res.combine$gene %in% c("DDX21", "AGRN", "SORT1", "ARHGAP1", "ATF3", "TXLNA", "MAP7D1",
                                        "FAHD1", "THRA", "EHD2", "TMEM214","TRIB3", "H3F3B", 
                                        "SORT1", "MAP7D1", "CAMK2N1", "C12orf43", "MOCOS", "PLCG1", "RHOBTB3", "ZNF696", "RNPEPL1") ,]
b = res.combine[res.combine$gene %in% c("XRCC3", "NEDD8", "H3F3B", "ISOC2", "MRPS25", "ETF1", "RXRB",
                                        "ADGRD1", "FCGRT", "DBI", "NBEAL2","DFFB") ,]
b = res.combine[res.combine$gene %in% c("ADGRD1","ADGRG1","ATP5F1D","CACFD1","CEMIP","DBI","DFFB","ETF1",
  "H3F3B","ISOC2","KLHL25","MRPS25","NBEAL2","RNPEPL1","RXRB","XRCC3","ZFAND1","ZNF696"), ]
d = res.combine[res.combine$gene %in% c("FKBP10","OSMR","NIBAN2","FADS2"), ]

res.combine.sub2.2 = res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),
                                 c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type", "gene")] %>% filter(UPF1.type %in% c(1,2,4))
  # rbind(res.combine.sub2.1[res.combine.sub2.1$UPF1.HeLa == "0",
  #                                             c("UPF1.type", "log2FoldChange.wt", "log2FoldChange.ko", "gene", "label.wt" )],
  #                          res.combine.sub2[, c("UPF1.type", "log2FoldChange.wt", "log2FoldChange.ko", "gene", "label.wt" )])
# res.combine.sub2.2$UPF1.type2 = factor(ifelse(res.combine.sub2.2$UPF1.type ==0, "N", "Y"), levels = c("N", "Y"))
d = res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, ]
dim(res.combine.sub2.2)
table(res.combine.sub2.2$UPF1.type)
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 0, "log2FoldChange.wt"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.wt"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2, "log2FoldChange.wt"])
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2 & res.combine.sub2.2$gene != "TRIB3", "log2FoldChange.wt"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4, "log2FoldChange.wt"])
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4  & !(res.combine.sub2.2$gene %in% c("SLC20A1", "TMEM214", "ATF3")), "log2FoldChange.wt"])

# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 0, "log2FoldChange.ko"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.ko"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2, "log2FoldChange.ko"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4, "log2FoldChange.ko"])
median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4 & !(res.combine.sub2.2$gene %in% c("SLC20A1", "TMEM214", "ATF3")), "log2FoldChange.ko"])

# res.combine.sub2.2[, c("UPF1.type2", "log2FoldChange.wt", "log2FoldChange.ko", "gene")] %>%
#   melt(id = c("gene", "UPF1.type2")) %>% 
#   ggplot(aes(UPF1.type2, value, color = UPF1.type2)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
#   theme_classic() + facet_grid(. ~ variable) +
#   scale_color_manual(values = c("grey","red"))
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "N", "log2FoldChange.wt"])
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "Y", "log2FoldChange.wt"])

table(res.combine.sub2.2$UPF1.type)
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "N", "log2FoldChange.ko"])
# median(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "Y", "log2FoldChange.ko"])
# res.combine.sub2.2 = res.combine.sub2.2[res.combine.sub2.2$log2FoldChange.wt <3, ]
wilcox.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.wt"], 
            res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2, "log2FoldChange.wt"])$p.value
wilcox.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.wt"], 
            res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4, "log2FoldChange.wt"])$p.value

wilcox.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.wt"], 
       res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2, "log2FoldChange.wt"])$p.value
t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.ko"], res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 4, "log2FoldChange.ko"])
t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 1, "log2FoldChange.ko"], res.combine.sub2.2[res.combine.sub2.2$UPF1.type == 2, "log2FoldChange.ko"])

t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "1", "log2FoldChange.wt"],
       res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "1", "log2FoldChange.ko"], paired = T)$p.value
t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "2", "log2FoldChange.wt"],
       res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "2", "log2FoldChange.ko"], paired = T)$p.value
t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "4", "log2FoldChange.wt"],
       res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "4", "log2FoldChange.ko"], paired = T)$p.value

t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "4" & !(res.combine.sub2.2$gene %in% c("SLC20A1", "TMEM214", "ATF3")), "log2FoldChange.wt"],
       res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "4" & !(res.combine.sub2.2$gene %in% c("SLC20A1", "TMEM214", "ATF3")), "log2FoldChange.ko"], paired = T)$p.value

# t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "N", "log2FoldChange.wt"],
#        res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "N", "log2FoldChange.ko"], paired = T)
# t.test(res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "Y", "log2FoldChange.wt"],
#        res.combine.sub2.2[res.combine.sub2.2$UPF1.type2 == "Y", "log2FoldChange.ko"], paired = T)
res.combine.sub2.2[res.combine.sub2.2$UPF1.type == "2", ]

######### final plot figs ########
#total
p1 <- res.combine[, c("log2FoldChange.wt", "UPF1.type")] %>% mutate(gene = rownames(res.combine)) %>% 
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ variable) +
  scale_color_manual(values = c("grey","#999933","#339933", "#006699")) + ylim(c(-3,3))

#not around
p2 <- res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V8)), rownames(res.combine)),
                  c("log2FoldChange.wt", "UPF1.type", "gene")] %>% filter(UPF1.type != 0) %>%
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ variable) +
  scale_color_manual(values = c("#999933","#339933", "#006699")) + ylim(c(-3,3))

#around
p3 <- res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),
                  c("log2FoldChange.wt", "UPF1.type", "gene")] %>% filter(UPF1.type != 0) %>%
  melt(id = c("gene", "UPF1.type")) %>% #head()
  ggplot(aes(UPF1.type, value, color = UPF1.type)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.5) + 
  theme_classic() + facet_grid(. ~ variable) + 
  scale_color_manual(values = c("#999933","#339933", "#006699")) + ylim(c(-3,3))

pdf("/Users/peihong 1/Desktop/projects/2020m5C/files/figures/fig3_UPF1kd_box1_v3.pdf", height = 5, width = 5)
ggplot2.multiplot(p1, p2, p3)
dev.off()


## total
p1 <- res.combine[, c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type")] %>% mutate(gene = rownames(res.combine)) %>% 
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(variable, value, color = UPF1.type, linetype = variable)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ UPF1.type) +
  scale_color_manual(values = c("grey","#999933","#339933", "#006699")) + ylim(c(-3,3))

## not around
p2 <- res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V8)), rownames(res.combine)),
                  c("log2FoldChange.wt", "log2FoldChange.ko", "UPF1.type", "gene")] %>% filter(UPF1.type != 0) %>%
  melt(id = c("gene", "UPF1.type")) %>% 
  ggplot(aes(variable, value, color = UPF1.type, linetype = variable)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.4) + 
  theme_classic() + facet_grid(. ~ UPF1.type) +
  scale_color_manual(values = c("#999933","#339933", "#006699")) + ylim(c(-3,3))

## around
p3 <- res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),
                  c("log2FoldChange.wt", "log2FoldChange.ko","UPF1.type", "gene")] %>% filter(UPF1.type != 0) %>% 
  melt(id = c("gene", "UPF1.type")) %>% #head()
  ggplot(aes(variable, value, color = UPF1.type, linetype = variable)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.5) + 
  theme_classic() + facet_grid(. ~ UPF1.type) + 
  scale_color_manual(values = c("#999933","#339933", "#006699")) + ylim(c(-3,3)) 

pdf("/Users/peihong 1/Desktop/projects/2020m5C/files/figures/fig3_UPF1kd_box2_v3.pdf", height = 6, width = 9)
ggplot2.multiplot(p1, p2, p3)
dev.off()

table(res.combine[intersect(unique(UPF1_overlap.file.N2$V8), rownames(res.combine)),
                  c("UPF1.type")])

table( res.combine[intersect(setdiff(unique(NSUN2_m5C.file$gene), unique(UPF1_overlap.file.N2$V8)), rownames(res.combine)),
                   c("UPF1.type")])

## not have N2-dependent m5C 
res.combine[setdiff(rownames(res.combine), unique(NSUN2_m5C.file$gene)),
            c("log2FoldChange.wt", "log2FoldChange.ko","UPF1.type", "gene")] %>% filter(UPF1.type != 0) %>% 
  melt(id = c("gene", "UPF1.type")) %>% #head()
  ggplot(aes(variable, value, color = UPF1.type, linetype = variable)) + geom_boxplot(width = 0.6, alpha = 0.5, outlier.size = 0.5) + 
  theme_classic() + facet_grid(. ~ UPF1.type) + 
  scale_color_manual(values = c("#999933","#339933", "#006699")) + ylim(c(-3,3)) 
table(res.combine$UPF1.type)
res.combine.sub2.3 = res.combine[setdiff(rownames(res.combine), unique(NSUN2_m5C.file$gene)),
                                 c("log2FoldChange.wt", "log2FoldChange.ko","UPF1.type", "gene")]
table(res.combine.sub2.3$UPF1.type)
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 1, "log2FoldChange.wt"])
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 2, "log2FoldChange.wt"])
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 4, "log2FoldChange.wt"])

table(res.combine.sub2.3$UPF1.type)
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 1, "log2FoldChange.ko"])
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 2, "log2FoldChange.ko"])
median(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == 4, "log2FoldChange.ko"])

t.test(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "1", "log2FoldChange.wt"],
       res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "1", "log2FoldChange.ko"])$p.value
t.test(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "2", "log2FoldChange.wt"],
       res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "2", "log2FoldChange.ko"])$p.value
t.test(res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "4", "log2FoldChange.wt"],
       res.combine.sub2.3[res.combine.sub2.3$UPF1.type == "4", "log2FoldChange.ko"])$p.value








