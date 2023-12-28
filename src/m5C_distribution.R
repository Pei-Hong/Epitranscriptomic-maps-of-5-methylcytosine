library(ggplot2)
library(dplyr)
library(data.table)

## C
select.C.dis = fread("cat \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/0_15sample/15sample_C.dis\" | cut -f 4,12-15 ",
                     header = F, sep = "\t",stringsAsFactors = F)
names(select.C.dis) = c("gene", "region", "total.len", "len")
select.C.dis$label = "C"
gene_len = select.C.dis[!duplicated(select.C.dis[, c("gene", "total.len")]), ]
l.3utr = mean(unlist(gene_len[gene_len$region == "3utr", "total.len"]))
l.5utr = mean(unlist(gene_len[gene_len$region == "5utr", "total.len"]))
l.cds = mean(unlist(gene_len[gene_len$region == "cds", "total.len"]))
l.ave = (l.3utr +l.5utr +l.cds)/3

select.C.dis$len.nor = ifelse(select.C.dis$region == "3utr", l.5utr/l.ave + l.cds/l.ave + (select.C.dis$len / select.C.dis$total.len)*l.3utr/l.ave,
                              ifelse(select.C.dis$region == "5utr", select.C.dis$len  / select.C.dis$total.len *l.5utr/l.ave,
                                     l.5utr/l.ave + (select.C.dis$len  / select.C.dis$total.len) *l.cds/l.ave))

dat<-with(density(select.C.dis$len.nor), data.frame(x,y), adjust = 0.2)

ggplot(select.C.dis, aes(len.nor)) + 
  geom_density(adjust = 0.2, size = 0.8, outline.type = "full") +
  geom_vline(xintercept = c(l.5utr/l.ave, l.cds/l.ave + l.5utr/l.ave), linetype="dashed", color = "grey")


## m5C
select.m5C.dis = fread("cat \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/0_15sample/15sample_m5C.dis\" | cut -f 4,12-15 ",
                            header = F, sep = "\t",stringsAsFactors = F)
names(select.m5C.dis) = c("gene", "region", "total.len", "len")
select.m5C.dis$label = "m5C"

select.m5C.dis$len.nor = ifelse(select.m5C.dis$region == "3utr", l.5utr/l.ave + l.cds/l.ave + (select.m5C.dis$len / select.m5C.dis$total.len)*l.3utr/l.ave,
                                     ifelse(select.m5C.dis$region == "5utr", select.m5C.dis$len  / select.m5C.dis$total.len *l.5utr/l.ave,
                                            l.5utr/l.ave + (select.m5C.dis$len  / select.m5C.dis$total.len) *l.cds/l.ave))

dat<-with(density(select.m5C.dis$len.nor), data.frame(x,y), adjust = 0.2)

ggplot(select.m5C.dis, aes(len.nor)) + 
  geom_density(adjust = 0.2, size = 0.8, outline.type = "full") +
  geom_vline(xintercept = c(l.5utr/l.ave, l.cds/l.ave + l.5utr/l.ave), linetype="dashed", color = "grey")

## combine
combine.C = rbind(select.m5C.dis,
                  select.C.dis)
combine.C$type = factor(combine.C$label, levels = c("C", "m5C"))
ggplot(combine.C, aes(len.nor, color = label)) + 
  geom_density(adjust = 0.2, size = 0.8, outline.type = "full") +
  geom_vline(xintercept = c(l.5utr/l.ave, l.cds/l.ave + l.5utr/l.ave), 
             linetype="dashed", color = "grey") +
  theme_classic() +
  scale_color_manual(values = c("grey", "#99BBEA"))



##################### enriched RBPs #############
HeLa_select.m5C.dis = fread("cat \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/eCLIP_HepG2/enriched_RBP_m5C.dis\" | cut -f 12-15 ",
                            header = F, sep = "\t",stringsAsFactors = F)
names(HeLa_select.m5C.dis) = c("region", "total.len", "len", "RBP")
HeLa_select.m5C.dis$label = "m5C"


HeLa_select.m5C.dis$len.nor = ifelse(HeLa_select.m5C.dis$region == "3utr", l.5utr/l.ave + l.cds/l.ave + (HeLa_select.m5C.dis$len / HeLa_select.m5C.dis$total.len)*l.3utr/l.ave,
                                     ifelse(HeLa_select.m5C.dis$region == "5utr", HeLa_select.m5C.dis$len  / HeLa_select.m5C.dis$total.len *l.5utr/l.ave,
                                            l.5utr/l.ave + (HeLa_select.m5C.dis$len  / HeLa_select.m5C.dis$total.len) *l.cds/l.ave))
HeLa_select.m5C.dis$RBP = factor(HeLa_select.m5C.dis$RBP, 
                                 levels = c("AQR", "BCLAF1", "BUD13", "GRWD1", "EIF3H", "PRPF8", "U2AF2", "IGF2BP1", "FXR2","PPIG",
                                            "DDX3X","UPF1","RPS3","LARP4","NCBP2","DHX30"))

dat<-with(density(HeLa_select.m5C.dis$len.nor), data.frame(x,y), adjust = 0.2)
ymax = 1.3
ggplot() + 
  geom_area(data = dat[dat$x <= (l.5utr/l.ave), ],
            aes(x=x,y=ymax),fill="#E6F5FF") +
  geom_area(data = dat[dat$x < (l.cds/l.ave + l.5utr/l.ave) & 
                         dat$x > (l.5utr/l.ave), ], 
            aes(x=x,y=ymax), fill="#FFF4F3") +
  geom_area(data = dat[dat$x >= (l.cds/l.ave + l.5utr/l.ave), ], 
            aes(x=x,y=ymax),fill="#FFF3F9") +
  geom_density(data = HeLa_select.m5C.dis, aes(len.nor), adjust = 0.2, size = 0.8, outline.type = "full") + 
  theme_classic() +
  geom_vline(xintercept = c(l.5utr/l.ave, l.cds/l.ave + l.5utr/l.ave), linetype="dashed", color = "grey") + 
  xlim(0,3) + facet_wrap(. ~ RBP,nrow = 1, scales = 'free_y')+ 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 15),
        strip.background = element_rect(fill = 'grey', color = 'grey'),
        strip.text = element_text(size = 14)) + ylim(0, ymax)

## RBP
HeLa_select.RBP.dis = fread("cat \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/eCLIP_HepG2/enriched_RBP_RBP.dis\" | cut -f 12-15 ",
                            header = F, sep = "\t",stringsAsFactors = F)
names(HeLa_select.RBP.dis) = c( "region", "total.len", "len","RBP")
HeLa_select.RBP.dis$label = "m5C"

HeLa_select.RBP.dis$len.nor = ifelse(HeLa_select.RBP.dis$region == "3utr", l.5utr/l.ave + l.cds/l.ave + (HeLa_select.RBP.dis$len / HeLa_select.RBP.dis$total.len)*l.3utr/l.ave,
                                     ifelse(HeLa_select.RBP.dis$region == "5utr", HeLa_select.RBP.dis$len  / HeLa_select.RBP.dis$total.len *l.5utr/l.ave,
                                            l.5utr/l.ave + (HeLa_select.RBP.dis$len  / HeLa_select.RBP.dis$total.len) *l.cds/l.ave))
HeLa_select.RBP.dis$RBP = factor(HeLa_select.RBP.dis$RBP, 
                                 levels = c("AQR", "BCLAF1", "BUD13", "GRWD1", "EIF3H", "PRPF8", "U2AF2", "IGF2BP1", "FXR2","PPIG",
                                            "DDX3X","UPF1","RPS3","LARP4","NCBP2","DHX30"))
dat<-with(density(HeLa_select.RBP.dis$len.nor), data.frame(x,y), adjust = 0.2)
ymax = 6.5

ggplot() + 
  geom_area(data = dat[dat$x <= (l.5utr/l.ave), ],
            aes(x=x,y=ymax),fill="#E6F5FF") +
  geom_area(data = dat[dat$x < (l.cds/l.ave + l.5utr/l.ave) & 
                         dat$x > (l.5utr/l.ave), ], 
            aes(x=x,y=ymax), fill="#FFF4F3") +
  geom_area(data = dat[dat$x >= (l.cds/l.ave + l.5utr/l.ave), ], 
            aes(x=x,y=ymax),fill="#FFF3F9") +
  geom_density(data = HeLa_select.RBP.dis, aes(len.nor), adjust = 0.2, size = 0.8,outline.type = "full") + theme_classic() +
  geom_vline(xintercept = c(l.5utr/l.ave, l.cds/l.ave + l.5utr/l.ave), linetype="dashed", color = "grey") + 
  xlim(0,3) + facet_wrap(. ~ RBP,nrow = 1, scales = 'free_y')+ 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 15),
        strip.background = element_rect(fill = 'grey', color = 'grey'),
        strip.text = element_text(size = 14)) + ylim(0,ymax)
dev.off()

################## enrich around RBPs #################
distance.C = 1
distance.C = fread("cut -f 7,9 \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/eCLIP_HepG2/enrichRBP.distanceToC.trim\"",
      header = F, sep = "\t", stringsAsFactors = F)
names(distance.C) = c("dis", "sample")
distance.C$type = "C"
# distance.C = fread("cut -f 5,6 \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/bsRNA/70_protein/RBP2/eCLIP_HepG2/sample_separate/DisFromRBP_C_m5C.2.txt\"",
#                    header = T, sep = "\t", stringsAsFactors = F)
# names(distance.C) = c("sample","dis")
# distance.C$type = "C"
# distance.C$region = ifelse(distance.C$dis >=0 , distance.C$dis %/% 5,
#                              ((-1 * distance.C$dis) %/% 5) * (-1))
# distance.C.percentage = distance.C %>% group_by(region, sample) %>% summarise(n = n())
# distance.C.percentage = distance.C.percentage %>% group_by(sample) %>% summarise(sum = sum(n)) %>% as.data.frame() %>%
#   inner_join(distance.C.percentage) %>% mutate(per = n/sum)
# distance.C.percentage$type = "C"
# distance.C.percentage$per = ifelse(distance.C.percentage$region == 0, distance.C.percentage$per /2, distance.C.percentage$per)

distance.m5C = fread("cut -f 7,9 \"/Users/peihong 1/ly-svr5-rnomic8/2020_m5C/61_data_organize/3_20RC3C_5C_wt_recall/3_RBP/1_HeLa/eCLIP_HepG2/selectRBP.distanceTom5C\"",
                     header = F, sep = "\t", stringsAsFactors = F)
names(distance.m5C) = c("dis", "sample")
distance.m5C$type = "m5C"
distance.m5C$sample2 = factor(distance.m5C$sample, 
                              levels = c("DDX3X", "UPF1","RPS3", "NCBP2", "DHX30"))
distance.C$sample2 = factor(distance.C$sample, 
                            levels = c("DDX3X", "UPF1","RPS3", "NCBP2", "DHX30"))
head(distance.C)
head(distance.m5C)

pdf("/Users/peihong 1/Desktop/projects/2020m5C/files/figures/fig2_density_HeLa_v4.pdf",
    height = 4, width = 22)
rbind(distance.C[,c(1,3,4)], distance.m5C[, c(1,3,4)]) %>% 
  ggplot(aes(dis, color = type)) + geom_density(size=1.2)+ 
  theme_classic() + facet_wrap(. ~ sample2, nrow = 1, scales = "free_y") + 
  xlim(c(-500,500)) + #ylim(c(0,0.0021)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey', color = 'grey'),
        strip.text = element_text(size = 14)) +
  scale_color_manual(values = c("grey", "#993333")) #+ coord_cartesian(xlim = c(-300,300))
dev.off()












