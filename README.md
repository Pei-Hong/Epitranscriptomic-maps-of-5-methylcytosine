# Epitranscriptomic-maps-of-5-methylcytosine
This is the bioinformatical analysis pipline associated with "Epitranscriptomic maps of 5-methylcytosine reveal substrate diversity of NSUN methyltransferases and links to mRNA translation and turnover

## 1. Establishment of union m5C set from public datasets
### 1.1 Requirements
```
FastQC v0.11.9
Trimmomatic v0.38
Bowtie2 v2.3.5
MeRanTK v1.2.1b
HISAT2 v2.1.0
Picard v2.7.1
Samtools v 1.9
FIMO v5.4.1
MEME v5.4.1
```
### 1.2 Data pre-processing
```
for i in *_1.fastq.gz ; do echo ${i%_1.fastq.gz}_2.fastq.gz ${i%_1.fastq.gz}.fq.gz; java -jar /picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/trimmomatic-0.38.jar PE  -threads 30 -phred33 ${i} ${i%_1.fastq.gz}_2.fastq.gz -baseout ./tmp ILLUMINACLIP:/picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 > ${i%_1.fastq.gz}.trim.log 2>&1 ; done
```
### 1.3 Mapping to genome
* Genome mapping of the bsRNA-seq data underwent two-step mapping strategy via script `2020_m5C.sh`. In brief, clean reads were mapped on the 45S rRNA, pre-tRNA, and predicted mature tRNA using meRanT tool (align bsRNA-seq reads to reference using Bowtie2 v2.3.5) in MeRanTK (v1.2.1b) with options (-k 10). Then, paired reads were obtained from unmapped reads using seqkit. Finally, unmapped reads pairs were mapped to the hg38 genome using the meRanGh tool (align bsRNA-seq reads to reference using HISAT2 v2.1.0) in MeRanTK with options (-fmo). Only uniquely mapped reads were retained and used for m5C sites detection.

* The library complexity was analyzed for each datasets and represented by PCR Bottlenecking Coefficient 1 (PBC1), calculated as the number of genomic locations where exactly one read maps uniquely divided by the number of distinct genomic locations to which read maps uniquely. For low complexity datasets (T24 cells, Chen et al., 2019), PCR duplicates were marked and removed by Picard MarkDuplicates (v2.7.1) with the following commands
```
picard MarkDuplicates REMOVE_DUPLICATES=true I=T24_sorted_3C.C2T.bam O=T24_sorted_3C.C2T.markdup.bam M=T24_sorted_3C.C2T.markdup.txt > T24_sorted_3C.C2T.markdup.log 2>&1 &
picard MarkDuplicates REMOVE_DUPLICATES=true I=T24_sorted_3C.G2A.bam O=T24_sorted_3C.G2A.markdup.bam M=T24_sorted_3C.G2A.markdup.txt > T24_sorted_3C.G2A.markdup.log 2>&1 &
```

### 1.4 m5C calling
The discovery of m5C sites from each dataset was performed similarly to Schumann et al., (2020), via a modified script `call_m5C.sh`. Read coverage at each cytosine position in the genome was obtained using the ‘mpileup’ function in samtools (v1.9). The C-to-T mismatches on forward reads and G-to-A mismatches on reverse reads were regarded as unconverted cytosines and called using a custom script with parameters ‘-minBQ 30 --overhang 6’, where reads were filtered by minimum base quality score 30 and removed 6 nt from both 5’ and 3’ terminals to avoid overestimation of non-conversion. 

### 1.5 m5C filtering
Reads containing more than three unconverted cytosines were considered as conversion failure and removed from the bam files (‘3C filter’). In the meantime, the RNA icSHAPE values determined by Sun et al., (2021) in HeLa, HEK293 and HepG2 cells were downloaded from GEO (GSE145805) to confirm that m5C sites flagged by the “3C filter” were in regions with low icSHAPE values, i. e., in highly-structured regions. Then candidate sites with a signal-to-noise ratio < 0.9 (3C/raw; ‘S/N ≥ 0.9’) were dropped to reduce false positives. To retain high-confidence non-conversion sites, the following criteria were applied: (1). Minimum total read coverage was set as 20 (‘20RC’); (2). non-converted C ≥ 3 (‘3C’); (3). C + T coverage ≥ 80% (‘80CT’); (4). Non-conversion of ≥ 10% (‘10MM’); (5). For replicates integration, we referred to Huang et al., (2019). We select m5C candidates detected in at least 2 replicates, and set non-converted C ≥ 5 for sites detected in only one replicate. 

### 1.6 Distribution and enrichment of m5C along mRNA. (R scripts)
Genome annotation of candidate m5C sites were performed via `RAPID` pipeline , which is based on the hg38 GENCODE v32 gene reference from UCSC. The annotation priority is: 5′UTR, CDS, 3′UTR, ncRNA_exonic, intronic and intergenic. Only m5C sites on exonic protein-coding transcripts were used in the following analysis.
```
RAPID -m Build -g hg38 --gtf hg38_GENCODEv32.gtf –iso ## build RAPID 
RAPID -m Annot -i 15sample_m5C.bed -g hg38 -o 15sample_m5C -s
ln -s 15sample_m5C/15sample_m5C.study.annot ./15sample_m5C.annot
egrep "cds|3utr|5utr" 15sample_m5C.annot > 15sample_m5C.annot.mRNA
```
mRNAs with m5C sites were divided into three segments (5′UTR, CDS, and 3′UTR) and normalized according to their average length. All ‘C’ positions on the transcripts with m5C candidates were selected as background control. The relative position of each candidate site and background ‘C’ in the corresponding segment were calculate separately via bellow commands. Figs were plot by R script `m5C_distribution.R`.
For m5C in each sample:
```
RAPID -m Distance --annotFiles 15sample_m5C.annot.mRNA -g hg38 -o 15sample_m5C_dis --start_site 5utr &
RAPID -m Distance --annotFiles 15sample_m5C.annot.mRNA -g hg38 -o 15sample_m5C_dis --start_site 3utr &
RAPID -m Distance --annotFiles 15sample_m5C.annot.mRNA -g hg38 -o 15sample_m5C_dis --start_site cds &
cat 15sample_m5C_dis/*dis > 15sample_m5C.dis
```
For background ‘C’:
```
cut -f 8-9 15sample_m5C.annot.mRNA | sort -u | awk 'NR==FNR{a[$1":"$2]++}NR!=FNR{if($4 in a) print}' - gencode.v32.basic.annotation.bed12 | bedtools getfasta -fi /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa -s -name -bed - -split | awk 'BEGIN{name=0; strand=""}{if($1 ~ /^>/) {split($1,b,">"); split(b[2],c,"("); name=c[1]; strand=c[2]} else {if(strand=="+)") {split($1,a,""); for(i=1;i<=length(a); i++) {if(toupper(a[i]) == "C") print name"\t"i}} else {split($1,a,""); n=1; for(i=length(a); i>=1; i--) {if(toupper(a[i]) == "C") print name"\t"length(a)-i+1; n=n+1}} }}' | awk 'NR==FNR{start[$4]=$2; lens[$4]=$11; starts[$4]=$12; chr[$4]=$1; chrom[$4]=$6}NR!=FNR{if($1 in start) {split(starts[$1], s,","); split(lens[$1], l, ","); s0=start[$1]; l1=0; l2=0; for(i=1; i<=length(s); i++) {s0=start[$1]+s[i]; l2=l1+l[i]; if($2 < l2) {print chr[$1]"\t"s0+$2-1-l1"\t"s0+$2-l1"\t"chr[$1]":"s0+$2-l1"\t0\t"chrom[$1]; break} else{l1=l1+l[i]} }}}' gencode.v32.basic.annotation.bed12 - | awk 'NR==FNR{a[$4]++}NR!=FNR{if(!($4 in a)) print}' 15sample_m5C.annot.mRNA - > 15sample_C.mRNA
RAPID -m Annot -i 15sample_C.mRNA -g hg38 -o 15sample_C -s
RAPID -m Distance --annotFiles 15sample_C.annot.mRNA -g hg38 -o 15sample_C_dis --start_site 3utr &
RAPID -m Distance --annotFiles 15sample_C.annot.mRNA -g hg38 -o 15sample_C_dis --start_site 5utr &
RAPID -m Distance --annotFiles 15sample_C.annot.mRNA -g hg38 -o 15sample_C_dis --start_site cds &
cat 15sample_C_dis/*dis > 15sample_C.dis
```

## 2. NSUN-enzyme dependence of m5C sites.
We collected NSUN2 (Huang et al., 2019; Yang et al., 2017) and NSUN6 depletion datasets (Liu et al., 2021a) for NSUN2/6-dependent m5C sites analysis. The analysis method was referred to in (Chen et al., 2019). In detail, m5C sites with more than 0.05 methylation ratio decrease and methylation level reduced to less than 0.1 in NSUN2 or NSUN6 depleted cells were considered to be catalysed by NSUN2 or NSUN6, respectively. In particular, for NSUN2-dependent m5C we used the union set of Huang et al., (2019), and Yang et al., (2017). 
The sequence of position 0 - +5 of NSUN-dependent m5C sites were then used to calculate position weight matrixes (PWM) of consensus motifs. 
NSUN-dependence of m5C sites in union set can be further predicted by PWM through FIMO (Find Individual Motif Occurrences) in the MEME suite (v 5.4.1). 
```
awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($8":"$9 in a) print a[$8":"$9]"\t"$3; else print $1"\t"$2-100"\t"$3+100"\tintergenic\t0\t"$6"\t"$2-100"\t"$3+100"\t"$2-100"\t"$3+100"\t"$3}' gencode.v32.annotation.temp 15sample_m5C.annot | cut -f1-8,11| Get_sequence_around.py -s - --type around --seq_len 10 --out_type bed12 | cut -f1-12 | bedtools getfasta -fi /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa -s -name -bed - -split  | awk '{if($1 !~ />chr/) {gsub("T|t", "U", $1); print toupper($1)} else print}' > 15sample_m5C.fa
cat 15sample_m5C.fa | paste - - | awk '{if(length($2) == 21) print $1"\n"$2}' > 15sample_m5C.21.fa

fimo --oc NSUN2_motif --verbosity 4 --thresh 0.05 NSUN2.meme 15sample_m5C.21.fa
fimo --oc NSUN6_motif --verbosity 4 --thresh 0.05 NSUN6.meme 15sample_m5C.21.fa

```
For potential NSUN5-dependent sites, we used NSUN5-overexpression and epigenetically silenced NSUN5 data in LN229 cells Janin et al., (2019), and regarded sites with 0.05 methylation ratio increase and methylation level higher than 0.1 in NSUN5-OE group are potential NSUN5-dependent sites. All the consensus motifs are plotted by ggseqlogo (Wagih, 2017) in R script `m5C_distribution.R`.

## 3. Proximity to RBP binding sites.
RBP footprints reported by (Van Nostrand et al., 2020) were downloaded from ENCODE; ALYREF footprints in HeLa cells line were obtained from CLIPdb database (Yang et al., 2015) as a positive control. For ENCODE datasets, the intersections of two biological replicates with fold-enrichment ≥ 4 and p ≤ 10-3, were filtered as significant peaks and the middle sites of the peaks were regarded as RBP binding sites. 
Genes with both RBP binding sites and m5C sites were considered for enrichment analysis. Bins were divided around RBP binding sites with the number of m5C sites and total C within the bins. Then Fisher’s exact test was applied to calculate the enrichment significance. We set up a gradient region as ±20 nt, ±30 nt, ±50 nt, and ±70 nt, to test the influence of bin size. Finally, ±50 nt was used for further analysis. The relative distance of m5C sites to RBP binding sites was identified as well as the background C within the same region to get the distribution of m5C around RBP binding sites.

For HepG2 eCLIP dataset and HepG2 m5C
```
awk 'NR==FNR{a[$1]++}NR!=FNR{if($4 in a) print}' Zhang_HepG2.union 15sample_m5C.annot.mRNA > Zhang_HepG2.mRNA

for region in 20 30 50 70; do RAPID -m Enrich -i Zhang_HepG2.mRNA --annotFiles HepG2_RBP_mRNA -g hg38 --fa /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa --mature --extend $region -o HepG2_eCLIP_HepG2_2 -p 5 --outOverlaps -b C; done

for i in 20 30 50 70; do for f in HepG2_eCLIP_HepG2_2/*_${i}_0.txt; do awk -v r=$i '{print $1"\t"$16"\t"$19"\t"r}' $f | sed '1d' ; done; done | grep -v UPF1..HeLa.combine.mid > HepG2.all.prop_2

```
For HepG2 eCLIP dataset and HeLa m5C
```
for region in 20 30 50 70; do RAPID -m Enrich -i HeLa5.annot.mRNA --annotFiles HepG2_RBP_mRNA -g hg38 --fa /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa --mature --extend $region -o HeLa5_eCLIP_HepG2_2 -p 10 --outOverlaps -b C; done

for i in 20 30 50 70; do for f in HeLa5_eCLIP_HepG2_2/*_${i}_0.txt; do awk -v r=$i '{print $1"\t"$16"\t"$19"\t"r}' $f | sed '1d' ; done; done | grep -v UPF1..HeLa.combine.mid > HeLa.all.prop_2

```

For HepG2/K562 eCLIP dataset and union m5C
```
for region in 20 30 50 70; do RAPID -m Enrich -i 15sample_m5C.annot.mRNA --annotFiles HepG2_RBP_mRNA -g hg38 --fa /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa --mature --extend $region -o all_eCLIP_HepG2 -p 10 --outOverlaps -b C; done &
for region in 20 30 50 70; do RAPID -m Enrich -i 15sample_m5C.annot.mRNA --annotFiles K562_RBP_mRNA -g hg38 --fa /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa --mature --extend $region -o all_eCLIP_K562 -p 10 --outOverlaps -b C; done &

for i in 20 30 50 70; do for f in all_eCLIP_HepG2/*_${i}_0.txt; do awk -v r=$i '{print $1"\t"$16"\t"$19"\t"r}' $f | sed '1d' ; done; done | grep -v UPF1..HeLa.combine.mid > all_eCLIP_HepG2.prop
for i in 20 30 50 70; do for f in all_eCLIP_K562/*_${i}_0.txt; do awk -v r=$i '{print $1"\t"$16"\t"$19"\t"r}' $f | sed '1d' ; done; done | grep -v UPF1..HeLa.combine.mid > all_eCLIP_K562.prop

```

## 4. UPF1 knockdown RNA-seq and differential expression analysis
The analysis pipeline can be found in `UPF1kd_script.sh`.
The raw reads were subjected to FastQC (v0.11.9). Low-quality bases and adaptor sequences were removed using Trimmomatic (v0.38) with options (ILLUMINCLIP:Adapter.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW: 4:20 MINLEN:50). Clean reads were mapped to 4S rRNA using bowtie2 (v2.3.5) with parameters (‘-q --sensitive --reorder --no-unal --un-conc-gz’). Unmapped reads were mapped to hg38 genome using HISAT2 (v2.1.0) with parameters (‘--no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k5 --dta’). The mapped reads were processed to featureCounts (v2.0.1) for feature counting. 
To quantify RNA, fragment abundance in genes with ≥50 mapped reads were selected and normalized using the size factors estimated by the median of all genes implemented in the DESeq2 (v1.30.1) Bioconductor package. Differential expression analysis was performed by DESeq2. Genes with log2(fold change) ≥ 1.2 and adjusted p-value ≤ 0.05 were regarded as differentially expressed genes.

## 5. License
Copyright (C) 2022 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use

