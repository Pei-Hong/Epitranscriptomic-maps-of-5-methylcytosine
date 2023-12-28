#!/bin/bash

# /Users/peihong 1/ly-svr5-rnomic4/0_script/2020_m5C.sh
# sample_list_PE=(00_Huang2019/Heart_N1 00_Huang2019/Heart_N5 00_Huang2019/Huang_HEK293T 00_Huang2019/Huang_Hela_ctl_r1 00_Huang2019/Huang_Hela_ctl_r2 
# 08_Chen2019/N0466 08_Chen2019/T0466 02_Janin2019/LN229_EV2 02_Janin2019/LN229_OE3 01_Yang2017/Yang_Hela_ctl_r1 01_Yang2017/Yang_Hela_ctl_r2)
# sample_list_PE=(00_Huang2019/Heart_N1 00_Huang2019/Heart_N5 00_Huang2019/Huang_HEK293T 00_Huang2019/Huang_Hela_ctl_r1 00_Huang2019/Huang_Hela_ctl_r2 08_Chen2019/N0466 08_Chen2019/T0466 02_Janin2019/LN229_EV2 02_Janin2019/LN229_OE3 01_Yang2017/Yang_Hela_ctl_r1 01_Yang2017/Yang_Hela_ctl_r2)
# sample_list_PE=(09_Thomas/Thomas_Hela_B1_S1)
# sample_list_PE=(01_Yang2017/Yang_Hela_ctl_r1 01_Yang2017/Yang_Hela_ctl_r2)
# sample_list_PE=(00_Huang2019/Huang_HEK293T 00_Huang2019/Huang_Hela_ctl_r1 00_Huang2019/Huang_Hela_ctl_r2 
# 00_Huang2019/Lung_N1 00_Huang2019/Lung_N7 00_Huang2019/Frontal_cortex_N4 00_Huang2019/Heart_N1 00_Huang2019/Heart_N5 00_Huang2019/Liver_N1 00_Huang2019/Liver_N6 
# 00_Huang2019/Muscle_N2 00_Huang2019/Muscle_N5 00_Huang2019/Spleen_N1 00_Huang2019/Spleen_N6 00_Huang2019/Testis_N2 00_Huang2019/Testis_N7 09_Thomas/Thomas_Hela_B1_S1 09_Thomas/B2_S2 09_Thomas/B3_S3)
# sample_list_PE=(08_Chen2019/N0156 08_Chen2019/N0453 08_Chen2019/N0474 08_Chen2019/N0531 08_Chen2019/T0156 08_Chen2019/T0453 08_Chen2019/T0474 08_Chen2019/T0531)
# sample_list_PE=(02_Janin2019/LN229_EV2 02_Janin2019/LN229_EV1 02_Janin2019/LN229_EV3 02_Janin2019/LN229_OE1 02_Janin2019/LN229_OE2 02_Janin2019/LN229_OE3 
# 08_Chen2019/N0466 08_Chen2019/T0466 08_Chen2019/N0156 08_Chen2019/N0453 08_Chen2019/N0474 08_Chen2019/N0531 08_Chen2019/T0156 08_Chen2019/T0453 08_Chen2019/T0474 08_Chen2019/T0531)
# sample_list_PE=(09_Thomas/C1_S5 09_Thomas/C2_S6 09_Thomas/C3_S7 09_Thomas/E1_S9 09_Thomas/E2_S10 09_Thomas/C4_S8 09_Thomas/E3_S11 09_Thomas/E4_S12)
# sample_list_SE=(05_Dai/Dai_293T_r2_1 07_Wei2018/Wei_MCF10A_1 03_Sajini2019/Sajini_H9_rep1_trimmed 03_Sajini2019/Sajini_HEK_rep1_trimmed 00_Huang2019/Heart_N5.ctl 
# 00_Huang2019/Huang_Hela_ctl_r1.ctl 07_Wei2018/Wei_MCF10A_1.ctl 05_Dai/Dai_293T_r3_1)
# sample_list_PE=(010_Zhang/Zhang_HEK293T_WT_rep2 010_Zhang/Zhang_HEK293T_cyto 
# 010_Zhang/Zhang_HEK293T_N6_ko_rep1 010_Zhang/Zhang_HEK293T_N6_ko_rep2 010_Zhang/Zhang_HEK293T_nul )  # 010_Zhang/Zhang_HEK293T_WT_rep1 
# sample_list_PE=(01_Huang/Lung_M1 01_Huang/Small_intestine_F1 01_Huang/Small_intestine_M1 01_Huang/Spleen_F1 01_Huang/Spleen_M1 01_Huang/Muscle_F2 01_Huang/Muscle_F3 01_Huang/Muscle_F4)
# 010_Zhang/Zhang_HepG2_rep1 010_Zhang/Zhang_HepG2_rep  08_Chen2019/T24 YBX1_RIP_bs_T24
sample_list_PE=(0_fastq/Sample_1_S1 0_fastq/Sample_4_S4)
#  00_fastq/Zhang_HEK293T_NSUN2ko_r2 00_fastq/Zhang_HEK293T_NSUN2ko6ko_r1
# sample_list_PE=(00_fastq/Huang_Hela_NSUN2ko_r1 00_fastq/Huang_Hela_NSUN2ko_r2) # 00_fastq/Huang_Hela_NSUN2ko_r1 00_fastq/Huang_Hela_NSUN2ko_r2 00_fastq/Yang_Hela_NSUN2ko_r1 00_fastq/Yang_Hela_NSUN2ko_r2
# sample_list_SE=(00_fastq/H9_wt1_rep1 00_fastq/H9_wt1_rep2 00_fastq/H9_ko1_rep1 00_fastq/H9_ko1_rep2)
# sample_list_SE=(02_Blanco/skinCancer_NSUN2ko1_trimmed 02_Blanco/skinCancer_NSUN2ko2_trimmed 02_Blanco/skinCancer_NSUN2ko3_trimmed 02_Blanco/skinCancer_NSUN2ko4_trimmed
#    02_Blanco/skinCancer_wt1_trimmed 02_Blanco/skinCancer_wt2_trimmed 02_Blanco/skinCancer_wt3_trimmed 02_Blanco/skinCancer_wt4_trimmed)
# sample_list_SE=(02_Blanco/skinCancer_NSUN2ko1_trimmed 02_Blanco/skinCancer_NSUN2ko2_trimmed 02_Blanco/skinCancer_NSUN2ko3_trimmed 02_Blanco/skinCancer_NSUN2ko4_trimmed
#     02_Blanco/skinCancer_wt1_trimmed 02_Blanco/skinCancer_wt2_trimmed 02_Blanco/skinCancer_wt3_trimmed 02_Blanco/skinCancer_wt4_trimmed)
# sample_list_SE=(02_Blanco/skinCancer_NSUN2ko3_trimmed 02_Blanco/skinCancer_NSUN2ko4_trimmed
#     02_Blanco/skinCancer_wt1_trimmed 02_Blanco/skinCancer_wt2_trimmed 02_Blanco/skinCancer_wt3_trimmed 02_Blanco/skinCancer_wt4_trimmed)
sample_list_SE=()
# sample_list_SE=(01_Tommaso/H9_ko1_rep3_trimmed 01_Tommaso/H9_ko1_rep4_trimmed 01_Tommaso/H9_wt1_rep3_trimmed 01_Tommaso/H9_wt1_rep4_trimmed)
# sample_list_SE=(01_Tommaso/H9_ko2_rep1_trimmed 01_Tommaso/H9_ko2_rep2_trimmed 01_Tommaso/H9_ko2_rep3_trimmed 01_Tommaso/H9_ko2_rep4_trimmed)

out_dir=103_mapping
# hg38
# index=/data/rnomics8/peihong/2020_m5C/fasta_chromsome/BSgenomeIDX
# gtf=/data/rnomics8/peihong/2020_m5C/fasta_chromsome/gencode.v32.annotation.gtf
# hg38
index=/data/rnomics8/peihong/2020_m5C/fasta_chromsome_ctl/BSgenomeIDX
gtf=/data/rnomics8/peihong/2020_m5C/fasta_chromsome_ctl/hg38_chromsomes_ctl.gtf

# ## mm10
# index=/data/rnomics8/peihong/2020_m5C_mouse/fasta_chrom/BSgenomeIDX
# gtf=/picb/rnomics4/peihong/mm10_index/gencode_v23/gencode.vM23.annotation.gtf

for ((j=0; j<${#sample_list_PE[@]}; j++)); do
    i=${sample_list_PE[$j]}
    echo ${i} ${out_dir}/${i#*/}_tRNA_rRNA/${i#*/}_2P_unmapped.fq.gz
    ## 用paired reads
    if [[ ! -f ${out_dir}/${i#*/}_tRNA_rRNA/fq_filted/${i#*/}_2P_unmapped.fq.gz ]]; then
        /picb/rnomics3/peihong/software/seqkit pair -1 ${out_dir}/${i#*/}_tRNA_rRNA/${i#*/}_1P_unmapped.fq.gz -2 ${out_dir}/${i#*/}_tRNA_rRNA/${i#*/}_2P_unmapped.fq.gz -O ${out_dir}/${i#*/}_tRNA_rRNA/fq_filted
    fi
    /picb/rnomics3/peihong/software/meRanTK-1.2.1b/meRanGh align -o ${out_dir}/${i#*/} -f ${out_dir}/${i#*/}_tRNA_rRNA/fq_filted/${i#*/}_2P_unmapped.fq.gz -r ${out_dir}/${i#*/}_tRNA_rRNA/fq_filted/${i#*/}_1P_unmapped.fq.gz -S ${i#*/}.sam -un -ud ${out_dir}/${i#*/} -MM -fmo \
    -id ${index} \
    -GTF ${gtf} -threads 35 > ${out_dir}/${i#*/}.meRanGh.log 2>&1

done

for ((j=0; j<${#sample_list_SE[@]}; j++)); do
    i=${sample_list_SE[$j]}
    if [[ ! -f 102_mapping/${i#*/}_tRNA_rRNA/trim/ ]]; then
        mkdir 102_mapping/${i#*/}_tRNA_rRNA/trim/
    fi
    zcat 102_mapping/${i#*/}_tRNA_rRNA/${i#*/}_unmapped.fq.gz | awk 'BEGIN{n=0}{n=n+1; if(n==2) {print substr($1,1,(length($1)-6)); n=0} else print $0}' | gzip > 102_mapping/${i#*/}_tRNA_rRNA/trim/${i#*/}_unmapped.fq.gz
    echo ${i} 102_mapping/${i#*/}_tRNA_rRNA/trim/${i#*/}_unmapped.fq.gz
    /picb/rnomics3/peihong/software/meRanTK-1.2.1b/meRanGh align -o ${out_dir}/${i#*/} -f 102_mapping/${i#*/}_tRNA_rRNA/trim/${i#*/}_unmapped.fq.gz -S ${i#*/}.sam -un -ud ${out_dir}/${i#*/} -MM -fmo \
    -id ${index} \
    -GTF ${gtf} -threads 35 > ${out_dir}/${i#*/}.meRanGh.2.log 2>&1

    
done



 
