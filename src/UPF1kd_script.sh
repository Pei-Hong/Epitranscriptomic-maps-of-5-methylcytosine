#! /bin/bash

file=$1

# 0_fastq
####### cutadaptor
if [[ ! -f 0_fastq/${file}.trim.log ]];then
    echo $file "cutadaptor"
    java -jar /picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/trimmomatic-0.38.jar PE  -threads 15 -phred33 0_fastq/${file}_1.fq.gz 0_fastq/${file}_2.fq.gz -baseout 0_fastq/${file}.fq.gz ILLUMINACLIP:/picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 > 0_fastq/${file}.trim.log 2>&1
else
    a=`grep "Input Read Pairs:" 0_fastq/${file}.trim.log`
    check1=$?
    if [[ $check1 -eq 1 ]]; then
        echo $file "cutadaptor"
        java -jar /picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/trimmomatic-0.38.jar PE  -threads 15 -phred33 0_fastq/${file}_1.fq.gz 0_fastq/${file}_2.fq.gz -baseout 0_fastq/${file}.fq.gz ILLUMINACLIP:/picb/rnomics4/rotation/fuzhican/download/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30 > 0_fastq/${file}.trim.log 2>&1
    fi
fi

####### bowtie2
if [[ ! -f 10_rRNA_bt2/${file}.bt2.log ]];then
    echo $file "bowtie2"
    echo "owtie2 -q --sensitive -a -p 25 --reorder -x /picb/rnomics4/peihong/hg38_index/RNA45S5 --no-unal -1 0_fastq/${file}_1P.fq.gz -2 0_fastq/${file}_2P.fq.gz /
    -S 10_rRNA_bt2/${file}.bt2.sam --un-conc-gz 10_rRNA_bt2/${file} > 10_rRNA_bt2/${file}.bt2.log 2>&1"
    bowtie2 -q --sensitive -a -p 15 --reorder -x /picb/rnomics4/peihong/hg38_index/RNA45S5 --no-unal -1 0_fastq/${file}_1P.fq.gz -2 0_fastq/${file}_2P.fq.gz -S 10_rRNA_bt2/${file}.bt2.sam --un-conc-gz 10_rRNA_bt2/${file} > 10_rRNA_bt2/${file}.bt2.log 2>&1
    ln -s ${file}.1 10_rRNA_bt2/${file}_1.rm.fq.gz
    ln -s ${file}.2 10_rRNA_bt2/${file}_2.rm.fq.gz
else
    a=a=`samtools quickcheck 10_rRNA_bt2/${file}.bt2.sam`
    check1=$?
    if [[ $check1 -eq 1 ]]; then
        echo $file "bowtie2"
        bowtie2 -q --sensitive -a -p 15 --reorder -x /picb/rnomics4/peihong/hg38_index/RNA45S5 --no-unal -1 0_fastq/${file}_1P.fq.gz -2 0_fastq/${file}_2P.fq.gz -S 10_rRNA_bt2/${file}.bt2.sam --un-conc-gz 10_rRNA_bt2/${file} > 10_rRNA_bt2/${file}.bt2.log 2>&1
        ln -s ${file}.1 10_rRNA_bt2/${file}_1.rm.fq.gz
        ln -s ${file}.2 10_rRNA_bt2/${file}_2.rm.fq.gz
    fi
fi

####### hisat2
if [[ ! -f 11_mapping/${file}.bam ]];then
    echo $file "hisat2"
    hisat2 --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k5 --dta -t -p 15 /picb/rnomics1/database/Human/hg38/genome_primary/hg38 --known-splicesite-infile /picb/rnomics4/peihong/hg38_index/genecode_v32/gencode.v32.spsites.txt -1 10_rRNA_bt2/${file}_1.rm.fq.gz -2 10_rRNA_bt2/${file}_2.rm.fq.gz -S 11_mapping/${file}.sam > 11_mapping/${file}.log 2>&1
    samtools view 11_mapping/${file}.sam -bS | samtools sort -@ 5 > 11_mapping/${file}.bam 
    samtools index 11_mapping/${file}.bam 
else
    a=`samtools quickcheck 11_mapping/${file}.bam`
    check1=$?
    if [[ $check1 -eq 1 ]]; then
        echo $file "hisat2"
        hisat2 --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k5 --dta -t -p 15 /picb/rnomics1/database/Human/hg38/genome_primary/hg38 --known-splicesite-infile /picb/rnomics4/peihong/hg38_index/genecode_v32/gencode.v32.spsites.txt -1 10_rRNA_bt2/${file}_1.rm.fq.gz -2 10_rRNA_bt2/${file}_2.rm.fq.gz -S 11_mapping/${file}.sam > 11_mapping/${file}.log 2>&1
    samtools view 11_mapping/${file}.sam -bS | samtools sort -@ 5 > 11_mapping/${file}.bam
    samtools index 11_mapping/${file}.bam
    fi
fi
