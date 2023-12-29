#!/bin/bash
helpdoc(){
    cat <<EOF
Description:
	software-Enrich.sh 
Usage:
    software -m Enrich [options]
    -h/---help  -- help information
Module:         Enrich 
Options: --annotFiles (-i) -g hg38 --fa --mature --model fisher --extend 50 --mid –o object --type around/left/right --start 0  
 --annotFiles -- Annotted files, which is the output for  "Annot" module.
 -i          -- Studied file which is supposed to be annotated by "Annot" module.
 -g          -- Genome reference for "software".
 --fa        -- Fasta file for genome.
 -o          -- Directory for outputs, the default is "software_out".
 --mature    -- "pre" or "mature". The relative distance would be calculated 
                on precursor or mature RNA, the default is "pre".
 --model     -- Statistal model, "chisq" for Chi-Squared test, "binom" for Binomial
                Test and "prop" for Z-test, the default is "prop".
 -b          -- Background nucleotide, if choose "fisher" as the Statistal 
                model, you should set one of "T|C|G|A" as the background 
                nucleotide.
 --extend    -- Int value, extend the region from the middle site, 
                default is 0.
 --type      -- When decide to "extend", you can choose around/left/right 
                to set the direction to extend, effective only when "extend" 
                >= 0. The default is "around".
 --start     -- Int value. When decide to "extend", you can set up the start 
                site of the extension. If set "start > 0", then the original 
                region didn't be considered.
 --outOverlaps -- Choose to output the sites overlapped by functional regions.
EOF
}

if [ $# = 0 ]; then
    helpdoc
    exit 1
fi
software_path=`which software`
export PATH=${software_path%/software}/src:$PATH

#get command line parameters
ARGS=`getopt -o i:g:o:b:h --long fa:,annotFiles:,mature:,model:,extend:,type:,start:,outOverlaps:,help -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-i) studyFile=$2 ; shift 2;;
		-g) Genome=$2 ; shift 2;;
		--annotFiles) annotFiles=$2 ; shift 2;;
		-o) OUTDIR=$2 ; shift 2;;
		--fa) FASTA=$2 ; shift 2;;
		--mature) MATURE=$2 ; shift 2;;
		--model) MODEL=$2 ; shift 2;;
		-b) BACKGROUND=$2 ; shift 2;;
		--extend) EXTEND=$2 ; shift 2;;
		--type) TYPE=$2 ; shift 2;;
		--start) START=$2 ; shift 2;;
		--outOverlaps) outOverlaps=$2 ; shift 2;;
		-h) helpdoc ; exit 1;;
		--help) helpdoc ; exit 1;;
		--)
			shift
			break
			;;
		*) echo "unknown parameter:" $1; helpdoc; exit 1;;
	esac
done
software_path=`which software`
software_path=${software_path%/software}
#### enrichment analysis
region=$EXTEND


if [ ! -d ${OUTDIR} ]; then
	mkdir ${OUTDIR}
fi

for i in ${annotFiles}; do
	## total functional A on the gene
	sample=${i##*/}; sample=${sample%.function.annot}
	awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type around --seq_len 20 --out_type bed6 --seq_start 0 | sort -k1,1 -k2,2n | awk -v OFS="\t" '{split($7,a,":");$4=a[2]; print}' |bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bed6ToBed12.py - | bedtools getfasta -fi $FASTA -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > ${OUTDIR}/${sample}_total_functional_A.tmp

	## total A on the gene
	egrep -v "intergenic|intron" $studyFile | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($4,t,":"); $4=t[2]; split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct |bed6ToBed12.py - | bedtools getfasta -fi $FASTA -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > ${OUTDIR}/${sample}_total_A.tmp 

	## length of gene
	egrep -v "intergenic|intron" $studyFile | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($4,t,":"); $4=t[2]; split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct |bed6ToBed12.py - | bedtools getfasta -fi $FASTA -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">");print t[2]"\t"length($3)}' > ${OUTDIR}/${sample}_gene_len.tmp 

	## overlapped m6A on the gene
	awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -s | sort -u | egrep -v "intergenic|intron" | cut -f 8 | sort | uniq -c | awk '{print $2"\t"$1}' > ${OUTDIR}/${sample}_overlap_m6A.tmp 

	## total m6A on the gene
	egrep -v "intergenic|intron" $studyFile| cut -f8 | sort | uniq -c | awk '{print $2"\t"$1}' > ${OUTDIR}/${sample}_total_m6A.tmp

	awk 'NR ==FNR{a[$1]=$0}NR!=FNR{if($1 in a) print a[$1]"\t"$2}' ${OUTDIR}/${sample}_total_m6A.tmp ${OUTDIR}/${sample}_total_A.tmp |  awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]; else print $0"\t0"}' ${OUTDIR}/${sample}_overlap_m6A.tmp - | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]; else print $0"\t0"}' ${OUTDIR}/${sample}_total_functional_A.tmp - |awk '{print $1"\t"$4"\t"$2-$4"\t"$5-$4"\t"$3-$2-$5+$4}' | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]}' ${OUTDIR}/${sample}_gene_len.tmp - > ${OUTDIR}/${sample}.allgene 

	# awk 'NR ==FNR{a[$1]=$0}NR!=FNR{if($1 in a) print a[$1]"\t"$2}' ${OUTDIR}/${sample}_total_A.tmp ${OUTDIR}/${sample}_total_functional_A.tmp | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]}' ${OUTDIR}/${sample}_total_m6A.tmp - | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]; else print $0"\t0"}' ${OUTDIR}/${sample}_overlap_m6A.tmp - | awk '{print $1"\t"$5"\t"$4-$5"\t"$3-$5"\t"$2-$3-$4+$5}' > ${OUTDIR}/${sample}.allgene
done

