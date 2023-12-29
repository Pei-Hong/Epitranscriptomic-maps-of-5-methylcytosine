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
	sample=${i##*/}; sample=${sample%.function.annot}
	study_all=`egrep -v "intergenic|intron" $studyFile | wc -l`
	function_all=`egrep -v "intergenic|intron" ${i} | wc -l`
	study_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -s | sort -u | egrep -v "intergenic|intron"| wc -l`
	function_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a - -b $studyFile -wa -s | sort -u | egrep -v "intergenic|intron"| wc -l`
	study_per=`awk -v study_all=$study_all -v study_overlap=$study_overlap 'BEGIN{print study_overlap/study_all}'`
	function_per=`awk -v function_all=$function_all -v function_overlap=$function_overlap 'BEGIN{print function_overlap/function_all}'`

	N_study_total=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{print}' ${i} $studyFile| egrep -v "intergenic|intron"|wc -l` #### YES!
	N_study=$study_overlap  #### NO!
	# echo $study_overlap
	if [ $study_overlap != 0 ];then #### NO!
		N_bg=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bedtools getfasta -fi $FASTA -s -name -bed - | grep -v '>' -|awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
	else
		N_bg=0
	fi
	# echo $N_study_total
	if [ $N_study_total != 0 ];then ### YES!
		N_bg_total=`cat $studyFile| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct | bedtools getfasta -fi $FASTA -s -name -bed - -split| grep -v '>' -|awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
	else
		N_bg_total=0
	fi
	echo ${sample}$'\t'${Genome}$'\t'mature$'\t'${MODEL}$'\t'${EXTEND}$'\t'${TYPE}$'\t'${START}$'\t'${study_overlap}$'\t'${study_per}$'\t'${function_overlap}$'\t'${function_per}$'\t'${N_study}$'\t'$(( N_study_total - N_study ))$'\t'$((N_bg - N_study))$'\t'$((N_bg_total - N_study_total - N_bg + N_study))

	## total functional A 
	awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type around --seq_len 20 --out_type bed6 --seq_start 0 | sort -k1,1 -k2,2n | awk -v OFS="\t" '{split($7,a,":");$4=a[2]; print}' |bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bed6ToBed12.py - | bedtools getfasta -fi $FASTA -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > ${OUTDIR}/${sample}_total_functional_A.tmp
	# awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' GSE73405..miCLIP..HepG2..Control.study.annot RBPs_annot/IGF2BP1.combine.mid.function.annot| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' /picb/rnomics4/peihong/0_script/software/genome/hg38.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type around --seq_len 20 --out_type bed6 --seq_start 0 | sort -k1,1 -k2,2n | awk -v OFS="\t" '{split($7,a,":");$4=a[2]; print}' |bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bed6ToBed12.py - | bedtools getfasta -fi /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > total_functional_A.txt
	## total A on the same RNA
	awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($4,t,":"); $4=t[2]; split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct |bed6ToBed12.py - | bedtools getfasta -fi $FASTA -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > ${OUTDIR}/${sample}_total_A.tmp 
	# awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' GSE73405..miCLIP..HepG2..Control.study.annot RBPs_annot/IGF2BP1.combine.mid.function.annot | egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' /picb/rnomics4/peihong/0_script/software/genome/hg38.temp -| awk '{split($4,t,":"); $4=t[2]; split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct |bed6ToBed12.py - | bedtools getfasta -fi /data/rnomics8/peihong/2020_m5C/fasta_chromsome/hg38_chromsomes.fa -s -name -bed - -split| paste - - |sed "s/(/\t/g" | awk -v B=A '{split($1,t,">"); count=0; $3=toupper($3); split($3,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1} print t[2]"\t"count}' > total_A.txt 
	## overlapped m6A
	awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -s | sort -u | egrep -v "intergenic|intron" | cut -f 8 | sort | uniq -c | awk '{print $2"\t"$1}' > ${OUTDIR}/${sample}_overlap_m6A.tmp 
	# awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' GSE73405..miCLIP..HepG2..Control.study.annot RBPs_annot/IGF2BP1.combine.mid.function.annot| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' /picb/rnomics4/peihong/0_script/software/genome/hg38.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type around --seq_len 20 --out_type bed6 --seq_start 0 |  intersectBed -a GSE73405..miCLIP..HepG2..Control.study.annot -b - -wa -s | sort -u | egrep -v "intergenic|intron" | cut -f 8 | sort | uniq -c | awk '{print $2"\t"$1}' > overlap_m6A.tmp
	## total m6A
	awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' ${i} $studyFile| egrep -v "intergenic|intron" | cut -f8 | sort | uniq -c | awk '{print $2"\t"$1}' > ${OUTDIR}/${sample}_total_m6A.tmp
	# awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' RBPs_annot/IGF2BP1.combine.mid.function.annot GSE73405..miCLIP..HepG2..Control.study.annot| egrep -v "intergenic|intron" | cut -f8 | sort | uniq -c | awk '{print $2"\t"$1}' > total_m6A.tmp
	awk 'NR ==FNR{a[$1]=$0}NR!=FNR{if($1 in a) print a[$1]"\t"$2}' ${OUTDIR}/${sample}_total_A.tmp ${OUTDIR}/${sample}_total_functional_A.tmp | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1 in a) print $0"\t"a[$1]}' ${OUTDIR}/${sample}_total_m6A.tmp - | awk 'NR==FNR{a[$1]=$2}NR!=FNR{if($1in a) print $0"\t"a[$1]; else print $0"\t0"}' ${OUTDIR}/${sample}_overlap_m6A.tmp - | awk '{print $1"\t"$5"\t"$4-$5"\t"$3-$5"\t"$2-$3-$4+$5}' > ${OUTDIR}/${sample}.allgene
	rm ${OUTDIR}/${sample}_total_A.tmp ${OUTDIR}/${sample}_total_functional_A.tmp ${OUTDIR}/${sample}_total_m6A.tmp ${OUTDIR}/${sample}_overlap_m6A.tmp
	
done