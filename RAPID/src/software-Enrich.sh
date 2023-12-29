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
# if [ $studyFile == "NULL" ]; then
# 	ls ${annotFiles}/*.study.annot > /dev/null 2>&1
# 	check=`echo $?`
# 	if [ $check == 0 ]; then
# 		studyFile=`ls ${annotFiles}/*.study.annot`
# 	else
# 		echo "The file.study.annot should be supplied by --annotFiles or -i."
# 		helpdoc; exit 1
# 	fi
# fi
# ls ${annotFiles}/*.function.annot > /dev/null 2>&1
# check=`echo $?`
# if [ $check != 0 ]; then
# 	echo "The file.function.annot should be contented in the directory supplied by --annotFiles."
# 	helpdoc; exit 1
# fi
# studySample=${studyFile%.study.annot}
# studySample=${studySample##*/}

if [ ! -d ${OUTDIR} ]; then
	mkdir ${OUTDIR}
fi

# if [ $MODEL != "NULL" ]; then
	# echo $BACKGROUND
	if [ $BACKGROUND == "N" ]; then
		echo "The -b should be set as one of T|C|G|A" ; helpdoc; exit 1
	else
		# echo sample$'\t'Genome$'\t'MATURE$'\t'MODEL$'\t'EXTEND$'\t'TYPE$'\t'START$'\t'study_overlap$'\t'study_per$'\t'function_overlap$'\t'function_per$'\t'N_study_overlap$'\t'N_study_not_overlap$'\t'N_bg_overlap$'\t'N_bg_not_overlap > ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
		if [ $EXTEND == 0 ]; then
			if [ $MATURE == "pre" ]; then
				# ls ${annotFiles}/*.function.annot| while read i; do 
				# for i in ${annotFiles}/*.function.annot; do
				for i in ${annotFiles}; do
					sample=${i##*/}; sample=${sample%.function.annot}
					study_all=`grep -v intergenic $studyFile | wc -l`
					function_all=`grep -v intergenic ${i} | wc -l`
					study_overlap=`intersectBed -a $studyFile -b ${i} -wa -s | sort -u| grep -v intergenic| wc -l`
					function_overlap=`intersectBed -a ${i} -b ${studyFile} -wa -s |grep -v intergenic| sort -u| wc -l`
					study_per=`awk -v study_all=$study_all -v study_overlap=$study_overlap 'BEGIN{print study_overlap/study_all}'`
					function_per=`awk -v function_all=$function_all -v function_overlap=$function_overlap 'BEGIN{print function_overlap/function_all}'`

					N_study_total=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' ${i} $studyFile| grep -v intergenic| wc -l`
					N_study=$study_overlap
					if [ $study_overlap != 0 ];then
						N_bg=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| grep -v intergenic | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bedtools getfasta -fi $FASTA -s -name -bed - | grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0
					fi
					if [ $N_study_total != 0 ];then
						N_bg_total=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| grep -v intergenic| awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]}' ${software_path}/genome/${Genome}.temp -| cut -f 1-6| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct | bedtools getfasta -fi $FASTA -s -name -bed - -split| grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0; N_bg_total=0
					fi
					echo ${sample}$'\t'${Genome}$'\t'precursor$'\t'${MODEL}$'\t'${EXTEND}$'\t'${TYPE}$'\t'${START}$'\t'${study_overlap}$'\t'${study_per}$'\t'${function_overlap}$'\t'${function_per}$'\t'${N_study}$'\t'$(( N_study_total - N_study ))$'\t'$((N_bg - N_study))$'\t'$((N_bg_total - N_study_total - N_bg + N_study))
					if [[ $outOverlaps == 1 && $study_overlap != 0 ]]; then
						intersectBed -a $studyFile -b ${i} -wa -wb -s | sort -u| grep -v intergenic > ${OUTDIR}/${sample}_${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.sites.txt
					fi
				# done >> ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
				done
			else
				for i in ${annotFiles}; do
					sample=${i##*/}; sample=${sample%.function.annot}
					study_all=`egrep -v "intergenic|intron" $studyFile | wc -l`
					function_all=`egrep -v "intergenic|intron" ${i} | wc -l`
					study_overlap=`intersectBed -a $studyFile -b ${i} -wa -s| egrep -v "intergenic|intron" | sort -u| wc -l`
					function_overlap=`intersectBed -a ${i} -b ${studyFile} -wa -s |egrep -v "intergenic|intron"| sort -u| wc -l`
					study_per=`awk -v study_all=$study_all -v study_overlap=$study_overlap 'BEGIN{print study_overlap/study_all}'`
					function_per=`awk -v function_all=$function_all -v function_overlap=$function_overlap 'BEGIN{print function_overlap/function_all}'`

					N_study_total=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' ${i} $studyFile| egrep -v "intergenic|intron"| wc -l`
					N_study=$study_overlap
					if [ $study_overlap != 0 ];then
						N_bg=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron"| sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bedtools getfasta -fi $FASTA -s -name -bed - | grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0
					fi
					if [ $N_study_total != 0 ];then
						N_bg_total=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron"| awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct | bedtools getfasta -fi $FASTA -s -name -bed - -split| grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0; N_bg_total=0
					fi
					echo ${sample}$'\t'${Genome}$'\t'mature$'\t'${MODEL}$'\t'${EXTEND}$'\t'${TYPE}$'\t'${START}$'\t'${study_overlap}$'\t'${study_per}$'\t'${function_overlap}$'\t'${function_per}$'\t'${N_study}$'\t'$(( N_study_total - N_study ))$'\t'$((N_bg - N_study))$'\t'$((N_bg_total - N_study_total - N_bg + N_study))
					if [[ $outOverlaps == 1 && $study_overlap != 0 ]]; then
						intersectBed -a $studyFile -b ${i} -wa -wb -s| egrep -v "intergenic|intron" | sort -u > ${OUTDIR}/${sample}_${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.sites.txt
					fi
				# done >> ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
				done
			fi
		else
			if [ $MATURE == "pre" ]; then
				for i in ${annotFiles}; do
					sample=${i##*/}; sample=${sample%.function.annot}
					study_all=`grep -v intergenic $studyFile | wc -l`
					function_all=`grep -v intergenic ${i} | wc -l`
					study_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | grep -v intergenic| awk -v OFS="\t" 'NR==FNR{$7=$2; $8=$3; a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -s | grep -v intergenic | sort -u | wc -l`
					function_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | grep -v intergenic| awk -v OFS="\t" 'NR==FNR{$7=$2; $8=$3; a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a - -b $studyFile -wa -s | grep -v intergenic | sort -u | wc -l`
					study_per=`awk -v study_all=$study_all -v study_overlap=$study_overlap 'BEGIN{print study_overlap/study_all}'`
					function_per=`awk -v function_all=$function_all -v function_overlap=$function_overlap 'BEGIN{print function_overlap/function_all}'`

					N_study_total=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' ${i} $studyFile | grep -v intergenic |wc -l`
					N_study=$study_overlap
					if [ $study_overlap != 0 ];then
						N_bg=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | grep -v intergenic| awk -v OFS="\t" 'NR==FNR{$7=$2; $8=$3; a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bedtools getfasta -fi $FASTA -s -name -bed - | grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0
					fi
					if [ $N_study_total != 0 ];then
						N_bg_total=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | grep -v intergenic| awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]}' ${software_path}/genome/${Genome}.temp -| cut -f 1-6| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct | bedtools getfasta -fi $FASTA -s -name -bed - -split| grep -v '>' -| awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0; N_bg_total=0
					fi
					echo ${sample}$'\t'${Genome}$'\t'precursor$'\t'${MODEL}$'\t'${EXTEND}$'\t'${TYPE}$'\t'${START}$'\t'${study_overlap}$'\t'${study_per}$'\t'${function_overlap}$'\t'${function_per}$'\t'${N_study}$'\t'$(( N_study_total - N_study ))$'\t'$((N_bg - N_study))$'\t'$((N_bg_total - N_study_total - N_bg + N_study))
					if [[ $outOverlaps == 1 && $study_overlap != 0 ]]; then
						awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | grep -v intergenic| awk -v OFS="\t" 'NR==FNR{$7=$2; $8=$3; a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -wb -s | grep -v intergenic | sort -u > ${OUTDIR}/${sample}_${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.sites.txt
					fi
				# done >> ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
				done
			else
				for i in ${annotFiles}; do
					sample=${i##*/}; sample=${sample%.function.annot}
					study_all=`egrep -v "intergenic|intron" $studyFile | wc -l`
					function_all=`egrep -v "intergenic|intron" ${i} | wc -l`
					study_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -s | sort -u | egrep -v "intergenic|intron"| wc -l`
					function_overlap=`awk 'NR==FNR{a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a - -b $studyFile -wa -s | sort -u | egrep -v "intergenic|intron"| wc -l`
					study_per=`awk -v study_all=$study_all -v study_overlap=$study_overlap 'BEGIN{print study_overlap/study_all}'`
					function_per=`awk -v function_all=$function_all -v function_overlap=$function_overlap 'BEGIN{print function_overlap/function_all}'`

					N_study_total=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' ${i} $studyFile| egrep -v "intergenic|intron"|wc -l`
					N_study=$study_overlap

					## The selection of background regions is relative to ${software_path}/genome/${Genome}.temp file
					## Change the ${software_path}/genome/${Genome}.temp file to set 3utr/5utr/cds/intron/start-codon/stop-codon as background
					if [ $study_overlap != 0 ];then
						N_bg=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START | sort -k1,1 -k2,2n | bedtools merge -i - -c 4,5,6 -o distinct,distinct,distinct -s | bedtools getfasta -fi $FASTA -s -name -bed - | grep -v '>' -|awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg=0
					fi
					# echo $N_study_total
					if [ $N_study_total != 0 ];then
						N_bg_total=`awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i} | egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$2}' ${software_path}/genome/${Genome}.temp -| awk '{split($7,start,","); split($8,end, ","); for(i=1; i<length(start); i++) {print $1"\t"start[i]"\t"end[i]"\t"$4"\t"$5"\t"$6}}'| sort -k1,1 -k2,2n | bedtools merge -s -i - -c 4,5,6 -o distinct,distinct,distinct | bedtools getfasta -fi $FASTA -s -name -bed - -split| grep -v '>' -|awk -v B=$BACKGROUND 'BEGIN{count=0}{$1=toupper($1); split($1,a,""); for(i=1; i<=length(a); i++){if(a[i] == B) count = count+1}}END{print count}'`
					else
						N_bg_total=0
					fi
					if [[ $outOverlaps == 1 && $study_overlap != 0 ]]; then
						awk 'NR==FNR{if($7 !~ /intergenic|intron/) a[$8]++}NR!=FNR{if($8 in a) print}' $studyFile ${i}| egrep -v "intergenic|intron" | awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"int(($2+1+$3)/2)}' ${software_path}/genome/${Genome}.temp -| cut -f 1-8,11 | Get_sequence_around.py -s - --type $TYPE --seq_len $region --out_type bed6 --seq_start $START |  intersectBed -a $studyFile -b - -wa -wb -s | sort -u | egrep -v "intergenic|intron" > ${OUTDIR}/${sample}_${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.sites.txt
					fi
					echo ${sample}$'\t'${Genome}$'\t'mature$'\t'${MODEL}$'\t'${EXTEND}$'\t'${TYPE}$'\t'${START}$'\t'${study_overlap}$'\t'${study_per}$'\t'${function_overlap}$'\t'${function_per}$'\t'${N_study}$'\t'$(( N_study_total - N_study ))$'\t'$((N_bg - N_study))$'\t'$((N_bg_total - N_study_total - N_bg + N_study))
				# done >> ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
				done

			fi
		fi
	fi
	# Rscript ${software_path}/src/chisqTest.r ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
	# rm ${OUTDIR}/${Genome}_${MATURE}_${MODEL}_${TYPE}_${EXTEND}_${START}.temp.txt
# else

# 	echo "nothing"
# fi
