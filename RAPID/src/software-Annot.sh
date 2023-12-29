#!/bin/bash
helpdoc(){
    cat <<EOF
Description:
	SCAPRUTE pipeline. A deep learning-embedded pipeline that captures polyadenylation information from 3 prime tag-based RNA-seq of single cells
Usage:
    software -m <build|Annot|Distance|Enrich> [options]
    -h/---help  -- help information
Module:         Annot
Options: -i studiedFile.bed -f FunctionalFile.bed -g hg38 -o object 
 -i          -- bed6 format file for sites studied.  
 -f          -- bed6 format file for regions with biological functions. You 
                can provide a directory content files.
 -g          -- The name of genome reference. You can build them use "software Build"
 -s			 -- Annoted sites by strand.
 -o          -- Directory for outputs
EOF
}

if [ $# = 0 ]; then
    helpdoc
    exit 1
fi

#get command line parameters
ARGS=`getopt -o i:f:g:o:h:s: --long dir:,help -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-i) StudyFile=$2 ; shift 2;;
		-f) FunctionFile=$2 ; shift 2;;
		-g) Genome=$2 ; shift 2;;
		-s) Strand=$2 ; shift 2;;
		-o) OUTDIR=$2 ; shift 2;;
		-h) helpdoc ; exit 1;;
		--help) helpdoc ; exit 1;;
		--)
			shift
			break
			;;
		*) echo "unknown parameter:" $1; helpdoc; exit 1;;
	esac
done

###### annotation for bed
software_path=`which software`
software_path=${software_path%/software}
mkdir -p $OUTDIR

if [ $FunctionFile != "NULL" ]; then
	type=function
	if [ -f $FunctionFile ]; then
		out_file=${FunctionFile##*/}
		d=${OUTDIR}/${out_file%.bed}
		${software_path}/src/Get_intersectRegion2.py -g ${software_path}/genome/${Genome}_region -f $FunctionFile -o ${d} -t $type -s $Strand
	else
		for file in ${FunctionFile}/*.bed; do
			out_file=${file##*/}
			d=${OUTDIR}/${out_file%.bed}
			${software_path}/src/Get_intersectRegion2.py -g ${software_path}/genome/${Genome}_region -f $file -o ${d} -t $type -s $Strand
		done
	# ls ${FunctionFile} | while read i; do 
	# 	# mkdir -p ${OUTDIR}/${i%.bed}
	# 	d=${OUTDIR}/${i%.bed}
	# 	if [ ! -f ${FunctionFile} ]; then
	# 		i=${FunctionFile}/${i}
	# 	fi
	# 	${software_path}/src/Get_intersectRegion2.py -g ${software_path}/genome/${Genome}_region -f $i -o ${d} -t $type -s $Strand
		# cut -f1-6 ${i} | intersectBed -a - -b ${software_path}/genome/genome_region/cds.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14 | sort -k1,4 -u > ${d}/cds.txt 
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/cds.txt ${i} | cut -f 1-6 > ${d}/cds.next
		# intersectBed -a ${d}/cds.next -b ${software_path}/genome/genome_region/5utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/5utr.txt  
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/5utr.txt ${d}/cds.next > ${d}/5utr.next
		# intersectBed -a ${d}/5utr.next -b ${software_path}/genome/genome_region/3utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/3utr.txt 
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/3utr.txt ${d}/5utr.next > ${d}/3utr.next
		# intersectBed -a ${d}/3utr.next -b ${software_path}/genome/genome_region/noncoding.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/noncoding.txt 
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/noncoding.txt ${d}/3utr.next > ${d}/noncoding.next
		# intersectBed -a ${d}/noncoding.next -b ${software_path}/genome/genome_region/intron.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/intron.txt
		# awk -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print $0"\tintergenic\t--\t--"}' ${d}/intron.txt ${d}/noncoding.next | sort -k1,4 -u > ${d}/intergenic.txt
		# cat ${d}/*.txt > ${d}.function.annot
		# rm -r ${d}
	fi
fi
if [ $StudyFile != "NULL" ]; then
	ls ${StudyFile} | while read i; do 
		# mkdir -p ${OUTDIR}/${i%.bed}
		type=study
		d=${OUTDIR}/${i%.bed}
		${software_path}/src/Get_intersectRegion2.py -g ${software_path}/genome/${Genome}_region -f $i -o ${d} -t $type -s $Strand
		# cut -f1-6 ${i} | intersectBed -a - -b ${software_path}/genome/genome_region/cds.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14 | sort -k1,4 -u > ${d}/cds.txt
		# cut -f1-6 ${i} | intersectBed -a - -b ${software_path}/genome/genome_region/cds.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14 | sort -k1,4 -u > ${OUTDIR}/${i%.bed}.study.annot
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/cds.txt ${i} | cut -f1-6 > ${d}/cds.next
		# intersectBed -a ${d}/cds.next -b ${software_path}/genome/genome_region/5utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/5utr.txt  
		# intersectBed -a ${d}/cds.next -b ${software_path}/genome/genome_region/5utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u >> ${OUTDIR}/${i%.bed}.study.annot 
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/5utr.txt ${d}/cds.next > ${d}/5utr.next
		# intersectBed -a ${d}/5utr.next -b ${software_path}/genome/genome_region/3utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/3utr.txt 
		# intersectBed -a ${d}/5utr.next -b ${software_path}/genome/genome_region/3utr.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u >> ${OUTDIR}/${i%.bed}.study.annot
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/3utr.txt ${d}/5utr.next > ${d}/3utr.next
		# intersectBed -a ${d}/3utr.next -b ${software_path}/genome/genome_region/noncoding.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/noncoding.txt 
		# intersectBed -a ${d}/3utr.next -b ${software_path}/genome/genome_region/noncoding.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u >> ${OUTDIR}/${i%.bed}.study.annot
		# awk 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print}' ${d}/noncoding.txt ${d}/3utr.next > ${d}/noncoding.next
		# intersectBed -a ${d}/noncoding.next -b ${software_path}/genome/genome_region/intron.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u > ${d}/intron.txt
		# intersectBed -a ${d}/noncoding.next -b ${software_path}/genome/genome_region/intron.txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u >> ${OUTDIR}/${i%.bed}.study.annot
		# awk -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print $0"\tintergenic\t--\t--"}' ${d}/intron.txt ${d}/noncoding.next| sort -k1,4 -u > ${d}/intergenic.txt
		# awk -v OFS="\t" 'NR==FNR{a[$1":"$2":"$3":"$4]++}NR!=FNR{if(!($1":"$2":"$3":"$4 in a)) print $0"\tintergenic\t--\t--"}' ${d}/intron.txt ${d}/noncoding.next| sort -k1,4 -u >> ${OUTDIR}/${i%.bed}.study.annot
		# cat ${d}/*.txt > ${OUTDIR}/${i%.bed}.study.annot
		# rm -r $d
	done
fi
