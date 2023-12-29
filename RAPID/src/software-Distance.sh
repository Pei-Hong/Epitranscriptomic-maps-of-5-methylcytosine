#!/bin/bash
helpdoc(){
    cat <<EOF
Description:
	SCAPRUTE pipeline. A deep learning-embedded pipeline that captures polyadenylation information from 3 prime tag-based RNA-seq of single cells
Usage:
    software -m <Build|Annot|Distance|Enrich> [options]
    -h/---help  -- help information

Module:         Distance
Options: -i studiedFile.bed -f FunctionalFile -g hg38 –mature –o OUTDIR
 --annotFiles -- The output for "Annot" module, which could be one file or
                a directory content the files.
 -g          -- Genome reference for "software".
 -o          -- Directory for outputs, the default is "software_dis".
 --start_site -- Distance is calculated between a start site and your studied 
                sites. You can specify the start site of a certain region 
                ('5utr', 'cds' and '3utr'). And can also specify a 'custom' 
                site as the start site. If you choose 'custom', please add 
                the sites on col10 of file.annot.

Version: 1.0 2022/04/06
Author: Pei-Hong Zhang Email: zhangpeihong@picb.ac.cn
EOF
}

if [ $# = 0 ]; then
    helpdoc
    exit 1
fi

#default parameters:
Genome="NULL" 
OUTDIR="NULL" 
annotFiles="NULL"
start_site="NULL"

#get command line parameters
ARGS=`getopt -o g:o: --long annotFiles:,start_site:,help -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-g) Genome=$2 ; shift 2;;
		--annotFiles) annotFiles=$2 ; shift 2;;
		-o) OUTDIR=$2 ; shift 2;;
		--start_site) start_site=$2 ; shift 2;;
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
export PATH=${software_path%/software}/src:$PATH

mkdir -p $OUTDIR
if [ $start_site != "custom" ]; then
	if [ -f $annotFiles ]; then
		out_file=${annotFiles##*/}
		awk -v r=$start_site 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a && $7 == r) print a[$9":"$8]"\t"$3"\t"$7}' ${software_path}/genome/${Genome}.temp $annotFiles | Get_sequence_around.py -s - --regionlength $start_site > ${OUTDIR}/${out_file%.annot}.${start_site}.dis
	else
		for file in ${annotFiles}/*.annot; do
			out_file=${file##*/}
			awk -v r=$start_site 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a && $7 == r) print a[$9":"$8]"\t"$3"\t"$7}' ${software_path}/genome/${Genome}.temp $file | Get_sequence_around.py -s - --regionlength $start_site > ${OUTDIR}/${out_file%.annot}.${start_site}.dis
		done
	fi
else
	if [ -f $annotFiles ]; then
		out_file=${annotFiles##*/}
		awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$3"\t"$7"\t"$10}' ${software_path}/genome/${Genome}.temp $annotFiles | awk -v OFS="\t" '{$9=$13; $10=$13; print}' | cut -f 1-12| Get_sequence_around.py -s - --regionlength $start_site > ${OUTDIR}/${out_file%.annot}.${start_site}.dis
	else
		for file in ${annotFiles}/*.annot; do
			out_file=${file##*/}
			awk 'NR==FNR{a[$4]=$0}NR!=FNR{if($9":"$8 in a) print a[$9":"$8]"\t"$3"\t"$7"\t"$10}' ${software_path}/genome/${Genome}.temp $file | awk -v OFS="\t" '{$9=$13; $10=$13; print}' | cut -f 1-12| Get_sequence_around.py -s - --regionlength $start_site > ${OUTDIR}/${out_file%.annot}.${start_site}.dis
		done
	fi
fi

