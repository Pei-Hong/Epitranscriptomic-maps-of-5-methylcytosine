#!/bin/bash
helpdoc(){
    cat <<EOF
Description:
	SCAPRUTE pipeline. A deep learning-embedded pipeline that captures polyadenylation information from 3 prime tag-based RNA-seq of single cells
Usage:
    software -m <build|Annot|Distance|Enrich> [options]
    -h/---help  -- help information
Module:         build
Options:
 -g          -- Prefix of output genome reference
 --gtf       -- gene annotation file (GTF format, recomannd GENOCODE annotation with "gene_name" and "gene_type" tags)
 
Module:         Annot
Options: -i studiedFile.bed -f FunctionalFile -g hg38 -o OUTDIR
 -i          -- bed6 format file for sites studied  
 -f          -- bed6 format file for regions with biological functions or a dictionary content bed6 format files.
 # --dir       -- Dictionary for functional region files, if set the "-f" would be ignored
 -g          -- Genome reference for ""
 -o          -- Directory for outputs, the default is "software_out"

Module:         Distance
Options: -i studiedFile.bed -f FunctionalFile -g hg38 –mature –o OUTDIR
 --annotFiles -- Directory for annotted files, which is the output for "Annot" module.
 --mature    -- If chose, the relative distance on mature RNA would be calculated, the default is precursor RNA.
 -g          -- Genome reference for ""
 -o          -- Directory for outputs, the default is "software_out"
 
Module:         Enrich 
Options: -i studiedFile.bed -f FunctionalFile -g hg38 --mature --model fisher --extend 50 --mid –o object --type around/left/right --start 0  
 --annotFiles -- Directory for annotted files, which is the output for "Annot" module, the default is "software_out"
 -g          -- Genome reference for "software"
 --fa        -- Fasta file for genome
 -o          -- Directory for outputs, the default is "software_out"
 --mature    -- If chose, the relative distance on mature RNA would be calculated, the default is precursor RNA.
 --model          -- Statistal model, "fisher" and "", the default is "fisher".
 -b          -- Background nucleotide, if choose "fisher" as the Statistal model, you should set one of "T|C|G|A" as the background nucleotide.
 --extend    -- Int value, extend the region, default is 0.
 --type      -- You can choose around/left/right to set the direction to extend, effective only when "extend" >= 0. The default is "around".
 --start     -- Int value, set up the start site of the extension.
Version: 1.0 2022/04/06
Author: Pei-Hong Zhang Email: zhangpeihong@picb.ac.cn
EOF
}

if [ $# = 0 ]; then
    helpdoc
    exit 1
fi

#default parameters:
MODULE="NULL" 
StudyFile="NULL"
FunctionFile="NULL" 
Genome="NULL" 
OUTDIR="software_out" 
# FunctionFile_Dir="NULL" 
BACKGROUND="N"
MATURE=0
MODEL="fisher"
EXTEND=0
# MID=0
TYPE="around" 
START=0
FASTA="NULL"
annotFiles="software_out"

#get command line parameters
ARGS=`getopt -o m:i:f:g:o:b:h --long gtf:,fa:,annotFiles:,mature,model:,extend:,type:,start:,help -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-m) MODULE=$2 ; shift 2;;
		-i) StudyFile=$2 ; shift 2;;
		-f) FunctionFile=$2 ; shift 2;;
		-g) Genome=$2 ; shift 2;;
		--annotFiles) annotFiles=$2 ; shift 2;;
		-o) OUTDIR=$2 ; shift 2;;
		--gtf) GTF=$2 ; shift 2;;
		--fa) FASTA=$2 ; shift 2;;
		# --dir) FunctionFile_Dir=$2 ; shift 2;;
		--mature) MATURE=1 ; shift ;;
		--model) MODEL=$2 ; shift 2;;
		-b) BACKGROUND=$2 ; shift 2;;
		--extend) EXTEND=$2 ; shift 2;;
		# --mid) MID=1 ; shift ;;
		--type) TYPE=$2 ; shift 2;;
		--start) START=$2 ; shift 2;;
		-h) helpdoc ; exit 1;;
		--help) helpdoc ; exit 1;;
		--)
			shift
			break
			;;
		*) echo "unknown parameter:" $1; helpdoc; exit 1;;
	esac
done

if [ $MODULE == "NULL" ];then
	helpdoc; exit 1
fi

#get path of software and package install 
if [ -z "$SCAPTUREPATH" ]; then
#	echo "scapture path is not set. Try to find scapture in current ENV"
	SCAPTUREPATH=$(which scapture)
	if [[ "$SCAPTUREPATH" == *scapture ]]; then
#		echo "scapture is found in: "$SCAPTUREPATH
		SCAPTUREPATH=${SCAPTUREPATH%scapture}
	else
		echo "scapture is not found!"
		exit 1
	fi
fi

#<build|Annot|Distance|Enrich>
if [ $MODULE == "build" ];then
	if [ $Genome != "NULL" && $GTF != "NULL" ];then
		software-build.sh -g $Genome --gtf $GTF
	else
		echo "The -gtf and -g options should be applied values."
		helpdoc; exit 1
	fi
fi

if [ $MODULE == "Annot" ]; then
	if [ $StudyFile != "NULL" && $FunctionFile != "NULL" && $Genome != "NULL" ];then
		software-Annot.sh -i $StudyFile -f $FunctionFile -g $Genome -o $OUTDIR
	else
		echo "The -i,-f and -g options should be applied values."
		helpdoc; exit 1
	fi
fi

if [ $MODULE == "Distance" ]; then
	if [ $StudyFile != "NULL" && $FunctionFile != "NULL" && $Genome != "NULL" ];then
		software-Distance.sh --annotFiles $annotFiles -g $Genome -o $OUTDIR
	else
		echo "The -i,-f and -g options should be applied values."
		helpdoc; exit 1
	fi
fi

if [ $MODULE == "Enrich" ]; then
	if [ $StudyFile != "NULL" && $FunctionFile != "NULL" && $Genome != "NULL" && $FASTA != "NULL" ];then
		software-Enrich.sh --annotFiles $annotFiles -g $Genome --fa $FASTA --mature $MATURE --model $MODEL -b $BACKGROUND --extend $EXTEND –o $OUTDIR --type $TYPE --start $START
	else
		echo "The -i,-f,--fa and -g options should be applied values."
		helpdoc; exit 1
	fi
fi
