#!/bin/bash
helpdoc(){
    cat <<EOF
Description:
	software-build.sh pipeline.
Usage:
    software-build.sh -g $Genome --gtf $GTF
    -h/---help  -- help information
Module:         build
Options:
 -g          -- Prefix of output genome reference
 --gtf       -- gene annotation file (GTF format, recomannd GENOCODE annotation with "gene_name" and "gene_type" tags)
EOF
}

if [ $# = 0 ]; then
    helpdoc
    exit 1
fi

#get command line parameters
ARGS=`getopt -o g:h --long gtf:,help -n "$0" -- "$@"`
eval set -- "$ARGS"
while true ; do
	case "$1" in
		-g) Genome=$2 ; shift 2;;
		--gtf) GTF=$2 ; shift 2;;
		--iso) ISO=T ; shift 1;;
		-h) helpdoc ; exit 1;;
		--help) helpdoc ; exit 1;;
		--)
			shift
			break
			;;
		*) echo "unknown parameter:" $1; helpdoc; exit 1;;
	esac
done

#### build library
software_path=`which software`
software_path=${software_path%/software}
mkdir ${software_path}/genome/
# gtfToGenePred -genePredExt $GTF ${software_path}/genome/${Genome}.genePred
if [[ $ISO == "T" ]]; then
	gtfToGenePred -genePredExt -geneNameAsName2 $GTF ${software_path}/genome/${Genome}.genePred
	awk -v OFS="\t" '{$1=$1"\t"$12; print $0}' ${software_path}/genome/${Genome}.genePred | cut -f1-11 > ${software_path}/genome/${Genome}.txt
else
	gtfToGenePred $GTF ${software_path}/genome/${Genome}.genePred
	awk -v OFS="\t" '{$1=$1"\t1"; print $0}' ${software_path}/genome/${Genome}.genePred | cut -f1-11 > ${software_path}/genome/${Genome}.txt
fi
awk '{print $3"\t"$5"\t"$6"\t"$1":"$2"\t0\t"$4"\t"$10"\t"$11"\t"$7"\t"$8}' ${software_path}/genome/${Genome}.txt > ${software_path}/genome/${Genome}.temp

mkdir ${software_path}/genome/${Genome}_region
genome_region_dir="${software_path}/genome/${Genome}_region"

echo ''' awk '\''{if($7!=$8) {split($10, start, ","); split($11, end, ","); for(i=1; i< length(start); i++) {if(end[i] > $7 && end[i] < $8){ if(start[i]<=$7) {print $3"\t"$7"\t"end[i]"\tcds\t0\t"$4"\t"$2"\t"$1} else {print $3"\t"start[i]"\t"end[i]"\tcds\t0\t"$4"\t"$2"\t"$1}} else if(start[i] <$7 && end[i] >=$8){ print $3"\t"$7"\t"$8"\tcds\t0\t"$4"\t"$2"\t"$1} else if(end[i]>= $8 && start[i]< $8) print $3"\t"start[i]"\t"$8"\tcds\t0\t"$4"\t"$2"\t"$1} }}'\'' ''' ${software_path}/genome/${Genome}.txt "> ${genome_region_dir}/cds.txt" > script.sh

echo ''' awk '\''{split($10, start, ","); split($11, end, ","); for(i=1; i< (length(start)-1); i++) {print $3"\t"end[i]"\t"start[i+1]"\tintron\t0\t"$4"\t"$2"\t"$1}}'\'' ''' ${software_path}/genome/${Genome}.txt "> ${genome_region_dir}/intron.txt" >> script.sh

echo ''' awk '\''{if($7==$8) {split($10, start, ","); split($11, end, ","); for(i=1; i< length(start); i++) {print $3"\t"start[i]"\t"end[i]"\tnoncoding\t0\t"$4"\t"$2"\t"$1}} }'\'' ''' ${software_path}/genome/${Genome}.txt "> ${genome_region_dir}/noncoding.txt" >> script.sh

echo ''' awk '\''{if($7!=$8 && $4 == "+") {split($10, start, ","); split($11, end, ","); for(i=1; i < length(start); i++) {if(end[i] > $8 && start[i] <= $8) print $3"\t"$8"\t"end[i]"\t3utr\t0\t"$4"\t"$2"\t"$1; else if(end[i] > $8 && start[i] > $8) print $3"\t"start[i]"\t"end[i]"\t3utr\t0\t"$4"\t"$2"\t"$1}} else if($7!=$8 && $4 == "-") {split($10, start, ","); split($11, end, ","); for(i=1; i < length(start); i++) {if(start[i] < $7 && end[i] < $7) print $3"\t"start[i]"\t"end[i]"\t3utr\t0\t"$4"\t"$2"\t"$1; else if(start[i] < $7 && end[i] >=$7) print $3"\t"start[i]"\t"$7"\t3utr\t0\t"$4"\t"$2"\t"$1 }}}'\'' ''' ${software_path}/genome/${Genome}.txt "> ${genome_region_dir}/3utr.txt" >> script.sh 

echo ''' awk '\''{if($7!=$8 && $4 == "+") {split($10, start, ","); split($11, end, ","); for(i=1; i < length(start); i++) {if(start[i] < $7 && end[i] < $7) print $3"\t"start[i]"\t"end[i]"\t5utr\t0\t"$4"\t"$2"\t"$1; else if(start[i] < $7 && end[i] >=$7) print $3"\t"start[i]"\t"$7"\t5utr\t0\t"$4"\t"$2"\t"$1} } else if($7 !=$8 && $4 == "-") {split($10, start, ","); split($11, end, ",");  for(i=1; i < length(start); i++) {if(end[i] > $8 && start[i] <= $8) print $3"\t"$8"\t"end[i]"\t5utr\t0\t"$4"\t"$2"\t"$1; else if(end[i] > $8 && start[i] >$8) print $3"\t"start[i]"\t"end[i]"\t5utr\t0\t"$4"\t"$2"\t"$1}}}'\'' ''' ${software_path}/genome/${Genome}.txt "> ${genome_region_dir}/5utr.txt" >> script.sh

cat script.sh | parallel -j 5

 
