#!/bin/python

import os
import sys
import commands
from optparse import OptionParser, SUPPRESS_HELP
usage = '''
Get_intersectRegion.py -g path/to/genome/genome_region -f $i -o Dir
'''

parser = OptionParser(usage=usage)
parser.add_option("-g", dest="genome", help = "")
parser.add_option("-f", dest="inputfile",help = "")
parser.add_option("-o", dest="outDir", help = "")

(options, args) = parser.parse_args()
input_file = options.inputfile
genome_region_path = options.genome
outDir = options.outDir

def saveFile(outfile, data):
	out = open(outfile, 'w')
	out.write(''.join([data, '\n']))
	out.close()

region_order = ["cds","5utr","3utr","noncoding","intron","intergenic"]
results = dict()
for i in region_order[0:(len(region_order)-1)]:
	# print(i)
	if region_order.index(i) == 0:
		cmd=["cut -f1-6 ", input_file, " | intersectBed -a - -b ", genome_region_path,"/", i ,".txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14 | sort -k1,4 -u "]

		results[i] = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", i ,".txt"]), results[i])
		cmd=["awk 'NR==FNR{a[$1\":\"$2\":\"$3\":\"$4]++}NR!=FNR{if(!($1\":\"$2\":\"$3\":\"$4 in a)) print}' ", outDir, "/", i ,".txt ", input_file, "| cut -f 1-6"]
		result_next = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", i ,".next"]), result_next)
		file_next=i
	elif region_order.index(i) == (len(region_order) - 2):
		cmd=["intersectBed -a ", outDir, "/", file_next, ".next -b ", genome_region_path, "/", i ,".txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u "]
		results[i] = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", i ,".txt"]), results[i])

		cmd=["awk -v OFS=\"\t\" 'NR==FNR{a[$1\":\"$2\":\"$3\":\"$4]++}NR!=FNR{if(!($1\":\"$2\":\"$3\":\"$4 in a)) print $0\"\tintergenic\t--\t--\"}' ", outDir, "/", i ,".txt ", outDir, "/", file_next, ".next | sort -k1,4 -u "]
		results[region_order[len(region_order)-1]] = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", region_order[len(region_order)-1] ,".txt"]), results[region_order[len(region_order)-1]])
	else:
		cmd=["intersectBed -a ", outDir, "/", file_next, ".next -b ", genome_region_path, "/", i ,".txt -wa -wb -s -f 0.5| cut -f 1-6,10,13-14| sort -k1,4 -u "]
		results[i] = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", i ,".txt"]), results[i])

		cmd=["awk 'NR==FNR{a[$1\":\"$2\":\"$3\":\"$4]++}NR!=FNR{if(!($1\":\"$2\":\"$3\":\"$4 in a)) print}' ", outDir, "/", i ,".txt ", outDir, "/", file_next, ".next "]
		result_next = commands.getoutput(''.join(cmd)).strip('\n')
		saveFile(''.join([outDir, "/", i ,".next"]), result_next)
		file_next=i

final_out = list()
for key in region_order:
	final_out.append(results[key])
final_out.append("\n")

final_out_file = open(''.join([outDir, ".function.annot"]), 'w')
final_out_file.write('\n'.join(final_out)) 
final_out_file.close()

