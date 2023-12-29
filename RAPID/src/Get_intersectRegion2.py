#!/bin/python

import os
import sys
import commands
from optparse import OptionParser, SUPPRESS_HELP
import pybedtools

usage = '''
Get_intersectRegion.py -g path/to/genome/genome_region -f $i -o Dir -t type
'''

parser = OptionParser(usage=usage)
parser.add_option("-g", dest="genome", help = "")
parser.add_option("-f", dest="inputfile",help = "")
parser.add_option("-o", dest="outDir", help = "")
parser.add_option("-t", dest="type", help = "")
parser.add_option("-s", dest="strand", help = "")

(options, args) = parser.parse_args()
input_file = options.inputfile
genome_region_path = options.genome
outDir = options.outDir
TYPE = options.type

def saveFile(outfile, data):
	out = open(outfile, 'w')
	out.write(''.join([data, '\n']))
	out.close()

region_order = ["3utr","noncoding", "5utr","intron","cds","intergenic"]
results = dict()
bed_file = pybedtools.BedTool(commands.getoutput(''.join(["cut -f1-6 ", input_file])), from_string=True).sort()
bed_file_rest = bed_file
for i in region_order[0:(len(region_order)-1)]:
	# print(i)
	region = pybedtools.BedTool(''.join([genome_region_path,"/", i ,".txt"]))
	if options.strand == "F":
		tmp = list(bed_file_rest.intersect(region, wa=True, wb=True, f=0.5))
	else:
		tmp = list(bed_file_rest.intersect(region, wa=True, wb=True, s=True, f=0.5))
	tmp = [[tmp[y][x] for x in [0,1,2,3,4,5,9,12,13]] for y in range(0,len(tmp))]
	results[i] = list()
	peak_dic = dict()
	for x in range(0,len(tmp)):
		if '\t'.join(tmp[x][0:4]) not in peak_dic:
			peak_dic['\t'.join(tmp[x][0:4])] = tmp[x]
	for key in peak_dic:
		results[i].append(peak_dic[key])

	# [results[i].append(x) for x in tmp if not x in results[i]]
	results[i] = pybedtools.BedTool(results[i]).sort()
	if options.strand == "F":
		bed_file_rest = bed_file_rest.intersect(region, v=True, f=0.5).sort()
	else:
		bed_file_rest = bed_file_rest.intersect(region, v=True, f=0.5, s=True).sort()

results["intergenic"] = list(bed_file_rest)
results["intergenic"] = [results["intergenic"][x][:] + ["intergenic", "--", "--"] for x in range(0,len(results["intergenic"]))]

final_out = list()
for key in region_order:
	final_out = final_out + list(results[key])
	# results[key] = ['\t'.join(results[key][x]) for x in range(0,len(results[key]))]
	# final_out.append('\n'.join(results[key]))
t = pybedtools.BedTool(final_out).saveas(''.join([outDir, ".", TYPE, ".annot"]))

# final_out.append("")

# final_out_file = open(''.join([outDir, ".function.annot"]), 'w')
# final_out_file.write('\n'.join(final_out)) 
# final_out_file.close()
