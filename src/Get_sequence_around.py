#!/picb/rnomics3/peihong/software/svr3/miniconda2/bin/python
#/Users/peihong 1/ly-svr5-rnomic4/0_script/Get_sequence_around.py
import os
import sys
import time
from optparse import OptionParser, SUPPRESS_HELP

usage = """
    Get sequence around certain sites in RNA.
    Get_sequence_around.py -s <sites.bed> 
    or
    Get_sequence_around.py -s <sites.bed> --regionlength  5utr 

input file format (without --regionlength):
chrome start   end     name     value   strand  exon    exon        site
chr21	8256780	8256933	5_8S_rRNA	0   +	8256780,	8256933,	8256933 
input file format (with --regionlength):
chrome start   end     name     value   strand  exon    exon        UTR    UTR    site
chr21	8256780	8256933	5_8S_rRNA	0   +	8256780,	8256933,	8256780  8256933    8256933 """
# input file format
#chrome start   end     name                      value strand  exon    exon        site
# chr21	8256780	8256933	5_8S_rRNA:ENST00000610460.1	0	+	8256780,	8256933,	8256933
# chr2	89600570	89600680	5S_rRNA:ENST00000612131.1	0	+	89600570,	89600680,	89600680
parser = OptionParser(usage=usage)
parser.add_option("--sites", "-s", dest="sites", help="A list of sites on RNA, bed formate")
parser.add_option("--out", "-o", dest="out_file", default="STD", help="a bed file output, ")
parser.add_option("--type", "-t", dest="sequence_type", default="around", help="Select the position of the sequence, 'right', 'left' or 'around', default=around. data 1 base")
parser.add_option("--seq_len", "-l", dest="sequence_length", type="int", default=50, help="Length of sequence you want, default=50")
parser.add_option("--seq_start", dest="sequence_start", type="int", default=0, help="Distance of the sequence you want to the sites, default=0")
parser.add_option("--out_type", dest="out_type", default="bed6", help="bed6 or bed12, default=bed6")
parser.add_option("--regionlength", dest="region", default="N", help="You want to get the sequence from regions or sites, defalut is 'N', you can set '5utr','3utr' or 'cds' or 'custom' and then '--seq_len' and '--seq_start' won't be used. And if you choose 'custom', please change both col9 and col10 as the site you want.")
parser.add_option("--if_half", dest="if_half", default="N", help="you can select 'Y' or 'N', and need cols which content sequence region, default='N")

(options, args) = parser.parse_args()

sites_file = options.sites
sequence_type = options.sequence_type

if sites_file == "-":
    sites_list = sys.stdin.readlines()
else:
    sites_list = open(sites_file).readlines()
seq_len = options.sequence_length
seq_start = options.sequence_start

def Get_Region_Len(region_start, region_end, start_list, end_list):
    start_list = [int(strand + i) for i in exon_starts.strip(",").split(",")]
    start_list.sort()
    end_list = [int(strand + i) for i in exon_ends.strip(",").split(",")]
    end_list.sort()
    if strand == "-":
        temp = start_list
        start_list = end_list
        end_list = temp
        # temp = region_start
        # region_start = region_end
        # region_end = temp
    start_exon = 0
    end_exon = 0
    # print(start_list)
    exon_length = [end_list[i] - start_list[i] for i in range(0, len(start_list))]
    # print(end_list[start_exon], int(region_start), end_list[start_exon] < int(region_start))
    while end_list[start_exon] < int(region_start) : # or end_list[start_exon] < int(region_start)
        start_exon = start_exon+1
    while end_list[end_exon] < int(region_end) : # or end_list[end_exon] < int(region_end)
        end_exon = end_exon+1
    # print(start_exon, end_exon)
    # print(region_start, region_end)
    
    if int(region_end) < int(region_start):
        return("5_extend") #### update 20211204
    elif int(region_end) < start_list[end_exon]:
        return("intron") #### update 20211204
    elif start_exon == end_exon: return(int(region_end) - int(region_start))
    elif start_exon != end_exon:
        sequence_len = end_list[start_exon] - int(region_start) + sum(exon_length[start_exon+1: end_exon]) + int(region_end) - start_list[end_exon]
        return(sequence_len)

record = []
if options.region == "N":
    for i in range(0, len(sites_list)):
        sites_list_i = sites_list[i]
        line = sites_list[i].strip("\n").split("\t")
        (RNA_start, RNA_end, name, strand, exon_starts, exon_ends, origsite) = line[1:4] + line[5:9]
        site = strand + origsite
        start_list = [int(strand + i) for i in exon_starts.strip(",").split(",")]
        start_list.sort()
        end_list = [int(strand + i) for i in exon_ends.strip(",").split(",")]
        end_list.sort()
        # print(line)
        if strand == "-":
            temp = start_list
            start_list = end_list
            end_list = temp

        ### half of sequence
        if options.if_half == 'Y':
            origsite=line[10]
            site = strand + origsite
            (utr5, utr3) =line[8:10]
            utr5 = strand + utr5
            utr3 = strand + utr3
            if strand == "-":   # need to exchange, because utr5 and utr3 are not appoint
                temp = utr5
                utr5 = utr3
                utr3 = temp
            start_exon = 0
            end_exon = 0
            while end_list[start_exon] < int(utr5) :
                start_exon = start_exon+1
            while end_list[end_exon] < int(utr3) : # or end_list[end_exon] < int(region_end)
                end_exon = end_exon+1
            # print(RNA_len)
            start_list = start_list[start_exon+1:end_exon+1]
            start_list.insert(0,int(utr5))
            end_list = end_list[start_exon:end_exon]
            end_list.append(int(utr3))
            # print(start_exon, end_exon,  start_list, end_list)
            gene_exon_length = [end_list[i] - start_list[i] for i in range(0,len(start_list))]
            RNA_len = sum(gene_exon_length)/2
            # print(RNA_len)
            temp_len = 0
            temp_i = 0
            while temp_len < RNA_len:
                temp_len = end_list[temp_i] - start_list[temp_i] + temp_len
                temp_i = temp_i +1
            cds_mid = end_list[temp_i-1] - (temp_len - RNA_len)
            # print(len(start_list), len(end_list), temp_len ,RNA_len, cds_mid, int(site))
            if int(site) < cds_mid:
                start_list = start_list[0:temp_i]
                end_list = end_list[0:(temp_i-1)]
                end_list.append(cds_mid)
                # print(1, start_list, end_list)
            elif int(site) > cds_mid:
                start_list = start_list[temp_i:len(start_list)]
                start_list.append(cds_mid)
                start_list.sort()
                end_list = end_list[(temp_i-1):len(end_list)]
                # print(2, start_list, end_list)

        seq_start_list = []
        seq_end_list = []
        if options.sequence_type == "left" or options.sequence_type == "around":
            i=0
            while i < len(start_list):
                if start_list[i] <= int(site) and end_list[i] >= int(site) :
                    start_length = seq_start
                    start_site = int(site)
                    while start_length >=0 and i >=0 :
                        # print(start_site,start_site - start_length, start_list[i])
                        if start_site - start_length < start_list[i]:
                            start_length = start_list[i] - (start_site - start_length)
                            # print(start_length)
                            start_site = end_list[i-1]
                            # print(start_site)
                            i = i-1
                        elif start_site - start_length == start_list[i]:
                            start_length = -1
                            start_site = end_list[i-1]
                            i = i-1
                            length = seq_len
                            while length >= 0 and i >=0:
                                # print(length)
                                # print(start_site - length , start_list[i])
                                if start_site - length < start_list[i]:
                                    seq_end_list.append(int(str(start_site).strip("-")))
                                    seq_start_list.append(int(str(start_list[i]).strip("-")))
                                    # print(start_list[i],start_site,length,start_list[i] - (start_site - length))
                                    length = start_list[i] - (start_site - length)
                                    start_site = end_list[i-1]
                                    i = i-1
                                elif length >0:
                                    seq_end_list.append(int(str(start_site).strip("-")))
                                    start_site = start_site - length
                                    length = -1
                                    seq_start_list.append(int(str(start_site).strip("-")))
                        else:
                            start_site = start_site - start_length
                            # seq_end_list.append(int(str(start_site).strip("-")))
                            start_length = -1
                            length = seq_len
                            while length >= 0 and i >=0:
                                # print(length)
                                # print(start_site - length , start_list[i])
                                if start_site - length < start_list[i]:
                                    seq_end_list.append(int(str(start_site).strip("-")))
                                    seq_start_list.append(int(str(start_list[i]).strip("-")))
                                    # print(start_list[i],start_site,length,start_list[i] - (start_site - length))
                                    length = start_list[i] - (start_site - length)
                                    start_site = end_list[i-1]
                                    i = i-1
                                elif length >0:
                                    seq_end_list.append(int(str(start_site).strip("-")))
                                    start_site = start_site - length
                                    length = -1
                                    seq_start_list.append(int(str(start_site).strip("-")))
                    
                    i=len(start_list)
                elif start_list[i] > int(site):
                    seq_start_list.append(int(str(int(site) - int(seq_len)).strip("-")))
                    seq_end_list.append(int(str(site).strip("-")))
                    i=len(start_list)
                else: i=i+1
        if options.sequence_type == "right" or options.sequence_type == "around":
            i=0
            while i < len(start_list):
                if start_list[i] <= int(site) and end_list[i] >= int(site):
                    start_length = seq_start
                    start_site = int(site)
                    while start_length >=0 and i < len(end_list):
                        if start_site + start_length > end_list[i]: 
                            start_length = start_site + start_length - end_list[i] ## new start_length != 0
                            if i+1 <len(end_list):
                                start_site = start_list[i+1]
                            i = i+1
                        elif start_site + start_length == end_list[i]:    ## new start_length == 0
                            start_length = -1
                            if i+1 <len(end_list):
                                start_site = start_list[i+1]
                            i = i+1
                            length = seq_len
                            while length >=0 and i <len(end_list):
                                if start_site + length > end_list[i]:
                                    seq_start_list.append(int(str(start_site).strip("-")))
                                    seq_end_list.append(int(str(end_list[i]).strip("-")))
                                    length = start_site + length - end_list[i]
                                    if(i+1 <len(end_list)):
                                        start_site = start_list[i+1]
                                    i = i+1
                                elif length >0:
                                    seq_start_list.append(int(str(start_site).strip("-")))
                                    start_site = start_site + length
                                    length = -1
                                    seq_end_list.append(int(str(start_site).strip("-")))
                        else:
                            start_site = start_site + start_length
                            start_length = -1
                            length = seq_len
                            while length >=0 and i <len(end_list):
                                if start_site + length > end_list[i]:
                                    seq_start_list.append(int(str(start_site).strip("-")))
                                    seq_end_list.append(int(str(end_list[i]).strip("-")))
                                    length = start_site + length - end_list[i]
                                    if(i+1 <len(end_list)):
                                        start_site = start_list[i+1]
                                    i = i+1
                                elif length >0:
                                    seq_start_list.append(int(str(start_site).strip("-")))
                                    start_site = start_site + length
                                    length = -1
                                    seq_end_list.append(int(str(start_site).strip("-")))
                    i=len(start_list)
                elif start_list[i] > int(site):
                    seq_start_list.append(int(str(site).strip("-")))
                    seq_end_list.append(int(str(int(site) + int(seq_len)).strip("-")))
                    i=len(start_list)
                else: i=i+1
        seq_start_list.sort()
        seq_end_list.sort()
        temp = seq_start_list
        # print(seq_start_list)
        # print(seq_end_list)
        seq_start_list = [i for i in seq_start_list if i not in seq_end_list]
        seq_end_list = [i for i in seq_end_list if i not in temp]
        # print(seq_end_list)
        if strand == "-":
            temp = seq_start_list
            seq_start_list = seq_end_list
            seq_end_list = temp
        if seq_start_list != []:
            # print(seq_start_list)
            seq_start_list[0] = seq_start_list[0] - 1
        exon_length = [str(seq_end_list[i] - seq_start_list[i]) for i in range(0, len(seq_start_list))]
        exon_start = [str(seq_start_list[i] - seq_start_list[0]) for i in range(0, len(seq_start_list))]
        exon_length.append("")
        exon_start.append("")
        seq_start_list.append("")
        seq_end_list.append("")
        seq_start_list = [str(i) for i in seq_start_list]
        seq_end_list = [str(i) for i in seq_end_list]
        # print(seq_end_list)
        if options.out_type == "bed12"and seq_end_list != [""]:
            record.append("\t".join([line[0], seq_start_list[0], seq_end_list[len(seq_end_list)-2], line[0] + ":"+ origsite, "0", strand, seq_start_list[0], seq_end_list[len(seq_end_list)-2],
            "0", str(len(seq_start_list)-1), 
            ",".join([str(i) for i in exon_length]), ",".join([str(i) for i in exon_start]), name]))
        elif options.out_type == "bed6" and seq_end_list != [""]:
            # print("bed6")
            for j in range(0, len(seq_start_list)-1):
                record.append("\t".join([line[0], seq_start_list[j], seq_end_list[j], line[0] + ":"+ origsite, "0", strand, name]))

else:
    for i in range(0, len(sites_list)):
        sites_list_i = sites_list[i]
        # print(sites_list_i)
        line = sites_list[i].strip("\n").split("\t")
        # print(line)
        (RNA_start, RNA_end, name, strand, exon_starts, exon_ends) = line[1:4] + line[5:8]
        start_list = [int(strand + i) for i in exon_starts.strip(",").split(",")]
        start_list.sort()
        end_list = [int(strand + i) for i in exon_ends.strip(",").split(",")]
        end_list.sort()
        if strand == "+":
            (utr5, utr3, origsite) = line[8:11]
        elif strand == "-": 
            # (utr3, utr5, origsite) = line[8:11] ## update 20211204
            (utr5, utr3, origsite) = line[8:11]
            # print(RNA_end)
            temp = RNA_start  ##
            RNA_start = RNA_end
            RNA_end = temp
            temp = utr5  ## update 20211204
            utr5 = utr3  ## update 20211204
            utr3 = temp  ## update 20211204
            temp = start_list  ##
            start_list = end_list
            end_list = temp
            
        if options.region == "5utr":  
            region_start = strand + RNA_start
            region_end = strand + utr5
            # print(region_start, region_end)
            total_len = Get_Region_Len(region_start, region_end, start_list, end_list)
            region_start = strand + RNA_start
            region_end = strand + origsite
            # print(region_start, region_end)
            site_len = Get_Region_Len(region_start, region_end, start_list, end_list)
        elif options.region == "cds":
            region_start = strand + utr5
            region_end = strand + utr3
            total_len = Get_Region_Len(region_start, region_end, start_list, end_list)
            region_start = strand + utr5
            region_end = strand + origsite
            site_len = Get_Region_Len(region_start, region_end, start_list, end_list)
        elif options.region == "3utr":
            region_start = strand + utr3
            region_end = strand + RNA_end
            total_len = Get_Region_Len(region_start, region_end, start_list, end_list)
            region_start = strand + utr3
            region_end = strand + origsite
            site_len = Get_Region_Len(region_start, region_end, start_list, end_list)
            # print(region_start, region_end, site_len)
        elif options.region == "custom":
            region_start = strand + RNA_start
            region_end = strand + utr5
            total_len = Get_Region_Len(region_start, region_end, start_list, end_list)
            region_start = strand + RNA_start
            region_end = strand + origsite
            site_len = Get_Region_Len(region_start, region_end, start_list, end_list)
        # total_len: length from the certain site you set as a region to mRNA start site;
        # site_len: length from the site you want to know to mRNA start site;
        # and site_len - total_len is the relative distance of the site to a certain site in mRNA, that is awk '{print $13-$12}' .
        if site_len in ["intron", "5_extend"]: 1
        elif site_len > total_len and options.region != "custom":
            site_len = "3_extend"
        record.append("\t".join([sites_list_i.strip("\n"), str(total_len), str(site_len)]))



record.append("")
if options.out_file == "STD" :
    sys.stdout.write("\n".join(record))
else:
    out = open(options.out_file,"w")
    out.write("\n".join(record))
    out.close()



        # start_list = [int(strand + i) for i in exon_starts.strip(",").split(",")]
        # start_list.sort()
        # end_list = [int(strand + i) for i in exon_ends.strip(",").split(",")]
        # end_list.sort()
        # if strand == "-":
        #     temp = start_list
        #     start_list = end_list
        #     end_list = temp
        #     # temp = region_start
        #     # region_start = region_end
        #     # region_end = temp
        # start_exon = 0
        # end_exon = 0
        # # print(start_list)
        # exon_length = [end_list[i] - start_list[i] for i in range(0, len(start_list))]
        # # print(end_list[start_exon], int(region_start), end_list[start_exon] < int(region_start))
        # while end_list[start_exon] < int(region_start) : # or end_list[start_exon] < int(region_start)
        #     start_exon = start_exon+1
        # while end_list[end_exon] < int(region_end) : # or end_list[end_exon] < int(region_end)
        #     end_exon = end_exon+1
        # # print(start_exon, end_exon)
        # if start_exon == end_exon: record.append("\t".join([sites_list_i.strip("\n"), str(int(region_end) - int(region_start))]))
        # elif start_exon != end_exon:
        #     sequence_len = end_list[start_exon] - int(region_start) + sum(exon_length[start_exon+1: end_exon]) + int(region_end) - start_list[end_exon]
        #     record.append("\t".join([sites_list_i.strip("\n"), str(sequence_len)]))