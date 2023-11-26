import os
import collections

import_align_dir = '/home2/yxwu/findOG/OGfast/OG_eqonewithintron_all_MACSE_removeKAIcds/'
output_dir = '/home2/yxwu/cut_fas/output_eqone_MACSE/'

boundary_table = open("loci_boundary.csv", "w")
boundary_table.write("sequence" + "\t" + "alignment_start_site" + "\t" + "alignment_end_site" + "\t" + "degap_length" + "\n")

alignment_sum = open("alignment_sum.csv", "w")
alignment_sum.write("alignment" + "\t" + "trim_alignment_length" + "\t" + "max_ref_legnth" + "\n")

fasta_list = []

for file in os.listdir(import_align_dir):
    if (' ' in file):
        continue
    name_length = len(file.split('.'))
    if (file.split('.')[name_length-1] == 'fas'): #file type filelsit and taxalist
        left_intronlen = (file.split('_')[4])
        right_intronlen = (file.split('_')[5])
        if (((int(left_intronlen[4:])>0) and (int(left_intronlen[4:])<600)) or ((int(right_intronlen[5:])>0) and (int(right_intronlen[5:])<600))):
            fasta_list.append(file)

print(fasta_list)

for fasta in fasta_list:
    output_alignment_name = output_dir + "trim_" + fasta
    output_alignment = open(output_alignment_name, "w")
    fasta_file = import_align_dir + fasta
    fa_ID = []
    fa_Seq = []
    fa_Num = -1
    for line in open(fasta_file,"r").readlines():
        line = line.rstrip()
        if line.startswith('>'):
            fa_ID.append(line.replace('>',''))
            fa_Num = fa_Num + 1
            fa_Seq.append("")
        else:
            fa_Seq[fa_Num] = fa_Seq[fa_Num] + line
    fastaID_seq_dic = dict(zip(fa_ID,fa_Seq))
    ref_list = []
    start_site_list = []
    end_site_list = []
    ref_degap_length_list = []
    for sample in fa_ID:
        if sample[0:3] == 'KAI': #choosing reference seq for boundaries
            start = 0
            end = len(str(fastaID_seq_dic[sample]))-1
            while(str(fastaID_seq_dic[sample])[start] == '-'):
                start = start + 1
            while(str(fastaID_seq_dic[sample])[end] == '-'):
                end = end - 1
            start_site = start + 1
            end_site = end + 1
            ref_list.append(sample)
            start_site_list.append(start_site)
            end_site_list.append(end_site)
            ref_degap_length_list.append(len(str(fastaID_seq_dic[sample].replace('-',''))))
            #print(sample + " " + str(start_site) + " to " + str(end_site))
    aln_left = min(start_site_list)
    aln_right = max(end_site_list)
    seq_max_legnth = max(ref_degap_length_list)
    aln_size = aln_right - aln_left + 1
    alignment_sum.write(fasta + "\t" + str(aln_size) + "\t" + str(seq_max_legnth) + "\n")
    print(fasta + " boundary " + str(aln_left) + " to " + str(aln_right))
    for ID in sorted(fa_ID):
        extract_seq = str(fastaID_seq_dic[ID])[aln_left-1:aln_right]
        real_len = len(extract_seq.replace('-',''))
        if ID[0:3] == 'KAI':
            if len(extract_seq) > 79:
                output_alignment.write(">" + ID + "_" + "size" + str(real_len) + "\n" + extract_seq + "\n")
        else:
            if real_len > 79: # control for real length
                output_alignment.write(">" + ID + "_" + "size" + str(real_len) + "\n" + extract_seq + "\n")
            #print("output " + ID)
    output_alignment.close
    
    item = 1
    while(item <= len(ref_list)):
        E = item -1
        boundary_table.write(ref_list[E] + "\t" + str(start_site_list[E]) + "\t" + str(end_site_list[E]) + "\t" + str(ref_degap_length_list[E]) + "\n")
        item = item + 1

for fasta in fasta_list:
    output_alignment_name = output_dir + "trim_" + fasta
    linecount = len(open(output_alignment_name, 'rb').readlines())
    if (linecount <= 2 ):
        output_alignment = open(output_alignment_name, "w")
        output_alignment.write("")
        output_alignment.close

alignment_sum.close
boundary_table.close


