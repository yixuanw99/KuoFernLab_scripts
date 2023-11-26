import os


import_align_dir = '/home2/yxwu/cut_fas/output_eqone_MACSE/fasta_list5/'
output_dir = '/home2/yxwu/test_codon12_compare/output_eqone_MACSE/fasta_list5/'
wrong1_dir = output_dir + "wrong1/"
sample_check_dir = output_dir + "sample_check/"
perfect_dir = output_dir + "perfect/"
wrong2_dir = output_dir + "wrong2/"
alignment_err_dir = output_dir + "alignment_err/"
frame_check_dir = output_dir + "frame_check/"

def get_ref_ID(ref_str,ID_list):
    for header in ID_list:
        if header[0:3] == ref_str: #choosing reference seq for boundaries
            return header

def basecomparing(str1,str2):
    if ((str1 == "-") and (str2 == "-")):
        return float(1)
    elif ((str1 == "-") and (str2 != "-")):
        return float(0.2)
    elif ((str1 != "-") and (str2 == "-")):
        return float(0.2)
    elif (str1 == str2):
        return float(1)
    elif ((str1 != "-") and (str2 != "-") and (str1 != str2)):
        return float(0)
    else:
        return "error"

def comparing_two_sequence(sequence1,sequence2):
    point = 0.0
    if len(sequence1) != len(sequence2):
        return "error: sequence length not equal"
    else:
        for nuc_index in range(len(sequence1)):
            point = point + basecomparing(sequence1[nuc_index],sequence2[nuc_index])
        return point

def get_degap_length_infloat(sequence_str):
    return float(len(str(sequence_str.replace('-',''))))

def is_list_equal(listName):
    if len(set(listName)) == 1:
        return True
    elif len(set(listName)) > 1:
        return False
    else:
        return "error"

def all_gt_60(ls):
    for i in ls:
        if i < 0.6:
            return False
    return True

def all_lt_60(ls):
    for i in ls:
        if i > 0.6:
            return False
    return True

def similarity_60_determine(ls):
    if all_gt_60(ls) == all_lt_60(ls):
        return 1
    elif ((all_gt_60(ls) == True) and (all_lt_60(ls) == False)):
        return 2
    elif ((all_gt_60(ls) == False) and (all_lt_60(ls) == True)):
        return 0
    else:
        return "error"


fasta_list = []

for file in os.listdir(import_align_dir):
    name_length = len(file.split('.'))
    if (file.split('.')[name_length-1] == 'fas'): #file type filelsit and taxalist
        fasta_list.append(file)

print(fasta_list)

for fasta in fasta_list:
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

    frame_1 = []
    for sample in fa_ID:
        frame_1_temp = ""
        for pos1 in range(0,len(fastaID_seq_dic[sample]),3):
            if (pos1+1) < len(fastaID_seq_dic[sample]):
                frame_1_temp += (str(fastaID_seq_dic[sample])[pos1] + str(fastaID_seq_dic[sample])[pos1+1])
        frame_1.append(frame_1_temp)
    fastaID_frame1_dict = dict(zip(fa_ID,frame_1))

    frame_2 = []
    for sample in fa_ID:
        frame_2_temp = ""
        for pos1 in range(1,len(fastaID_seq_dic[sample]),3):
            if (pos1+1) < len(fastaID_seq_dic[sample]):
                frame_2_temp += (str(fastaID_seq_dic[sample])[pos1] + str(fastaID_seq_dic[sample])[pos1+1])
        frame_2.append(frame_2_temp)
    fastaID_frame2_dict = dict(zip(fa_ID,frame_2))

    frame_3 = []
    for sample in fa_ID:
        frame_3_temp = ""
        for pos1 in range(2,len(fastaID_seq_dic[sample]),3):
            if (pos1+1) < len(fastaID_seq_dic[sample]):
                frame_3_temp += (str(fastaID_seq_dic[sample])[pos1] + str(fastaID_seq_dic[sample])[pos1+1])
        frame_3.append(frame_3_temp)
    fastaID_frame3_dict = dict(zip(fa_ID,frame_3))

    fastaID_allframe_dictlist = [fastaID_frame1_dict,fastaID_frame2_dict,fastaID_frame3_dict]

    print("~~~~~~~~",fasta,"~~~~~~~~")
    max_val_list = []
    max_posision_list = []
    for sample in fa_ID:
        point = []
        similarity = []
        readingFrame = 0
        for codon12_dict in fastaID_allframe_dictlist:
            ref_header = get_ref_ID('KAI',fa_ID)
            if sample[0:3] == 'KAI':
                continue
            point.append(comparing_two_sequence(codon12_dict[ref_header],codon12_dict[sample]))
            similarity.append(point[readingFrame]/float(len((codon12_dict[ref_header]))))
            readingFrame += 1
        if sample[0:3] == 'KAI':
                continue
        max_val = max(similarity)
        max_val_list.append(max_val)
        max_posision = similarity.index(max_val)
        max_posision_list.append(max_posision)
    
    with open(fasta_file, "r") as input:
        if (is_list_equal(max_posision_list) == True):
            # print("good frame")
            if (similarity_60_determine(max_val_list)) == 0:
                print("wrong1")
                output_alignment_name = wrong1_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            elif (similarity_60_determine(max_val_list)) == 1:
                print("sample_check")
                output_alignment_name = sample_check_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            elif (similarity_60_determine(max_val_list)) == 2:
                print("perfect")
                output_alignment_name = perfect_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            else:
                print("error")
        elif (is_list_equal(max_posision_list) == False):
            # print("bad frame")
            if (similarity_60_determine(max_val_list)) == 0:
                print("wrong2")
                output_alignment_name = wrong2_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            elif (similarity_60_determine(max_val_list)) == 1:
                print("alignment_err")
                output_alignment_name = alignment_err_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            elif (similarity_60_determine(max_val_list)) == 2:
                print("frame_check")
                output_alignment_name = frame_check_dir + "final_" + fasta
                with open(output_alignment_name, "w") as output_alignment:
                    for line in input:
                        output_alignment.write(line)
                output_alignment.close
            else:
                print("error")
        else:
            print("error")




# for fasta in fasta_list:
#     print("--------",fasta,"---------")
#     for sample in fa_ID:
#         point = 0.0
#         similarity = 0.0
#         readingFrame = 0
#         for dict in fastaID_allframe_dictlist:
#             ref_header = get_ref_ID('KAI',fa_ID)
#             if sample[0:3] == 'KAI':
#                 continue
#             point = comparing_two_sequence(dict[ref_header],dict[sample])
#             similarity = point/get_degap_length_infloat(dict[ref_header])
#             readingFrame += 1
#             print(readingFrame,sample,"similarity",similarity)
