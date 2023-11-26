#! /bin/bash

##################################################################
# Function :generate a table containing exon and intron length   #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-02-15                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_ref_dir="/home2/yxwu/findOG/OGfast/OG_eqone_Lomariopsis_KAI_MACSE_sedgap"
input_seq_dir="/home2/yxwu/findOG/OGfast/OG_eqone_KAI_exon_split/fas_split"
output1_dir="/home2/yxwu/findOG/OGfast/OG_eqonewithintron_all_MACSE"
# output2_dir="/home2/yxwu/findOG/OGfast/OG_all_mafft"
#-------------------

for file in "$input_seq_dir"/*.fas
do
	file_name=${file##*/}
	OG=${file_name%%_*}
	echo "=> OG: "$OG
	echo "處理檔案：${file_name}"
	echo "其絕對路徑為: ${file}"
	# echo "${file_name%%.fas}"
	# mafft --auto "${input_ref_dir}/${OG}.fna" > "${output1_dir}/${OG}.fas"
	# mafft --auto --add "${input_seq_dir}/${OG}.fna" "${output1_dir}/${OG}.fas" > "${output2_dir}/${OG}_added.fas"
	mafft --localpair --maxiterate 1000 --addfragments "${file}" "${input_ref_dir}"/"${OG}"_MACSE.fas > "${output1_dir}"/"${file_name%%.fas}"_MACSE.fas
done

