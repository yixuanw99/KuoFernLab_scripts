#! /bin/bash

##################################################################
# Function :align OG sequence and add Lom sequence into alignment#
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-02-15                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_ref_dir="/home2/yxwu/findOG/OGfast/OG_ref"
input_seq_dir="/home2/yxwu/findOG/OGfast/OG_Lomariopsis"
output1_dir="/home2/yxwu/findOG/OGfast/OG_ref_mafft"
output2_dir="/home2/yxwu/findOG/OGfast/OG_all_mafft"
#-------------------

while read OG_list
do
	echo "OG: ${OG_list}"
	mafft --auto "${input_ref_dir}/${OG_list}.fna" > "${output1_dir}/${OG_list}.fas"
	mafft --auto --add "${input_seq_dir}/${OG_list}.fna" "${output1_dir}/${OG_list}.fas" > "${output2_dir}/${OG_list}_added.fas"
done < OG_list.txt
