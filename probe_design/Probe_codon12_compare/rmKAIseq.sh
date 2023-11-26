#! /bin/bash

#################################################################
# Function :rm KAI sequence in each fasta file                  #
# Platform :All Linux Based Platform                            #
# Version  :1.0                                                 #
# Date     :2023-03-16                                          #
# Author   :Yihsuan Wu                                          #
# Contact  :yixuan.w99@gmail.com                                #
#################################################################

#-------------------
input_dir="/home2/yxwu/Probe_codon12_compare/old_perfect_eqone_MACSE"
output_dir="/home2/yxwu/Probe_codon12_compare/old_perfect_eqone_MACSE_rmKAI"
# input_dir="/home/yixuan/OrthoFinder/check_Lom_in_dir/input"
# output_dir="/home/yixuan/OrthoFinder/check_Lom_in_dir/output"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  echo "processing: "${file_name}
  sed '/KAI/,+1 d' ${file} > ${output_dir}/${file_name}
done
