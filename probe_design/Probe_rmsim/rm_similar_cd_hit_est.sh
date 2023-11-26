#! /bin/bash

#################################################################
# Function :remove similar sequences by cd-hit-est clustering   #
# Platform :All Linux Based Platform                            #
# Version  :1.0                                                 #
# Date     :2023-03-16                                          #
# Author   :Yihsuan Wu                                          #
# Contact  :yixuan.w99@gmail.com                                #
#################################################################

#-------------------
input_dir="/home2/yxwu/codon12_compare/output_eqone_MACSE/1_perfect_all_rmKAI"
output_dir="/home2/yxwu/test_rmsim/1_perfect_all_rmKAI_rmSimilar"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  echo ${file_name}
  sed 's/-/N/g' ${file} > "${output_dir}/temp"${file_name}
  cd-hit-est -i "${output_dir}/temp"${file_name} -o "${output_dir}/"${file_name} -c 0.95 -n 10 -d 0 -gap 0 -gap-ext 0
  sed -i 's/N/-/g' "${output_dir}/"${file_name}
  rm "${output_dir}/temp"${file_name}
done
