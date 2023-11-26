#! /bin/bash

#################################################################
# Function :check whether each fas containing all Lom sequence  #
# Platform :All Linux Based Platform                            #
# Version  :1.0                                                 #
# Date     :2023-03-03                                          #
# Author   :Yihsuan Wu                                          #
# Contact  :yixuan.w99@gmail.com                                #
#################################################################

#-------------------
input_dir="/home2/yxwu/test_codon12_compare/perfect_eqone_all"
output_dir="/home2/yxwu/OrthoFinder/20221230run/probe_sequence"
# input_dir="/home/yixuan/OrthoFinder/check_Lom_in_dir/input"
# output_dir="/home/yixuan/OrthoFinder/check_Lom_in_dir/output"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  echo ${file_name}
  grep -q "K014639" ${file}
  check_K014639=$?
  grep -q "KTHU1995" ${file}
  check_KTHU1995=$?
  grep -q "RS_27" ${file}
  check_RS_27=$?
  if [[ $check_K014639 -eq 0 ]] && [[ $check_KTHU1995 -eq 0 ]] && [[ $check_RS_27 -eq 0 ]]; then
    echo "all pass"
	sed '/KAI/,+1 d' ${file} > ${output_dir}/${file_name}
  else
    echo "not all"
  fi
done
