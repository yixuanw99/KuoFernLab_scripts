#! /bin/bash

####################################################################
# Function :transfer fasta sequence in multiple line into one line #
# Platform :All Linux Based Platform                               #
# Version  :1.0                                                    #
# Date     :2023-01-30                                             #
# Author   :Yihsuan Wu                                             #
# Contact  :yixuan.w99@gmail.com                                   #
####################################################################

#-------------------
input_dir="/home2/yxwu/test_fas"
output_dir="/home2/yxwu/test_fas/output"
#-------------------

#副檔名是fasta, fas, fa, fna, faa才能做
for file in "$input_dir"/*
do
  if [[ ${file##*.} == "fasta" ]] || [[ ${file##*.} == "fas" ]] || [[ ${file##*.} == "fa" ]] || [[ ${file##*.} == "fna" ]] || [[ ${file##*.} == "faa" ]];then
    file_name=${file##*/}
	echo "處理檔案：${file##*/}"
	echo "其絕對路徑為: ${file}"
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${file} > "${output_dir}"/${file_name%.*}_temp.fas
	tail -n +2 "${output_dir}"/${file_name%.*}_temp.fas > "${output_dir}"/${file_name%.*}_oneline.fas
	rm "${output_dir}"/${file_name%.*}_temp.fas
  else
    echo ${file}"不是檔案"
  fi
done
