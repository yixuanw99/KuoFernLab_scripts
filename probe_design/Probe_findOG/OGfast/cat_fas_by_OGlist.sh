#! /bin/bash

##################################################################
# Function :cat fasta file in particular directory               #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-02-15                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_dir1="/home2/yxwu/findOG/OGfast/OG_Lomariopsis"
input_dir2="/home2/yxwu/findOG/OGfast/OG_KAI_intronadd"
output_dir="/home2/yxwu/findOG/OGfast/OG_narrow_all"
#-------------------

# for file in "$input_dir1"/*
while read OG_name
do
  file_name=${OG_name}.fna
  echo $file_name
  if [[ -f "$input_dir1"/"$file_name" ]] && [[ -f "$input_dir2"/"$file_name" ]]
  then
    echo "yooo"
	cat "$input_dir1"/"$file_name" "$input_dir2"/"$file_name" > "$output_dir"/${file_name%%.fna}.fas
  else
    echo "files not exist in one dir"
  fi
done < narrow_OGlist.txt
