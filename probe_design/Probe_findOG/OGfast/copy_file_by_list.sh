#! /bin/bash

###########################################
# Function :copy file by list             #
# Platform :All Linux Based Platform      #
# Version  :1.0                           #
# Date     :2023-02-24                    #
# Author   :Yihsuan Wu                    #
# Contact  :yixuan.w99@gmail.com          #
###########################################

#-------------------
list_dir="/home2/yxwu/findOG/OGfast/eqone_OGlist.txt"
input_dir="/home2/yxwu/findOG/OGfast/OG_Lomariopsis_KAI"
output_dir="/home2/yxwu/findOG/OGfast/OG_eqone_Lomariopsis_KAI"
#-------------------

while read line
do
  echo "copying ${input_dir}/${line} to ${output_dir}"
  cp "${input_dir}"/${line}.fna "${output_dir}"/
done < "$list_dir"