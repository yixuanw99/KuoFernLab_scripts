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
input_dir="/home2/yxwu/findOG/OGfast/OG_narrow_all"
output_dir="/home2/yxwu/findOG/OGfast/OG_narrow_all_mafft"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  mafft --genafpair --maxiterate 1000 $file > "$output_dir"/${file_name%%.fna}_all_mafft.fas
done
