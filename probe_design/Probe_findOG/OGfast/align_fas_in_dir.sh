#! /bin/bash

##################################################################
# Function :align_fas_in_dir                                     #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-02-15                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_dir="/home2/yxwu/findOG/OGfast/OG_Lomariopsis"
output_dir="/home2/yxwu/findOG/OGfast/OG_Lomariopsis_mafft"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  mafft --genafpair --maxiterate 1000 $file > "$output_dir"/${file_name%%.fna}_mafft.fas
done
