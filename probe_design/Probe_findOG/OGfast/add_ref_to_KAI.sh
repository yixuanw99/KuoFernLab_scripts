#! /bin/bash

##################################################################
# Function :add "ref" to KAI header                              #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-03-06                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_dir="/home2/yxwu/findOG/OGfast/OG_eqone_Lomariopsis_KAI"
# output_dir="/home2/yxwu/findOG/OGfast/OG_Lomariopsis_mafft"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  sed -i 's/KAI/ref_KAI/g' $file
done

# sed 's/KAI/ref_KAI/g' $file