#! /bin/bash

##################################################################
# Function :MACSE_removeKAIcds                                   #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-03-09                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_dir="/home2/yxwu/findOG/OGfast/OG_eqonewithintron_all_MACSE/oneline"
output_dir="/home2/yxwu/findOG/OGfast/OG_eqonewithintron_all_MACSE_removeKAIcds"
#-------------------

for file in "${input_dir}"/*
do
  file_name=${file##*/}
  echo ${file_name}
  sed '/ref/,+1 d' ${file} > "${output_dir}"/${file_name}
done
