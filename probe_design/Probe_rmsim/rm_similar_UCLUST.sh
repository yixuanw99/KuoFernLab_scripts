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
input_dir="/home2/yxwu/Probe_rmsim/10gene_0327"
output_dir="/home2/yxwu/Probe_rmsim/10gene_0327_UCLUST_output"
#-------------------

for file in "$input_dir"/*
do
  file_name=${file##*/}
  echo ${file_name}
  awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${file} > "${output_dir}/temp_"${file_name}
  /home2/yxwu/software/usearch -cluster_fast "${output_dir}/temp_"${file_name} -id 0.90 -centroids "${output_dir}/"${file_name} -uc "${output_dir}/"${file_name%.fas}".uc"
  rm "${output_dir}/temp_"${file_name}
done
