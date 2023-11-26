#! /bin/bash

####################################################################
# Function :count probe number                                     #
# Platform :All Linux Based Platform                               #
# Version  :1.0                                                    #
# Date     :2023-03-24                                             #
# Author   :Yihsuan Wu                                             #
# Contact  :yixuan.w99@gmail.com                                   #
####################################################################

#-------------------
input_dir="/home2/yxwu/Probe_rmsim/replaced_1_perfect_all_rmKAI_rmSimilar"
output_dir="/home2/yxwu/Probe_rmsim"
#-------------------

#if previous output exist, remove it
[[ -f $output_dir/"result_tab2.txt" ]] && rm $output_dir/result_tab2.txt

#make title of the table
printf "%s\t%s\t%s\n" "file" "seq number" "total length" >> $output_dir/result_tab2.txt

for file in "$input_dir"/*.fas
do
  echo "processing: $file"
  file_name=${file##*/}
  seq_number=$(grep -c ">" $file)
  total_length_temp=$(grep -v ">" $file | sed 's/-//g' | wc -c)
  total_length=$(($total_length_temp-1))
  printf "%s\t%s\t%s\n" "$file_name" "$seq_number" "$total_length" >> $output_dir/result_tab2.txt
done
