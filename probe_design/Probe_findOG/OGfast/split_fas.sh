#! /bin/bash

##################################################################
# Function :based on exon# to assign each sequence into fasta    #
# Platform :All Linux Based Platform                             #
# Version  :1.0                                                  #
# Date     :2023-02-24                                           #
# Author   :Yihsuan Wu                                           #
# Contact  :yixuan.w99@gmail.com                                 #
##################################################################

#-------------------
input_dir="/home2/yxwu/findOG/OGfast/OG_gtzero_KAI_exon_split"
output_dir="/home2/yxwu/findOG/OGfast/OG_gtzero_KAI_exon_split/fas_split"
#-------------------

#找到有exon的fasta檔案
#在裡面找exon1~n，分派到不同fasta裡，檔名叫做OGxxxxxxx_exon#.fasta

grep -l -r "_exon" ${input_dir} > with_exon_OGlist.txt

while read file
do
  file_name=${file##*/}
  number=$(awk 'END{print NR}' $file)
  echo ${file##*/}" line: "$number
  for ((i=1; i<=$number/2; i++))
  do
    i2=$(($i*2))
	out_file_name="${file_name%.*}"_$(sed -n "$(($i2-1))p" $file | sed 's/>//')
    sed -n "$(($i2-1)),${i2}p" $file > "${output_dir}"/"${out_file_name}".fas
  done
done < with_exon_OGlist.txt

echo "done"
rm with_exon_OGlist.txt
