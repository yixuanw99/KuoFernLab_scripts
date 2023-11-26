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
input_dir="/home2/yxwu/findOG/OGfast/OG_ref"
#-------------------

#if previous output exist, remove it
[[ -f "result_tab.txt" ]] && rm result_tab.txt

#make title of the table
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "orthogroup" "header" "exon" "length" "intron" "length" >> result_tab.txt

#add information into table
for file in "$input_dir"/*
do
  file_name=${file##*/} #remove the address string
  printf "%s\n" ${file_name} >> result_tab.txt #add file_name into table
  sequence_number=$(grep ">" -c $file) #get the sequence number of file
  for ((i=1; i<=${sequence_number}; i++))
  do
    header_index=$((${i}*2-1))
    echo "processing"$(sed -n "${header_index}p" $file)
    printf "\t%s\n" $(sed -n "${header_index}p" $file) >> result_tab.txt #add header into table
    sequence_index=$((${i}*2))
    sed -n "${sequence_index}p;${sequence_index}q" $file | awk 'BEGIN{count=1} \
	{if(match($0,/N+/))\
	{while(match($0,/N+/)) \
	{printf("\t\t%s%s\t%d\t%s%s\t%d\n","exon",count,RSTART-1,"intron",count,RLENGTH); $0=substr($0,RSTART+RLENGTH);count++;}\
	match($0,/[ATCG]$/); printf("\t\t%s%s\t%d\n","exon",count,RSTART)}else{printf("\t\t%s\t%d\n","total_len",length($0)-1);exit;}}'\
	>> result_tab.txt #add exon and intron length into table
  done
done
