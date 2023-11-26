#! /bin/bash

OG_list=($(awk '{if (NR!=1) {print $1}}' Orthogroups.tsv))
sample_list=($(awk 'NR==1{for(i=2;i<=NF;i++){printf("%s\n",$i)}}' Orthogroups.tsv))
echo ${#sample_list[@]}
echo ${#OG_list[@]}

for ((j=2; j<5; j++))
do
	for ((k=2; k<$((${#OG_list[@]}+2)); k++))
	do
		sed -n "${k}p" Orthogroups.tsv | awk -v j_index=$(($j+2)) 'BEGIN{FS="\t"} {print $j_index}' > pattern_file_Lom.txt
		filesize=$(ls -l pattern_file_Lom.txt | awk '{print $5}')
		echo $filesize
		if [[ ${filesize} -gt 2 ]];then
		  sed -i "s/, /\\n/g" pattern_file_Lom.txt
		  cat pattern_file_Lom.txt
		  grep -A 1 -w -f pattern_file_Lom.txt --no-group-separator /home2/yxwu/findOG/cds_with_intron/${sample_list[j]}.fas >> OG_outgroup/${OG_list[$(($k-2))]}.fna
		else
		  echo "-----------"
		fi
		rm pattern_file_Lom.txt
	done
done
