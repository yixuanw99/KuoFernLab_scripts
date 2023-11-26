#! /bin/bash

awk '{gsub(/:/,"",$1);print $1 > "OG/"$1".fna";}' Orthogroups_tab.txt
awk '{print $1}' Orthogroups_tab.txt

awk 'BEGIN{system("echo 'start'");} {gsub(/:/,"",$1);print $1 > "OG/"$1".fna";} END{system("echo 'end'");}' Orthogroups_tab.txt 
awk 'BEGIN{system("echo 'start'");} {print $1;gsub(/:/,"",$1);for(i=2; i< NF+1; i++){system("grep -A 1 -h -w "$i" /home2/yxwu/test/cds/*")> "OG/txt.txt";};printf("\n");} END{system("echo 'end'");}' Orthogroups_tab.txt 

awk 'BEGIN{system("echo 'start'");} NR==31492{print $1;gsub(/:/,"",$1);for(i=2; i< NF+1; i++){system("grep -A 1 -h -w "$i" /home2/yxwu/test/cds/*");};printf("\n");} END{system("echo 'end'");}' Orthogroups_tab.txt 
line_number=$(awk 'END{print FNR}' Orthogroups_tab.txt)
for ((j=66666; j<66669; j++))
do
	echo $j
	awk 'NR==""{for(i=2; i< NF+1; i++){system("grep -A 1 -h -w "$i" /home2/yxwu/test/cds/*");};}' Orthogroups_tab.txt 
done

#only do this line
awk '{gsub(/:/,"",$1);for(i=2; i< NF+1; i++){system("grep -A 1 -h -w "$i" /home2/yxwu/test/cds/* >> OG/"$1".fna");};}' Orthogroups_tab.txt 
