#! /bin/bash

#################################################################
# Function :use MACSE to align fasta in dir                     #
# Platform :All Linux Based Platform                            #
# Version  :1.0                                                 #
# Date     :2023-03-06                                          #
# Author   :Li-yaung kuo                                        #
# Contact  :yixuan.w99@gmail.com                                #
#################################################################

#-------------------
input_dir="/home2/yxwu/findOG/OGfast/OG_eqone_Lomariopsis_KAI"
output_dir="/home2/yxwu/findOG/OGfast/OG_eqone_Lomariopsis_KAI_MACSE"
#-------------------

for orthogroup in "${input_dir}"/*.fna
do
  file_name=${orthogroup##*/}
  OG=${file_name%%\.fna}
  echo $file_name
  # echo $OG
  java -jar /home/lykuo/MACSE_V2_PIPELINES/UTILS/macse_v2.03.jar -prog alignSequences -seq ${orthogroup} -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2 -out_NT "${output_dir}"/${OG}_MACSE.fas
done
