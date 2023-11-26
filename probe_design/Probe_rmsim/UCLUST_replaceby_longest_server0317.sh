#! /bin/bash

####################################################################
# Function :UCLUST_replaceby_longest                               #
# Platform :All Linux Based Platform                               #
# Version  :1.0                                                    #
# Date     :2023-03-17                                             #
# Author   :Yihsuan Wu                                             #
# Contact  :yixuan.w99@gmail.com                                   #
####################################################################

#-------------------
input_fasta_dir="/home2/yxwu/Probe_codon12_compare/output_eqone_MACSE/1_perfect_all_rmKAI"
input_uc_dir="/home2/yxwu/Probe_rmsim/1_perfect_all_rmKAI_rmSimilar"
input_centroids_dir="/home2/yxwu/Probe_rmsim/oneline_1_perfect_all_rmKAI_rmSimilar"
output_dir="/home2/yxwu/Probe_rmsim/replaced_1_perfect_all_rmKAI_rmSimilar"
#-------------------

get_Max_index_in_Array(){
arr=("$@")
index=0
max=${arr[-1]}
max_pos=""
for ((i=0; i<${#arr[@]}; i++))
do
  if [[ "${arr[i]}" -gt "${max}" ]]; then
    max="${arr[i]}"
    max_pos=$index
  fi
  index=$(($index+1))
done
echo "$max_pos"
}

for file in "$input_uc_dir"/*.uc
do
  file_name=${file##*/}
  echo "processing: "${file_name%\.uc}
  # cp $input_centroids_dir/${file_name%\.uc}.fas ${output_dir}/${file_name%\.uc}.fas
  sed 's/\*$/cluster/g' $file > $file"_temp"
  SHC_column=($(awk '/^H\t/{print $1}' $file"_temp"))
  similarSeq_column=($(awk '/^H\t/{print $9}' $file"_temp"))
  centroidSeq_column=($(awk '/^C\t/{print $9}' $file"_temp"))
  
  for centroid in ${centroidSeq_column[@]}
  do
    grep -A1 "${centroid}" $input_fasta_dir/${file_name%\.uc}.fas >> $output_dir/${file_name%\.uc}.fas
  done
  
  for centroid in ${centroidSeq_column[@]}
  do
    unset length_list
    similarSeq=""
    echo "==============centroid==============="${centroid}
    similarSeq=($(awk "/^H/&&/${centroid}/{print \$9}" $file"_temp"))
    for ((i=0; i<${#similarSeq[@]}; i++))
    do
      length_list+=($(printf "%d" "${similarSeq[$i]##*size}"))
    done
    length_list+=($(echo "${centroid##*size}"))
    echo "length of Seq in cluster: "${length_list[@]}
    similarSeq_index=$(get_Max_index_in_Array "${length_list[@]}")
    if [[ -z "${similarSeq_index}" ]]; then
    echo "not modified"
    else
      echo "replaced by: ${similarSeq[$similarSeq_index]}"
      sed -i "/${centroid}/,+1 d" ${output_dir}/${file_name%\.uc}.fas
      grep -A1 "${similarSeq[$similarSeq_index]}" $input_fasta_dir/${file_name%\.uc}.fas >> ${output_dir}/${file_name%\.uc}.fas
    fi
  done

  rm $file"_temp"
done
