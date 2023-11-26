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
input_dir="/home/yixuan/OrthoFinder/UCLUST_replaceby_longest/input"
output_dir="/home/yixuan/OrthoFinder/UCLUST_replaceby_longest/output"
#-------------------

get_Max_index_in_Array(){
arr=("$@")
index=0
max=${arr[0]}
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

for file in "$input_dir"/*
do
  file_name=${file##*/}
  echo "processing: "${file_name}
  sed '/KAI/,+1 d' ${file} > ${output_dir}/${file_name}
  cp $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.cen ${output_dir}/al_L100_OG0006743_combined.fasta.trim.keep.cen
  sed 's/\*$/cluster/g' $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc > $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp
  SHC_column=($(awk '/^H\t/{print $1}' $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp))
  similarSeq_column=($(awk '/^H\t/{print $9}' $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp))
  centroidSeq_column=($(awk '/^C\t/{print $9}' $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp))

  # echo ${#centroidSeq_column[@]}


  for centroid in ${centroidSeq_column[@]}
  do
    unset length_list
    similarSeq=""
    echo "============================="${centroid}
    similarSeq=($(awk "/^H/&&/${centroid}/{print \$9}" $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp))
    for ((i=0; i<${#similarSeq[@]}; i++))
    do
      length_list+=($(printf "%d" "${similarSeq[$i]##*size.}"))
    done
    length_list+=($(echo "${centroid##*size.}"))
    echo "length list: "${length_list[@]}
    similarSeq_index=$(get_Max_index_in_Array "${length_list[@]}")
    if [[ -z "${similarSeq_index}" ]]; then
    echo "not modified"
    else
      echo "replaced"
      sed -i "/${centroid}/,+1 d" ${output_dir}/al_L100_OG0006743_combined.fasta.trim.keep.cen
      grep -A1 "${similarSeq[$similarSeq_index]}" $input_dir/al_L100_OG0006743_combined.fasta.trim.keep >> ${output_dir}/al_L100_OG0006743_combined.fasta.trim.keep.cen
    fi
  done

  rm $input_dir/al_L100_OG0006743_combined.fasta.trim.keep.uc_temp
done