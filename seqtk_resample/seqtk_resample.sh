#!/bin/bash

# 定義函數
seqtk_sample() {
    echo "================= subsample preparing =================="
    local trimmed_fq_dir=$1
    local sample=$2
    local project=$3
    local repetition=$4
    local subsample_reads_number=$5

    output_dir="$working_dir/sample${subsample_reads_number}_fqgz"
    mkdir -p $output_dir

    for ((i=1; i<=$repetition; i++))
    do
        echo "ramdon seed $i generating ${subsample_reads_number} reads from ${sample}_*_trimmed.fastq.gz"
        seqtk sample -s$i $trimmed_fq_dir/${sample}_R1_trimmed.fastq.gz $subsample_reads_number | gzip > $output_dir/${sample}_${project}_rep${i}_R1.fq.gz
        seqtk sample -s$i $trimmed_fq_dir/${sample}_R2_trimmed.fastq.gz $subsample_reads_number | gzip > $output_dir/${sample}_${project}_rep${i}_R2.fq.gz
    done
}

run_hybpiper() {
    echo "================= run hybpiper =================="
    local sample=$1
    local project=$2
    local repetition=$3

    for ((i=1; i<=$repetition; i++))
    do
        echo "run hybpiper for $project rep${i}"
        hybpiper assemble -r $working_dir/sample1500000_fqgz/${sample}_${project}_rep${i}_R*.fq.gz -t_dna targets.fas --prefix ${sample}_${project}_rep${i} --run_intronerate
        echo ${sample}_${project}_rep${i} >> $working_dir/namelist.txt
    done
}


# 主程式

# 檢查是否提供了正確的參數
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

working_dir="/home/yixuan/Hybpiper/sample_3m"


# 讀取配置文件
config_file="$1"

# 檢查配置文件是否存在
if [ ! -f "$config_file" ]; then
    echo "Config file not found: $config_file"
    exit 1
fi


# 遍歷配置文件的每一行
while IFS= read -r line || [ -n "$line" ]; do
    # 解析行中的參數
    trimmed_fq_dir=$(echo "$line" | awk '{print $1}')
    sample=$(echo "$line" | awk '{print $2}')
    project=$(echo "$line" | awk '{print $3}')
    repetition=$(echo "$line" | awk '{print $4}')

    seqtk_sample $trimmed_fq_dir $sample $project $repetition 1500000
    run_hybpiper $sample $project $repetition
    rm -r $working_dir/sample1500000_fqgz



done < "$config_file"
