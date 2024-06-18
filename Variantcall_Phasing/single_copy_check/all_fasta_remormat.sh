#!/bin/bash
# 主程式
# 檢查是否提供了正確的參數
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

working_dir="/home2/yxwu/Phasing"

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
    sample=$(echo "$line" | awk '{print $1}')
    sample_dir="$working_dir/$sample"


    # 遍歷 sample 目錄下的所有資料夾
    for dir in "$sample_dir"/*; do
        # 檢查是否是資料夾
        if [ -d "$dir" ]; then
            # 
            dir_parts=(${dir//\// })
            locus="${dir_parts[-1]}"


            # 檢查是否存在 gtf，沒有的話就放一個nogtf的fasta
            if [ -s "$dir/phased/"*"_withref.fasta" ]; then
                echo "phased oneline fasta for locus: $locus"
                seqtk seq "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta" > "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_oneline.fasta"
            else
                echo "no phased fasta for locus: $locus of $sample in $dir/phased/_withref.fasta"
            fi
        fi
    done
done < "$config_file"


