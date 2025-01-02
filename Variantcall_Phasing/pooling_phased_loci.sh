#!/bin/bash

# 檢查是否提供了正確的參數
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

working_dir="/home2/yxwu/Phasing_MS"
single_copy_list=("ApPEFP_C" "OG0000764" "OG0000806" "OG0001134" "OG0001287" "OG0001587" "OG0003109" "OG0003385" "OG0003734" "OG0003748" "OG0003780" "OG0004009" "OG0004248" "OG0004360" "OG0004561" "OG0004566" "OG0004614" "OG0004903" "OG0004943" "OG0004969" "OG0005419" "OG0005470" "OG0005474" "OG0005881" "OG0006029" "OG0006334" "OG0006337" "OG0006377" "OG0006406" "OG0006606" "OG0006650" "OG0006781" "OG0006950" "OG0006968" "OG0007034" "OG0007059" "OG0007255" "OG0007337" "OG0007419" "OG0007472" "OG0007491" "OG0007526" "OG0007580" "OG0007584" "OG0007607" "OG0007648" "OG0007814" "OG0007872" "OG0008201" "OG0008464" "OG0008477" "OG0008526" "OG0008551" "OG0008618" "OG0008661" "OG0008757" "OG0008877" "OG0008889" "OG0008981" "OG0008997" "OG0009191" "OG0009405" "OG0009430" "OG0009437" "OG0009708" "OG0009761" "OG0009864" "OG0009975" "OG0010004" "OG0010124" "OG0010148" "OG0010266" "OG0010273" "OG0010275" "OG0010383" "OG0010439" "OG0010506" "OG0010519" "OG0010559" "OG0010562" "OG0010691" "OG0010995" "OG0011134" "OG0011178" "OG0011263" "OG0011365" "OG0011418" "OG0011530" "OG0011599" "OG0011616" "OG0011638" "OG0011683" "OG0011739" "OG0011740" "OG0011920" "OG0011936" "OG0011960" "OG0011994" "OG0012575" "OG0012619" "OG0012642" "OG0012644" "OG0012674" "OG0012719" "OG0012796" "OG0012887" "OG0012889" "OG0012896" "OG0012930" "OG0013120" "OG0013744" "OG0013794" "OG0013810" "OG0013859" "OG0013866" "OG0013892" "OG0014169" "OG0014204" "OG0014205" "OG0014282" "OG0015475")

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
    project=$(echo "$line" | awk '{print $2}')
    ploidy=$(echo "$line" | awk '{print $3}')

    output_dir="${working_dir}/pooling_allsample"
    mkdir -p "$output_dir"

    for singlecopygene in "${single_copy_list[@]}"; do
        singlecopy_dir="$output_dir/$singlecopygene"
        mkdir -p "$singlecopy_dir"
        
        swaped_file="${working_dir}/${sample}/${singlecopygene}/swap/pooled.fasta"
        supercontig_file="${working_dir}/${sample}/${singlecopygene}/bwa/${sample}_${singlecopygene}_supercontig.fasta"
        haploidy_file="${working_dir}/${sample}/${singlecopygene}/phased/${sample}_${singlecopygene}_haploidy.fasta"

        if [ -f "$swaped_file" ]; then
            echo "coping swaped_file for $sample and $singlecopygene"
            cat "$swaped_file" | sed "s/seq/${sample}_${singlecopygene}_phased/" >> "$singlecopy_dir/allsample_pooled.fasta"
        elif [ -f "$supercontig_file" ]; then
            echo "coping supercontig_file for $sample and $singlecopygene"
            cat "$supercontig_file" >> "$singlecopy_dir/allsample_pooled.fasta"
        elif [ -f "$haploidy_file" ]; then
            echo "coping haploidy_file for $sample and $singlecopygene"
            cat "$haploidy_file" >> "$singlecopy_dir/allsample_pooled.fasta"
        else
            echo "Error: $swaped_file or $supercontig_file or $haploidy_file not found for $sample and $singlecopygene"
        fi
    done

done < "$config_file"

for singlecopygene in "${single_copy_list[@]}"; do
    concat_dir="$working_dir/concat"
    mkdir -p "$concat_dir"

    echo "align ${singlecopygene} pooled fasta"
    mafft "$working_dir/pooling_allsample/$singlecopygene/allsample_pooled.fasta" > "$concat_dir/${singlecopygene}.fasta"
done