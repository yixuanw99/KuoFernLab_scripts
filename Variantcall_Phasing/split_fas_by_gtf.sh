#!/bin/bash

# 找到序列中第 x 個不是 "-" 的字元在算上 "-" 之後會是在第幾個位置
find_nth_non_dash() {
    local seq="$1"
    local n="$2"
    
    local count=0
    local pos=0
    
    while IFS= read -n1 char; do
        pos=$((pos+1))
        if [ "$char" != "-" ]; then
            count=$((count+1))
            if [ $count -eq $n ]; then
                echo "$pos"
                return
            fi
        fi
    done <<<"$seq"
    
    echo "Error: The $n-th non-dash character not found in the sequence."
}

# 函數：統計SPAdes序列中gap的數量
calculate_gaps() {
    local spades_seq="$1"
    local start="$2"
    local end="$3"
    
    # local gaps_before_start=$(echo "${spades_seq:0:start}" | grep -o "-" | wc -l)
    # local gaps_before_end=$(echo "${spades_seq:0:end}" | grep -o "-" | wc -l)
    
    # local actual_start=$((start + gaps_before_start))
    # local actual_end=$((end + gaps_before_end))

    local actual_start=$(find_nth_non_dash "$spades_seq" "$start")
    if [ $end -eq -1 ]; then
        local actual_end=${#spades_seq}
    else
        local actual_end=$(find_nth_non_dash "$spades_seq" "$end")
    fi

    echo "$actual_start,$actual_end"
}

# 函數：切割FASTA文件並根據GTF文件生成子序列
split_fasta() {
    local fasta_file="$1"
    local gtf_file="$2"
    

    # 讀取並解析SPAdes序列
    local spades_seq=$(grep -A 1 "SPAdes" "$fasta_file" | tail -n 1 | tr -d '\n')

    total_lines=$(wc -l < $gtf_file)
    local prev_end=0
    local interval_count=1

    # 逐行讀取 gtf 文件
    while read -r line; do
        # 篩選出包含 "exon" 的行
        if echo "$line" | grep -q "exon"; then
            # start=$(echo "$line" | awk '{print $4}')
            # end=$(echo "$line" | awk '{print $5}')

            # 提取 start 和 end 位置
            start=$((prev_end+1))
            end=$(echo "$line" | awk '{print $5}')
            output_file="$dir/phased/${sample}_${locus}_${start}to${end}_trimmed.fasta"

            if [ "$interval_count" -eq "$total_lines" ]; then
                end=-1
            fi
            interval_count=$((interval_count + 1))
            prev_end=$end

            
            # 計算SPAdes序列中的實際位置
            gap_info=$(calculate_gaps "$spades_seq" "$start" "$end")
            actual_start=$(echo "$gap_info" | cut -d, -f1)
            actual_end=$(echo "$gap_info" | cut -d, -f2)
            actual_length=$((actual_end - actual_start + 1))

            # 切分FASTA文件
            awk 'NR % 2 == 1 {print $0} NR % 2 == 0 {print substr($0, start, len)}' start=$actual_start len=$actual_length "$fasta_file" | sed '/SPAdes/,+1d' > "$output_file"

        fi
    done < "$gtf_file"


}


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
    project=$(echo "$line" | awk '{print $2}')
    ploidy=$(echo "$line" | awk '{print $3}')

    # 設置路徑
    sample_dir="$working_dir/$sample"

    # 遍歷 sample 目錄下的所有資料夾
    for dir in "$sample_dir"/*; do
        # 檢查是否是資料夾
        if [ -d "$dir" ]; then
            # 
            dir_parts=(${dir//\// })
            locus="${dir_parts[-1]}"

            if [ ! -d "$dir/phased" ]; then
                input_file="/home/yixuan/Hybpiper/LMRPS_capture/$project/$sample/$locus/$sample/sequences/intron/${locus}_supercontig.fasta"
                if [ ! -f $input_file ]; then
                    echo "no supercontig for locus: $locus and sample: $sample"
                    continue
                fi
                mkdir -p "$dir/phased"
                cp $input_file "$dir/phased/${sample}_${locus}_novariants.fasta"
                for ((i=1; i<=$ploidy; i++))
                do
                    sed "s/>.*$/>${sample}_${locus}_haplotype${i}/" "$dir/phased/${sample}_${locus}_novariants.fasta" > "$dir/phased/${sample}_${locus}_haplotype${i}.fasta"
                done
                echo "creating ${sample}_${locus}_novariants.fasta"
                cat "$dir/phased/${sample}_${locus}_haplotype"*".fasta" > "$dir/phased/${sample}_${locus}_novariants.fasta"
                continue
            else
                if [ -s "$dir/phased/"*_novariants.fasta ]; then
                    echo "$locus of $sample in $dir/phased/_novariants.fasta existed"
                    continue
                # 檢查是否存在 gtf，沒有的話就放一個nogtf的fasta
                elif [ ! -s "$dir/whatshap/"*.gtf ]; then
                    echo "no gtf for locus: $locus of $sample in $dir/whatshap/.gtf"
                    seqtk seq "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta" > "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta"
                    sed '/SPAdes/,+1d' "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta" > "$dir/phased/${sample}_${locus}_nogtf.fasta"
                    mv "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta" "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta"
                    continue
                fi
            fi

            # 創建臨時 FASTA 文件
            seqtk seq "$dir/phased/aligned_pooled_"*"_withref.fasta" > "$dir/phased/temp.fasta"

            # main
            split_fasta "$dir/phased/temp.fasta" "$dir/whatshap/"*".gtf"
            mv "$dir/phased/temp.fasta" "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta"

        fi
    done
done < "$config_file"