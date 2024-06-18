#!/bin/bash

# 函數：統計SPAdes序列中gap的數量
calculate_gaps() {
    local spades_seq="$1"
    local start="$2"
    local end="$3"
    
    local gaps_before_start=$(echo "${spades_seq:0:start}" | grep -o "-" | wc -l)
    local gaps_before_end=$(echo "${spades_seq:0:end}" | grep -o "-" | wc -l)
    
    local actual_start=$((start + gaps_before_start))
    local actual_end=$((end + gaps_before_end))

    echo "$actual_start,$actual_end"
}

# 函數：切割FASTA文件並根據GTF文件生成子序列
split_fasta() {
    local fasta_file="$1"
    local gtf_file="$2"

    # 讀取並解析SPAdes序列
    local spades_seq=$(grep -A 1 "SPAdes" "$fasta_file" | tail -n 1 | tr -d '\n')

    # 逐行讀取 gtf 文件
    while read -r line; do
        # 篩選出包含 "exon" 的行
        if echo "$line" | grep -q "exon"; then
            start=$(echo "$line" | awk '{print $4}')
            end=$(echo "$line" | awk '{print $5}')
            # 檢查start和end是否為數字
            if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]]; then
                echo "Invalid start or end position: $line"
                continue
            fi
            start=$((start - 1))  # 將起始位置轉換為0-based索引
            end=$((end - 1))      # 將結束位置轉換為0-based索引
            
            # 計算SPAdes序列中的實際位置
            gap_info=$(calculate_gaps "$spades_seq" "$start" "$end")
            actual_start=$(echo "$gap_info" | cut -d, -f1)
            actual_end=$(echo "$gap_info" | cut -d, -f2)

            # 切分FASTA文件
            output_file="$dir/phased/${sample}_${locus}_${start}to${end}_alignment.fasta"
            awk 'NR % 2 == 1 {print $0} NR % 2 == 0 {print substr($0, start, len)}' start=$actual_start len=$actual_end "$fasta_file" | sed '/SPAdes/,+1d' > "$output_file"

            # # 提取實際序列
            # actual_sequence=$(echo "$spades_seq" | cut -c$((actual_start+1))-$((actual_end+1)))

            # # 生成輸出FASTA文件
            # output_file="${seqname}_${start}_${end}.fasta"
            # echo ">${seqname}_${start}_${end}" > "$output_file"
            # echo "$actual_sequence" >> "$output_file"
        fi
    done < "$gtf_file"
}



# 函數: read_largest_block
# 功能: 從指定的 gtf 文件中讀取最大的phased block區間長度
# 用法: read_largest_block <phased.gtf_file>
read_largest_block() {
    # 檢查是否提供了正確的參數
    if [ $# -ne 1 ]; then
        echo "Usage: $0 <phased.gtf_file>"
        exit 1
    fi

    # 讀取 gtf 文件
    phased_file="$1"

    # 檢查配置文件是否存在
    if [ ! -f "$phased_file" ]; then
        echo "gtf file not found: $phased_file"
        exit 1
    fi

    max_start=0
    max_len=0

    # 逐行讀取 gtf 文件
    while read -r line; do
        # 篩選出包含 "exon" 的行
        if echo "$line" | grep -q "exon"; then
            start=$(echo "$line" | awk '{print $4}')
            end=$(echo "$line" | awk '{print $5}')
            # 檢查start和end是否為數字
            if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]]; then
                echo "Invalid start or end position: $line"
                continue
            fi
            length=$((end - start + 1))

            # 更新最大長度
            if [ "$length" -gt "$max_len" ]; then
                max_len=$length
                max_start=$start
            fi
        fi
    done < $phased_file
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
            # 執行 process_sample_locus 命令
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
                # 檢查是否存在 gtf
                if [ ! -s "$dir/whatshap/"*.gtf ]; then
                    echo "no gtf for locus: $locus of $sample in $dir/whatshap/.gtf"
                    seqtk seq "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta" > "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta"
                    sed '/SPAdes/,+1d' "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta" > "$dir/phased/${sample}_${locus}_nogtf.fasta"
                    mv "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref_temp.fasta" "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta"
                    continue
                fi
            fi

            # 創建臨時 FASTA 文件
            seqtk seq "$dir/phased/aligned_pooled_"*"_withref.fasta" > "$dir/phased/temp.fasta"

            # 重置 max_start 和 max_len 變量
            unset max_start max_len

            # 讀取最大的phased block區間長度
            read_largest_block "$dir/whatshap/"*".gtf"

            # 根據最大區間長度截取序列
            awk 'NR % 2 == 1 {print $0} NR % 2 == 0 {print substr($0, start, len)}' start=$max_start len=$max_len "$dir/phased/temp.fasta" | sed '/SPAdes/,+1d' > "$dir/phased/${sample}_${locus}_trimmed_alignment.fasta"
            mv "$dir/phased/temp.fasta" "$dir/phased/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta"
        fi
    done
done < "$config_file"