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
    local seq="$1"
    local start="$2"
    local end="$3"
    

    local actual_start=$(find_nth_non_dash "$seq" "$start")
    if [ $end -eq -1 ]; then
        local actual_end=${#seq}
    else
        local actual_end=$(find_nth_non_dash "$seq" "$end")
    fi

    echo "$actual_start,$actual_end"
}


swap_block_by_gtf() {
    GTF_FILE=$1
    FASTA_FILE=$2
    REF_FILE=$3
    sample=$4

    # 提取 GTF 中的區間
    INTERVALS=$(awk '{print $4"-"$5}' $GTF_FILE)
    if [[ -z "$INTERVALS" ]]; then
        # INTERVALS=$(awk 'NR==2 {print "1-" length($0)}' $FASTA_FILE)
        INTERVALS=$(grep -A 1 $sample $REF_FILE | awk '/^>/ {getline; gsub(/-/, "", $0); print "1-" length($0)}')

    fi

    # 建立臨時目錄
    TMP_DIR=$(mktemp -d)

    # 拆分 FASTA 文件到單獨的序列文件
    csplit -s -z -f $TMP_DIR/seq_ $FASTA_FILE '/^>/' '{*}'

    # 處理 FASTA 文件，移除標題行並合併成單行
    for SEQ_FILE in $TMP_DIR/seq_*; do
        mv $SEQ_FILE $SEQ_FILE.tmp
        tail -n +2 $SEQ_FILE.tmp | tr -d '\n' > $SEQ_FILE
        rm $SEQ_FILE.tmp
    done

    # 計算每個區間的相似性並排序
    for INTERVAL in $INTERVALS; do
        START=$(echo $INTERVAL | cut -d'-' -f1)
        END=$(echo $INTERVAL | cut -d'-' -f2)


        supercontig_seq=$(grep -A 1 $sample $REF_FILE | awk '/^>/ {getline; print}' | tr -d '\n')
        echo $supercontig_seq $START $END

        # 計算SPAdes序列中的實際位置
        gap_info=$(calculate_gaps "$supercontig_seq" "$START" "$END")
        START=$(echo "$gap_info" | cut -d, -f1)
        END=$(echo "$gap_info" | cut -d, -f2)
        echo $START $END

        > $TMP_DIR/similarities_$INTERVAL.txt  # 清空文件
        
        for SEQ_FILE in $TMP_DIR/seq_*; do
            SEQ=$(cat $SEQ_FILE)
            SEGMENT=${SEQ:$START-1:$END-$START+1}
            REF_SEQUENCE=$(grep -A 1 "YXD136" $REF_FILE | awk '/^>/ {getline; print}' | tr -d '\n')
            # echo "REF_SEQUENCE:$REF_SEQUENCE"
            REF_SEGMENT=${REF_SEQUENCE:$START-1:$END-$START+1}
            
            # 計算相似性（簡單計算共有多少匹配的字符）
            SIMILARITY=0
            for (( i=0; i<${#SEGMENT}; i++ )); do
                if [ "${SEGMENT:$i:1}" == "${REF_SEGMENT:$i:1}" ]; then
                    SIMILARITY=$((SIMILARITY+1))
                fi
            done
            
            # 計算相似性比例
            SIMILARITY_RATIO=$(echo "$SIMILARITY / ${#SEGMENT}" | bc -l)
            
            echo "$SEQ_FILE $SIMILARITY_RATIO" >> $TMP_DIR/similarities_$INTERVAL.txt
        done
        
        # 按相似性排序
        sort -k2,2nr $TMP_DIR/similarities_$INTERVAL.txt -o $TMP_DIR/similarities_sorted_$INTERVAL.txt
        
        # 交換序列
        SORTED_SEQ_FILES=($(awk '{print $1}' $TMP_DIR/similarities_sorted_$INTERVAL.txt))
        INDEX=0
        for SEQ_FILE in $TMP_DIR/seq_*; do
            SEQ=$(cat $SEQ_FILE)
            NEW_SEGMENT=$(cat ${SORTED_SEQ_FILES[$INDEX]} | cut -c $START-$END)
            echo "START-END:$START-$END,NEW_SEGMENT:$NEW_SEGMENT,${SORTED_SEQ_FILES[$INDEX]}"
            INDEX=$((INDEX + 1))
            # 替換區間內的序列
            NEW_SEQ=${SEQ:0:$START-1}$NEW_SEGMENT${SEQ:$END}
            # echo "$TMP_DIR/${INDEX}_$INTERVAL"
            echo $NEW_SEQ > $TMP_DIR/${INDEX}_$INTERVAL
        done

        # 用新的序列替換舊的序列
        INDEX=1
        for SEQ_FILE in $TMP_DIR/seq_*; do
            # echo $SEQ_FILE
            mv $TMP_DIR/${INDEX}_$INTERVAL $SEQ_FILE
            INDEX=$((INDEX + 1))
        done
    done


    # 用新的序列替換舊的序列
    for SEQ_FILE in $TMP_DIR/seq_*; do
        # echo "${SEQ_FILE##*/}"
        SEQ=$(cat $SEQ_FILE)
        printf ">${SEQ_FILE##*/}\n${SEQ}\n" > "$output_dir/${SEQ_FILE##*/}"
    done

    cat "$output_dir/seq_"* > "$output_dir/pooled.fasta"
    # 清理臨時目錄
    rm -r $TMP_DIR

}


prepare_swap_fasta() {
    echo "============= 準備下一步要swap的兩個fasta ============="
    local sample="$1"
    local locus="$2"

    # 創建目標目錄(如果不存在)
    output_dir="$working_dir/$sample/$locus/swap"
    mkdir -p $output_dir
    
    mafft --auto --addlong "$working_dir/YXD136/$locus/phased/"*"_haploidy.fasta" --keeplength "$working_dir/$sample/$locus/phased/"*"_withref_oneline.fasta" > "$output_dir/${sample}_${locus}_addlong.fasta"
    seqtk seq "$output_dir/${sample}_${locus}_addlong.fasta" > "$output_dir/${sample}_${locus}_addlong_oneline.fasta"
    rm "$output_dir/${sample}_${locus}_addlong.fasta"
    sed '/SPAdes/,+1d' "$output_dir/${sample}_${locus}_addlong_oneline.fasta" > "$output_dir/${sample}_${locus}_onlyhaplotypes.fasta"
    sed "/haplotype/,+1d" "$output_dir/${sample}_${locus}_addlong_oneline.fasta" > "$output_dir/${sample}_${locus}_reference.fasta"
            
}

is_single_copy() {
    echo "================= 檢查是否是single_copy =================="
    local single_copy_list=("ApPEFP_C" "OG0000764" "OG0000806" "OG0001134" "OG0001287" "OG0001587" "OG0003109" "OG0003385" "OG0003734" "OG0003748" "OG0003780" "OG0004009" "OG0004248" "OG0004360" "OG0004561" "OG0004566" "OG0004614" "OG0004903" "OG0004943" "OG0004969" "OG0005419" "OG0005470" "OG0005474" "OG0005881" "OG0006029" "OG0006334" "OG0006337" "OG0006377" "OG0006406" "OG0006606" "OG0006650" "OG0006781" "OG0006950" "OG0006968" "OG0007034" "OG0007059" "OG0007255" "OG0007337" "OG0007419" "OG0007472" "OG0007491" "OG0007526" "OG0007580" "OG0007584" "OG0007607" "OG0007648" "OG0007814" "OG0007872" "OG0008201" "OG0008464" "OG0008477" "OG0008526" "OG0008551" "OG0008618" "OG0008661" "OG0008757" "OG0008877" "OG0008889" "OG0008981" "OG0008997" "OG0009191" "OG0009405" "OG0009430" "OG0009437" "OG0009708" "OG0009761" "OG0009864" "OG0009975" "OG0010004" "OG0010124" "OG0010148" "OG0010266" "OG0010273" "OG0010275" "OG0010383" "OG0010439" "OG0010506" "OG0010519" "OG0010559" "OG0010562" "OG0010691" "OG0010995" "OG0011134" "OG0011178" "OG0011263" "OG0011365" "OG0011418" "OG0011530" "OG0011599" "OG0011616" "OG0011638" "OG0011683" "OG0011739" "OG0011740" "OG0011920" "OG0011936" "OG0011960" "OG0011994" "OG0012575" "OG0012619" "OG0012642" "OG0012644" "OG0012674" "OG0012719" "OG0012796" "OG0012887" "OG0012889" "OG0012896" "OG0012930" "OG0013120" "OG0013744" "OG0013794" "OG0013810" "OG0013859" "OG0013866" "OG0013892" "OG0014169" "OG0014204" "OG0014205" "OG0014282" "OG0015475")
    # local single_copy_list=("OG0009708" "OG0009761" "OG0009864" "OG0009975")
    local item=$1
    local found=1

    # 遍歷 single_copy_list 檢查 item 是否在 single_copy_list 中
    for list_item in "${single_copy_list[@]}"; do
        if [ "$item" == "$list_item" ]; then
            found=0
            break
        fi
    done

    return $found
}




# 主程式

# 檢查是否提供了正確的參數
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

working_dir="/home2/yxwu/Phasing_MS"


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

        # 設置sample目錄的路徑
    sample_dir="$working_dir/$sample"

    # 遍歷sample目錄下的所有資料夾
    for dir in "$sample_dir"/*; do
        # 檢查是否是資料夾
        if [ -d "$dir" ]; then
            # 解析目錄名稱以獲取 locus
            dir_parts=(${dir//\// })
            locus="${dir_parts[-1]}"

            if is_single_copy "$locus"; then
                echo "$locus is in single copy list, processing..."
            else
                echo "$locus is not in single copy list, skipping..."
                continue
            fi
            
            if [ ! -f "$working_dir/$sample/$locus/phased/"*"_withref_oneline.fasta" ]; then
                echo "no phased fasta: $locus and sample: $sample"
                continue
                
            fi

            prepare_swap_fasta $sample $locus
            swap_block_by_gtf $working_dir/$sample/$locus/whatshap/*.gtf "$output_dir/${sample}_${locus}_onlyhaplotypes.fasta" "$output_dir/${sample}_${locus}_reference.fasta" $sample
        
        fi
    done
done < "$config_file"
