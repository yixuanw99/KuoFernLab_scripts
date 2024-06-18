#!/bin/bash

# 遍歷配置文件的每一行
while IFS= read -r line || [ -n "$line" ]; do
    # 解析行中的參數
    sample=$(echo "$line" | awk '{print $1}')
    project=$(echo "$line" | awk '{print $2}')
    ploidy=$(echo "$line" | awk '{print $3}')

    # 設置raw_fq目錄的路徑
    raw_fq_dir="/home2/yxwu/sequence_raw/Miseq_raw/${project}_FASTQ/trimmed"

    # 設置sample目錄的路徑
    sample_dir="/home/yixuan/Hybpiper/LMRPS_capture/$project/$sample"

    # 遍歷sample目錄下的所有資料夾
    for dir in "$sample_dir"/*; do
        # 檢查是否是資料夾
        if [ -d "$dir" ]; then
            # 解析目錄名稱以獲取 locus
            dir_parts=(${dir//\// })
            locus="${dir_parts[-1]}"
            
            # 構建樣本與基因座位的路徑
            input_file="$dir/$sample/sequences/intron/${locus}_supercontig.fasta"
            RG_ID="$project"


            if [ ! -f $input_file ]; then
                echo "no supercontig for locus: $locus and sample: $sample"
                continue
            else
                if [ $ploidy -eq 1 ];then
                    mkdir -p "$working_dir/$sample/$locus/phased"
                    cp $input_file "$working_dir/$sample/$locus/phased/${sample}_${locus}_haploidy.fasta"
                    echo "copied $input_file of haploid sample(${sample})"
                    continue
                fi
                # 創建目標目錄(如果不存在)
                mkdir -p "$working_dir/$sample/$locus"
                
                # 複製並索引參考基因組序列
                copy_and_index_reference "$input_file" "$sample" "$locus"
                
                # 將讀段與參考基因組比對mapping
                bwa_mapping_reads "$RG_ID" "$sample" "$locus"
                
                # 標記重複讀段
                gatk_mark_duplicates "$sample" "$locus"
                
                # 呼叫變異位點
                gatk_call_variants "$ploidy" "$sample" "$locus"
                
                # 過濾變異位點
                gatk_filter_variants "$sample" "$locus"
                
                # 對變異位點進行分型
                whatshap_phase_variants "$ploidy" "$sample" "$locus"
                
                # whatshap 可視化
                whatshap_visualization "$sample" "$locus"
                
                # 生成haplotype並建立一致性序列
                bcftools_concensus "$input_file" "$sample" "$locus" "$ploidy"
            fi
        fi
    done
done < "$config_file"
