#!/bin/bash

# 定義函數
check_and_create_dir() {
    echo "================= 檢查並創建目錄 =================="
    dir_name="$1"

    if [ -d "$dir_name" ]; then
        legacy_dir="${dir_name}_legacy"
        if [ -d "$legacy_dir" ]; then
            # 如果 legacy 目錄存在,可以選擇覆蓋或退出
            echo "Error: Directory '$legacy_dir' already exists."
            echo "Please remove or rename it before continuing."
            return 1
        fi
        mv "$dir_name" "$legacy_dir"
        echo "Directory '$dir_name' already exists, renamed to '$legacy_dir'"
    fi
    mkdir -p "$dir_name"
    echo "Directory '$dir_name' created"
    echo "======================================================="

}

process_sample_locus() {
    echo "================= 處理樣本與基因座位 =================="
    # ref from hybpiper intronerate (supercontig.fasta)
    # example: process_sample_locus /home/yixuan/Hybpiper/LMRPS_capture/MS23395/YXD136/transducin/YXD136/sequences/intron/transducin_supercontig.fasta MS23395 1
    # 從參數獲取輸入文件路徑
    input_file="$1"
    RG_ID="$2"
    ploidy="$3"

    # 從輸入路徑中解析出 sample 和 locus
    dir_parts=(${input_file//\// })
    sample="${dir_parts[-6]}"
    locus="${dir_parts[-5]}"

    check_and_create_dir "$working_dir"
    mkdir "$working_dir/$sample"
    
    copy_and_index_reference "$input_file" "$sample" "$locus"
    bwa_mapping_reads "$RG_ID" "$sample" "$locus"
    gatk_mark_duplicates "$sample" "$locus"
    gatk_call_variants "$ploidy" "$sample" "$locus"
    gatk_filter_variants "$sample" "$locus"
    whatshap_phase_variants "$sample" "$locus"
    whatshap_visualization "$sample" "$locus"
    echo "======================================================="
}

copy_and_index_reference() {
    echo "============= 複製並索引參考基因組序列 ============="
    local input_file="$1"
    local sample="$2"
    local locus="$3"

    # 構建輸出目錄路徑
    output_dir="$working_dir/$sample/$locus/bwa"

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_dir"

    # 將 ref.fas 複製到目標目錄
    cp "$input_file" "$output_dir"/"${sample}_${locus}_supercontig.fasta"
    echo "======================================================="
}

bwa_mapping_reads() {
    echo "============= 將讀段與參考基因組比對mapping ============="
    local RG_ID="$1"
    local sample="$2"
    local locus="$3"
    # 將讀段與參考基因組(supercontig.fasta)進行比對mapping，只存mapping到的reads(samtools -F 4)
    # 從ref+trimmed_fastq.gz>_mapped.sam

    # 構建輸出目錄路徑
    output_bwa_dir="$working_dir/$sample/$locus/bwa"
    # 進入目標目錄並索引ref
    cd "$output_bwa_dir"
    bwa index -p "${sample}_${locus}" "${sample}_${locus}_supercontig.fasta"

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_bwa_dir"
    mkdir -p "$output_bwa_dir/output"

    #bwa | samtools
    bwa mem -M -t 20 -R "@RG\tID:${RG_ID}\tSM:${sample}_${locus}" \
    "$output_bwa_dir/${sample}_${locus}" \
    $raw_fq_dir/${sample}*"R1_trimmed.fastq.gz" $raw_fq_dir/${sample}*"R2_trimmed.fastq.gz" \
    2> "$output_bwa_dir/bwa.err" \
    | samtools view -h -F 4 -o "$output_bwa_dir/output/${sample}_${locus}_mapped.sam"
    echo "======================================================="
}

gatk_mark_duplicates() {
    echo "============= 標記重複讀段 ============="
    cd "$working_dir/$sample/$locus"
    # 構建輸出目錄路徑
    output_gatk_dir="$working_dir/$sample/$locus/gatk"
    # 標記重複讀段，從_mapped.sam>_marked_duplicates.bam
    local sample="$1"
    local locus="$2"
    
    mkdir -p "$output_gatk_dir"

    gatk MarkDuplicatesSpark \
    -I "$working_dir/$sample/$locus/bwa/output/${sample}_${locus}_mapped.sam" \
    -O "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam" \
    -M "$output_gatk_dir/${sample}_${locus}_marked_dup_metrics.txt"
    echo "======================================================="
}

gatk_call_variants() {
    echo "============= 呼叫變異位點 ============="
    # 呼叫變異位點，從_marked_duplicates.bam>.vcf
    local ploidy="$1"
    local sample="$2"
    local locus="$3"
    cd "$working_dir/$sample/$locus"
    # 構建輸出目錄路徑
    output_gatk_dir="$working_dir/$sample/$locus/gatk"
    gatk CreateSequenceDictionary -R "$working_dir/$sample/$locus/bwa/${sample}_${locus}_supercontig.fasta"
    
    samtools faidx "$working_dir/$sample/$locus/bwa/${sample}_${locus}_supercontig.fasta"
    
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -ploidy $ploidy \
    -R "$working_dir/$sample/$locus/bwa/${sample}_${locus}_supercontig.fasta" \
    -I "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam" \
    -O "$output_gatk_dir/${sample}_${locus}.vcf"
    echo "======================================================="

}

gatk_filter_variants() {
    echo "============= 過濾變異位點 ============="
    # 從.vcf>_filtered.vcf
    local sample="$1"
    local locus="$2"
    cd "$working_dir/$sample/$locus"
    # 構建輸出目錄路徑
    output_gatk_dir="$working_dir/$sample/$locus/gatk"
    # 過濾變異位點
    gatk VariantFiltration \
    --variant "$output_gatk_dir/${sample}_${locus}.vcf" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "AF < 0.05 || AF > 0.95" --filter-name "AF5-95" \
    --output "$output_gatk_dir/${sample}_${locus}_filtered.vcf"
    echo "======================================================="
}

whatshap_phase_variants() {
    echo "============= 對變異位點進行分型 ============="
    # 從_filtered.vcf+_marked_duplicates.bam>_phased.bcf
    local sample="$1"
    local locus="$2"

    cd "$working_dir/$sample/$locus"
    # 構建輸出目錄路徑
    output_gatk_dir="$working_dir/$sample/$locus/gatk"

    cp -r "$working_dir/$sample/$locus/bwa/." "$output_gatk_dir"
    #Index the variants file (skip if a .csi or .tbi file already exists):
    bcftools view "$output_gatk_dir/${sample}_${locus}_filtered.vcf" -Oz -o "$output_gatk_dir/${sample}_${locus}_filtered.vcf.gz"
    bcftools index "$output_gatk_dir/${sample}_${locus}_filtered.vcf.gz"
    
    #Index the alignment file (skip if a .csi or .bai file already exists):
    samtools index "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam"
    
    cd "$output_gatk_dir"

    # 構建輸出目錄路徑
    output_whatshap_dir="$working_dir/$sample/$locus/whatshap"
    # 創建目標目錄(如果不存在)
    mkdir -p "$output_whatshap_dir"

    #Run read-based phasing.
    whatshap phase \
    --reference="$output_gatk_dir/${sample}_${locus}_supercontig.fasta" \
    -o "$output_whatshap_dir/${sample}_${locus}_phased.bcf" \
    "$output_gatk_dir/${sample}_${locus}_filtered.vcf" "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam"
    echo "======================================================="
}

whatshap_visualization() {
    echo "============= whatshap 可視化 ============="
    local sample="$1"
    local locus="$2"


    # 構建輸出目錄路徑
    output_whatshap_dir="$working_dir/$sample/$locus/whatshap"

    cd "$output_whatshap_dir"


    # 用IGV視覺化前置index
    bcftools convert -O v -o "$output_whatshap_dir/${sample}_${locus}_phased.vcf" "$output_whatshap_dir/${sample}_${locus}_phased.bcf"
    whatshap stats --gtf="$output_whatshap_dir/${sample}_${locus}_phased.gtf" "$output_whatshap_dir/${sample}_${locus}_phased.bcf"
    bcftools index "$output_whatshap_dir/${sample}_${locus}_phased.bcf"

    # 用IGV視覺化
    whatshap haplotag -o "$output_whatshap_dir/${sample}_${locus}_haplotagged.bam" \
    --output-haplotag-list="$output_whatshap_dir/${sample}_${locus}_haplotag_table.tsv" \
    --reference="$working_dir/$sample/$locus/bwa/${sample}_${locus}_supercontig.fasta" \
    "$output_whatshap_dir/${sample}_${locus}_phased.bcf" \
    "$working_dir/$sample/$locus/gatk/${sample}_${locus}_marked_duplicates.bam"

    samtools index "$output_whatshap_dir/${sample}_${locus}_haplotagged.bam"
    echo "======================================================="
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

    # 設置raw_fq目錄的路徑
    raw_fq_dir="/home2/yxwu/sequence_raw/Miseq_raw/${project}_FASTQ/trimmed"

    # 設置sample目錄的路徑
    sample_dir="/home/yixuan/Hybpiper/LMRPS_capture/$project/$sample"

    # 遍歷sample目錄下的所有資料夾
    for dir in "$sample_dir"/*; do
        # 檢查是否是資料夾
        if [ -d "$dir" ]; then
            # 執行process_sample_locus命令\
            dir_parts=(${dir//\// })
            locus="${dir_parts[-1]}"
            if [ ! -f "$dir/$sample/sequences/intron/${locus}_supercontig.fasta" ]; then
                echo "no supercontig for locus: $locus"
                continue
            fi
            process_sample_locus "$dir/$sample/sequences/intron/${locus}_supercontig.fasta" "$project" "$ploidy"
        fi
    done
done < "$config_file"
