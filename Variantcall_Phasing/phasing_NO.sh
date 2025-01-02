#!/bin/bash

# 定義函數
is_single_copy() {
    echo "================= 檢查是否是single_copy =================="
    local single_copy_list=("ApPEFP_C" "OG0000764" "OG0000806" "OG0001134" "OG0001287" "OG0001587" "OG0003109" "OG0003385" "OG0003734" "OG0003748" "OG0003780" "OG0004009" "OG0004248" "OG0004360" "OG0004561" "OG0004566" "OG0004614" "OG0004903" "OG0004943" "OG0004969" "OG0005419" "OG0005470" "OG0005474" "OG0005881" "OG0006029" "OG0006334" "OG0006337" "OG0006377" "OG0006406" "OG0006606" "OG0006650" "OG0006781" "OG0006950" "OG0006968" "OG0007034" "OG0007059" "OG0007255" "OG0007337" "OG0007419" "OG0007472" "OG0007491" "OG0007526" "OG0007580" "OG0007584" "OG0007607" "OG0007648" "OG0007814" "OG0007872" "OG0008201" "OG0008464" "OG0008477" "OG0008526" "OG0008551" "OG0008618" "OG0008661" "OG0008757" "OG0008877" "OG0008889" "OG0008981" "OG0008997" "OG0009191" "OG0009405" "OG0009430" "OG0009437" "OG0009708" "OG0009761" "OG0009864" "OG0009975" "OG0010004" "OG0010124" "OG0010148" "OG0010266" "OG0010273" "OG0010275" "OG0010383" "OG0010439" "OG0010506" "OG0010519" "OG0010559" "OG0010562" "OG0010691" "OG0010995" "OG0011134" "OG0011178" "OG0011263" "OG0011365" "OG0011418" "OG0011530" "OG0011599" "OG0011616" "OG0011638" "OG0011683" "OG0011739" "OG0011740" "OG0011920" "OG0011936" "OG0011960" "OG0011994" "OG0012575" "OG0012619" "OG0012642" "OG0012644" "OG0012674" "OG0012719" "OG0012796" "OG0012887" "OG0012889" "OG0012896" "OG0012930" "OG0013120" "OG0013744" "OG0013794" "OG0013810" "OG0013859" "OG0013866" "OG0013892" "OG0014169" "OG0014204" "OG0014205" "OG0014282" "OG0015475")
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


copy_and_index_reference() {
    echo "============= 複製並索引參考基因組序列 ============="
    local input_file="$1"
    local sample="$2"
    local locus="$3"

    # 構建輸出目錄路徑
    output_bwa_dir="$working_dir/$sample/$locus/bwa"

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_bwa_dir"

    # 將 ref.fas 複製到目標目錄並把fasta header加上locus name
    sed "s/${sample}/${sample}_${locus}/" "$input_file" > "$output_bwa_dir"/"${sample}_${locus}_supercontig.fasta"

    # 進入目標目錄並索引ref
    cd "$output_bwa_dir"
    bwa index -p "${sample}_${locus}" "${sample}_${locus}_supercontig.fasta"

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

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_bwa_dir"
    mkdir -p "$output_bwa_dir/output"

    #bwa | samtools
    bwa mem -M -t 30 -R "@RG\tID:${RG_ID}\tSM:${sample}_${locus}" \
    "$output_bwa_dir/${sample}_${locus}" \
    "$raw_fq_dir/${sample}"*"R1_trimmed.fastq.gz" "$raw_fq_dir/${sample}"*"R2_trimmed.fastq.gz" \
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
    --read-filter SoftClippedReadFilter --invert-soft-clip-ratio-filter --soft-clipped-leading-trailing-ratio 0.3 \
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
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "DP < 5" --filter-name "DP5" \
    --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" \
    --output "$output_gatk_dir/${sample}_${locus}_filtered.vcf"

    bcftools view -f PASS "$output_gatk_dir/${sample}_${locus}_filtered.vcf" > "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf"
    echo "======================================================="
}

whatshap_phase_variants() {
    echo "============= 對變異位點進行分型 ============="
    # 從_filtered_PASS.vcf+_marked_duplicates.bam>_phased.bcf
    local ploidy="$1"
    local sample="$2"
    local locus="$3"

    cd "$working_dir/$sample/$locus"
    # 構建輸出目錄路徑
    output_gatk_dir="$working_dir/$sample/$locus/gatk"

    cp -r "$working_dir/$sample/$locus/bwa/." "$output_gatk_dir"
    #Index the variants file (skip if a .csi or .tbi file already exists):
    bcftools view "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf" -Oz -o "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf.gz"
    bcftools index "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf.gz"
    
    #Index the alignment file (skip if a .csi or .bai file already exists):
    samtools index "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam"
    
    cd "$output_gatk_dir"

    # 構建輸出目錄路徑
    output_whatshap_dir="$working_dir/$sample/$locus/whatshap"
    # 創建目標目錄(如果不存在)
    mkdir -p "$output_whatshap_dir"

    if [ $ploidy -eq 2 ];then
        #Run read-based phasing.
        whatshap phase \
        --reference="$output_gatk_dir/${sample}_${locus}_supercontig.fasta" \
        -o "$output_whatshap_dir/${sample}_${locus}_phased.bcf" \
        "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf" "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam"
        echo "======================================================="
    else
        #Run read-based polyphasing.
        whatshap polyphase \
        --ploidy $ploidy --threads 8 \
        --reference "$output_gatk_dir/${sample}_${locus}_supercontig.fasta" \
        -o "$output_whatshap_dir/${sample}_${locus}_phased.bcf" \
        "$output_gatk_dir/${sample}_${locus}_filtered_PASS.vcf" "$output_gatk_dir/${sample}_${locus}_marked_duplicates.bam"
        echo "======================================================="
    fi
}

whatshap_visualization() {
    echo "============= whatshap 可視化 ============="
    local sample="$1"
    local locus="$2"


    # 構建輸出目錄路徑
    output_whatshap_dir="$working_dir/$sample/$locus/whatshap"

    cd "$output_whatshap_dir"


    # 用IGV視覺化前置index
    bcftools index -c -f "$output_whatshap_dir/${sample}_${locus}_phased.bcf"
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

bcftools_concensus() {
    echo "============= bcftools將bcf的haplotype生成出來 ============="
    local input_file="$1"
    local sample="$2"
    local locus="$3"
    local ploidy="$4"

    whatshap_dir="$working_dir/$sample/$locus/whatshap"
    output_dir="$working_dir/$sample/$locus/phased"


    # 檢查bcf文件是否有任何variant
    num_variants=$(bcftools view -H "$whatshap_dir/${sample}_${locus}_phased.bcf" | wc -l)

    if [ "$num_variants" -eq 0 ]; then
        echo "No variants found in $whatshap_dir/${sample}_${locus}_phased.bcf. Skipping..."
        return
    fi

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_dir"

    for ((i=1; i<=$ploidy; i++))
    do
        echo "generating haplotype${i} of ${sample}_${locus}"
        bcftools consensus -H $i -f $input_file \
        "$whatshap_dir/${sample}_${locus}_phased.bcf" \
        | sed "s/>.*$/>${sample}_${locus}_haplotype${i}/" \
        > "${output_dir}/${sample}_${locus}_haplotype${i}.fasta"
    done
    
    cat "${output_dir}/${sample}_${locus}_haplotype"*".fasta" > "${output_dir}/pooled_${sample}_${locus}_haplotypes.fasta"

    mafft "${output_dir}/pooled_${sample}_${locus}_haplotypes.fasta" > "${output_dir}/aligned_pooled_${sample}_${locus}_haplotypes.fasta"
    mafft --auto --add $input_file "${output_dir}/aligned_pooled_${sample}_${locus}_haplotypes.fasta" > "${output_dir}/aligned_pooled_${sample}_${locus}_haplotypes_withref.fasta"
}

# 主程式

# 檢查是否提供了正確的參數
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

working_dir="/home2/yxwu/Phasing_No"


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
    raw_fq_dir="/home2/yxwu/sequence_raw/Novaseq_raw/${project}_FASTQ/trimmed"

    # 設置sample目錄的路徑
    sample_dir="/home/yixuan/Hybpiper/LMRPS_capture/$project/$sample"

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
                bcftools_concensus "$working_dir/$sample/$locus/bwa/${sample}_${locus}_supercontig.fasta" "$sample" "$locus" "$ploidy"
            fi
        fi
    done
done < "$config_file"
