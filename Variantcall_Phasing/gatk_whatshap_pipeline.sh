#!/bin/bash

working_dir="/home2/yxwu/Phasing"
raw_fq_dir="/home2/yxwu/sequence_raw/Miseq_raw/MS23395_FASTQ/trimmed"
ploidy=1

# 定義函數
check_and_create_dir() {
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

}

process_sample_locus() {
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

    check_and_create_dir "$working_dir/$sample"
    
    copy_and_index_reference "$input_file" "$sample" "$locus"
    bwa_mapping_reads "$RG_ID" "$sample" "$locus"
    gatk_mark_duplicates "$sample" "$locus"
    gatk_call_variants "$ploidy" "$sample" "$locus"
    gatk_filter_variants "$sample" "$locus"
    whatshap_phase_variants "$sample" "$locus"
    whatshap_visualization "$sample" "$locus"
}


copy_reference() {
    local input_file="$1"
    local sample="$2"
    local locus="$3"

    # 構建輸出目錄路徑
    output_dir="$working_dir/$sample/$locus/ref"

    # 創建目標目錄(如果不存在)
    mkdir -p "$output_dir"

    # 將 ref.fas 複製到目標目錄
    cp "$input_file" "$output_dir"/"${sample}_${locus}_supercontig.fasta"
}

bwa_mapping_reads() {
    local RG_ID="$1"
    local sample="$2"
    local locus="$3"
    # 將讀段與參考基因組(supercontig.fasta)進行比對mapping，只存mapping到的reads(samtools -F 4)
    # 從ref+trimmed_fastq.gz>_mapped.sam

    # 構建輸出目錄路徑
    output_dir="$working_dir/$sample/$locus/ref"
    # 進入目標目錄並索引ref
    cd "$output_dir"
    bwa index -p "${sample}_${locus}" "${sample}_${locus}_supercontig.fasta"

    # 創建目標目錄(如果不存在)
    mkdir -p "$working_dir/${sample}_bwa"
    mkdir -p "$working_dir/${sample}_bwa/output"

    #bwa | samtools
	bwa mem -M -t 20 -R "@RG\tID:${RG_ID}\tSM:${sample}_${locus}" \
	"$working_dir/${sample}_bwa/${sample}_${locus}" \
	"$raw_fq_dir/${sample}*R1_trimmed.fastq.gz" "$raw_fq_dir/${sample}*R2_trimmed.fastq.gz" \
	2> "$working_dir/${sample}_bwa/bwa.err" \
	| samtools view -h -F 4 -o "$working_dir/${sample}_bwa/output/${locus}/${sample}_${locus}_mapped.sam"

}

gatk_mark_duplicates() {
    # 標記重複讀段，從_mapped.sam>_marked_duplicates.bam
    local sample="$1"
    local locus="$2"
	
    mkdir -p "$working_dir/${sample}_gatk"

    gatk MarkDuplicatesSpark \
    -I "$working_dir/${sample}_bwa/output/${locus}/${sample}_${locus}_mapped.sam" \
    -O "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_duplicates.bam" \
    -M "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_dup_metrics.txt"
}

gatk_call_variants() {
    # 呼叫變異位點，從_marked_duplicates.bam>.vcf
    local ploidy="$1"
    local sample="$2"
    local locus="$3"
	gatk CreateSequenceDictionary -R "$working_dir/$sample/$locus/ref/${sample}_${locus}_supercontig.fasta"
    
	samtools faidx "$working_dir/$sample/$locus/ref/${sample}_${locus}_supercontig.fasta"
    
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -ploidy $ploidy \
    -R "$working_dir/$sample/$locus/ref/${sample}_${locus}_supercontig.fasta" \
    -I "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_duplicates.bam" \
    -O "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}.vcf"

}

gatk_filter_variants() {
    # 從.vcf>_filtered.vcf
    local sample="$1"
    local locus="$2"

    # 過濾變異位點
	gatk VariantFiltration \
    --variant "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}.vcf" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --filter-expression "AF < 0.05 || AF > 0.95" --filter-name "AF5-95" \
    --output "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_filtered.vcf"
}

whatshap_phase_variants() {
    # 從_filtered.vcf+_marked_duplicates.bam>_phased.bcf
    local sample="$1"
    local locus="$2"

    cp -r "$working_dir/$sample/$locus/ref/." "$working_dir/${sample}_gatk/${locus}"
    #Index the variants file (skip if a .csi or .tbi file already exists):
    bcftools view "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_filtered.vcf" -Oz -o "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_filtered.vcf.gz"
    bcftools index "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_filtered.vcf.gz"
    
    #Index the alignment file (skip if a .csi or .bai file already exists):
    samtools index "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_duplicates.bam"
	
	cd "$working_dir/${sample}_gatk/${locus}"

    # 創建目標目錄(如果不存在)
    mkdir -p "$working_dir/${sample}_whatshap/${locus}"

	#Run read-based phasing.
    whatshap phase \
    --reference="$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_supercontig.fasta" \
    -o "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.bcf" \
    "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_filtered.vcf" "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_duplicates.bam"

}

whatshap_visualization() {
    local sample="$1"
    local locus="$2"

    cd "$working_dir/${sample}_whatshap/${locus}"
    # 創建目標目錄(如果不存在)
    mkdir -p "$working_dir/${sample}_whatshap/${locus}/visualization"

    # 用IGV視覺化前置index
    bcftools convert -O v -o "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.vcf" "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.bcf"
    whatshap stats --gtf="$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.gtf" "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.bcf"
    bcftools index "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.bcf"

    # 用IGV視覺化
    whatshap haplotag -o "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_haplotagged.bam" \
    --output-haplotag-list="$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_haplotag_table.tsv" \
    --reference="$working_dir/$sample/$locus/ref/${sample}_${locus}_supercontig.fasta" \
    "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_phased.bcf" \
    "$working_dir/${sample}_gatk/${locus}/${sample}_${locus}_marked_duplicates.bam"

    samtools index "$working_dir/${sample}_whatshap/${locus}/${sample}_${locus}_haplotagged.bam"
}

# 主程式
process_sample_locus /home/yixuan/Hybpiper/LMRPS_capture/MS23395/YXD136/transducin/YXD136/sequences/intron/transducin_supercontig.fasta MS23395 1
