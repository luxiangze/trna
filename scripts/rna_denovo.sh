#!/bin/bash

# 批量运行rna_denovo命令处理输入文件夹下的所有文件
# 用法: bash scripts/rna_denovo.sh <输入目录> <PDB目录>
# 示例: bash scripts/rna_denovo.sh work/rosetta/Bm_in_Bm work/colabfold/Bm_aaRSs_output_pdb

# 检查参数
if [ $# -ne 2 ]; then
    echo "用法: bash scripts/rna_denovo.sh <输入目录> <PDB目录>"
    echo "示例: bash scripts/rna_denovo.sh work/rosetta/Bm_in_Bm work/colabfold/Bm_aaRSs_output_pdb"
    exit 1
fi

# 设置输入和输出目录
INPUT_DIR="$1"
PDB_DIR="$2"
OUTPUT_DIR="$INPUT_DIR/results"

# 检查输入目录是否存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "错误: 输入目录 '$INPUT_DIR' 不存在"
    exit 1
fi

if [ ! -d "$PDB_DIR" ]; then
    echo "错误: PDB目录 '$PDB_DIR' 不存在"
    exit 1
fi

# 创建输出目录和日志目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$INPUT_DIR/logs"

echo "======================================================="
echo "RNA Denovo 批量处理脚本"
echo "输入目录: $INPUT_DIR"
echo "PDB目录: $PDB_DIR"
echo "输出目录: $OUTPUT_DIR"
echo "开始时间: $(date)"
echo "======================================================="

# 定义处理单个样本的函数
process_sample() {
    local fasta_file="$1"
    local input_dir="$2"
    local pdb_dir="$3"
    local output_dir="$4"

    # 提取基本文件名（不含.fasta后缀）
    local base_name=$(basename "$fasta_file" .fasta)

    # 提取PDB ID（第一个下划线前的部分）
    local pdb_id=$(echo "$base_name" | cut -d '_' -f 1)

    # 构建对应的文件路径
    local secstruct_file="$input_dir/${base_name}.txt"
    local pdb_file="$pdb_dir/${pdb_id}.pdb"

    # 创建样本专用输出目录
    local sample_dir="$output_dir/$base_name"
    mkdir -p "$sample_dir"

    # 确保日志目录存在（在并行环境中很重要）
    mkdir -p "$input_dir/logs"

    # 创建日志文件路径
    local log_file="$input_dir/logs/${base_name}.log"

    # 初始化日志文件
    touch "$log_file"

    # 日志函数
    log_message() {
        local message="$1"
        # 确保日志文件存在
        if [ ! -f "$log_file" ]; then
            touch "$log_file"
        fi
        echo "[$(date)] $message" | tee -a "$log_file"
    }

    # 开始记录日志
    log_message "======================================================="
    log_message "开始处理样本: $base_name"
    log_message "FASTA文件: $fasta_file"
    log_message "结构文件: $secstruct_file"
    log_message "PDB文件: $pdb_file"
    log_message "输出目录: $sample_dir"
    log_message "日志文件: $log_file"
    log_message "======================================================="

    # 检查必要文件是否存在
    if [ ! -f "$secstruct_file" ]; then
        log_message "错误: 缺少结构文件 $secstruct_file，跳过处理"
        return 1
    fi

    if [ ! -f "$pdb_file" ]; then
        log_message "错误: 缺少PDB文件 $pdb_file，跳过处理"
        return 1
    fi

    log_message "所有必要文件检查通过，开始执行rna_denovo命令"

    # 进入样本目录执行命令
    cd "$sample_dir"

    # 构建完整的命令用于日志记录
    local cmd="rna_denovo -fasta $fasta_file -secstruct_file $secstruct_file -s $pdb_file -minimize_rna false -nstruct 10"
    log_message "执行命令: $cmd"

    # 执行rna_denovo命令，将标准输出和错误输出都重定向到日志文件
    {
        echo "[$(date)] === RNA Denovo 输出开始 ==="
        rna_denovo \
            -fasta "$fasta_file" \
            -secstruct_file "$secstruct_file" \
            -s "$pdb_file" \
            -minimize_rna false \
            -nstruct 10
        echo "[$(date)] === RNA Denovo 输出结束 ==="
    } >> "$log_file" 2>&1

    local result=$?

    if [ $result -eq 0 ]; then
        log_message "处理成功完成，返回码: $result"
        log_message "输出文件保存在: $sample_dir"
    else
        log_message "处理失败，返回码: $result"
    fi

    log_message "样本 $base_name 处理完成"
    log_message "======================================================="

    return $result
}

# 导出函数以供parallel使用
export -f process_sample

# 查找所有fasta文件并使用parallel处理
echo "正在查找FASTA文件..."
fasta_files=$(find "$INPUT_DIR" -name "*.fasta" -type f)
file_count=$(echo "$fasta_files" | wc -l)

if [ -z "$fasta_files" ]; then
    echo "错误: 在 '$INPUT_DIR' 中未找到任何FASTA文件"
    exit 1
fi

echo "找到 $file_count 个FASTA文件"
echo "使用 $(nproc) 个CPU核心进行并行处理"
echo "======================================================="

# 使用GNU parallel并行处理所有样本
echo "$fasta_files" | parallel -j $(nproc) process_sample {} "$INPUT_DIR" "$PDB_DIR" "$OUTPUT_DIR"

echo "======================================================="
echo "所有文件处理完毕！"
echo "结束时间: $(date)"
echo "结果保存在: $OUTPUT_DIR"
echo "======================================================="
