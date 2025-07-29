#!/bin/bash
#SBATCH -J rna_denovo_%a
#SBATCH -o logs/rna_denovo_%A_%a.out
#SBATCH -e logs/rna_denovo_%A_%a.err
#SBATCH -p comput
#SBATCH -t 500:00:00
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --array=1-716  # 实际任务总数

# 设置工作目录
WORK_DIR="$HOME/trna/work/rosetta/Bm_in_Sf"
OUTPUT_DIR="$HOME/trna/work/rosetta/Bm_in_Sf/results"

# 创建输出和日志目录
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# 进入工作目录
cd "$WORK_DIR"

# 获取所有符合条件的fasta文件并创建任务列表
find "$WORK_DIR" -name "*.fasta" -type f | while read fasta_file; do
    # 提取基本文件名（不含.fasta后缀）
    base_name=$(basename "$fasta_file" .fasta)

    # 提取PDB ID（第一个下划线前的部分）
    pdb_id=$(echo "$base_name" | cut -d '_' -f 1)

    # 构建对应的secstruct文件路径
    secstruct_file="${WORK_DIR}/${base_name}.txt"
    pdb_file="${WORK_DIR}/${pdb_id}.pdb"

    # 检查文件是否存在
    if [ -f "$secstruct_file" ] && [ -f "$pdb_file" ]; then
        echo "$fasta_file,$secstruct_file,$pdb_file"
    fi
done > task_list.txt

# 获取任务总数
OFFSET=1000
TOTAL_COUNT=$(( $(wc -l < task_list.txt) - OFFSET ))

# 更新脚本中的任务数组设置
# sed -i "s/TOTAL_COUNT/$TOTAL_COUNT/" "$0"

# 如果是数组作业中的一个任务，则执行对应的命令
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    # 获取当前任务要处理的文件
    # 计算实际任务ID（考虑偏移量）
    ACTUAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))
    TASK_INFO=$(sed -n "${ACTUAL_ID}p" task_list.txt)

    # 解析任务信息
    FASTA_FILE=$(echo "$TASK_INFO" | cut -d',' -f1)
    SECSTRUCT_FILE=$(echo "$TASK_INFO" | cut -d',' -f2)
    PDB_FILE=$(echo "$TASK_INFO" | cut -d',' -f3)

    # 提取基本文件名（用于输出）
    BASE_NAME=$(basename "$FASTA_FILE" .fasta)

    # 创建任务专用目录
    TASK_DIR="$OUTPUT_DIR/$BASE_NAME"
    mkdir -p "$TASK_DIR"

    # 进入任务目录
    cd "$TASK_DIR"

    # 加载所需模块（根据你的集群环境可能需要调整）
    # module load rosetta

    echo "======================================================="
    echo "任务 ${SLURM_ARRAY_TASK_ID}/${TOTAL_COUNT}: 处理 $BASE_NAME"
    echo "FASTA文件: $FASTA_FILE"
    echo "结构文件: $SECSTRUCT_FILE"
    echo "PDB文件: $PDB_FILE"
    echo "输出目录: $TASK_DIR"
    echo "开始时间: $(date)"
    echo "======================================================="

    # 执行rna_denovo命令
    rna_denovo \
        -fasta "$FASTA_FILE" \
        -secstruct_file "$SECSTRUCT_FILE" \
        -s "$PDB_FILE" \
        -minimize_rna false \
        -nstruct 10

    RESULT=$?

    echo "======================================================="
    echo "任务 ${SLURM_ARRAY_TASK_ID}/${TOTAL_COUNT}: $BASE_NAME 完成"
    echo "结束时间: $(date)"
    echo "返回状态: $RESULT"
    echo "======================================================="

    # 返回原始目录
    cd "$WORK_DIR"
else
    # 这部分仅在初次提交时执行，不会在计算节点上运行
    echo "已创建任务列表，共 $TOTAL_COUNT 个任务"
    echo "任务列表保存在: $WORK_DIR/task_list.txt"
    echo "使用命令 'sbatch $0' 提交作业"
fi
