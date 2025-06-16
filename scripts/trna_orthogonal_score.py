from dataclasses import dataclass
from Bio import Align
import argparse
import pandas as pd
import numpy as np

@dataclass
class tRNARecord:
    tRNA_id: str
    seq: str
    anticodon: str
    amino_acid_type: str
    identity_elements: list[int]
    indices: list[int]

def phase_identity_elements(identity_elements_file) -> {str: list[int]}:
    """读取identity_elements.txt文件, 解析identity elements"""
    identity_elements_dict = {}
    with open(identity_elements_file, "r") as f:
        for line in f:
            line = line.strip()
            identity_elements_dict[line.split("\t")[0]] = list(map(int, line.split("\t")[1].split(", ")))
    return identity_elements_dict

def parse_alignment(alignment, identity_elements_dict) -> list[tRNARecord]:
    """解析alignment"""
    tRNARecords = []
    for l in range(len(alignment)):
        tRNA_id = alignment.sequences[l].id
        if "iMet" in tRNA_id:
            amino_acid_type = "Met"
        else:
            amino_acid_type = tRNA_id.split("-")[1]
        tRNARecords.append(
            tRNARecord(
                tRNA_id = tRNA_id,
                seq = alignment[l],
                anticodon = tRNA_id.split("-")[2],
                amino_acid_type = amino_acid_type,
                identity_elements = identity_elements_dict[amino_acid_type + "RS"],
                indices = alignment.indices[l]
                )
            )
    return tRNARecords

def calculate_orthogonal_score(query_tRNARecords, target_tRNARecords) -> pd.DataFrame:
    """
    计算tRNA的正交得分
    根据文献方法：
    1. 比较每个候选tRNA和E. coli的各个tRNA异构体
    2. 在关键元素位置上，匹配得+1分，不匹配得-1分
    3. 为每对比较计算平均得分
    4. 对同一种氨基酸类型的所有异构体计算平均得分
    5. 过滤掉与其自身氨基酸类型的E. coli tRNA得分>0的tRNA
    """
    # 创建结果DataFrame
    result_data = []

    # 对每个查询tRNA进行处理
    for query_tRNA in query_tRNARecords:
        # 按照目标氨基酸类型分组目标tRNA
        target_by_aa = {}
        for target_tRNA in target_tRNARecords:
            if target_tRNA.amino_acid_type not in target_by_aa:
                target_by_aa[target_tRNA.amino_acid_type] = []
            target_by_aa[target_tRNA.amino_acid_type].append(target_tRNA)

        # 对每种氨基酸类型计算正交得分
        for aa_type, targets in target_by_aa.items():
            isoacceptor_scores = []

            # 计算与每个异构体tRNA的得分
            for target_tRNA in targets:
                fixed_identity_elements = [x - 1 for x in target_tRNA.identity_elements]
                individual_scores = []

                # 在每个身份元素位置比较并计算得分
                for identity_element in fixed_identity_elements:
                    query_idx_arr = np.where(query_tRNA.indices == identity_element)[0]
                    target_idx_arr = np.where(target_tRNA.indices == identity_element)[0]
                    if len(query_idx_arr) == 0 or len(target_idx_arr) == 0:
                        continue
                    query_index = query_idx_arr[0]
                    target_index = target_idx_arr[0]
                    query_pos = query_tRNA.seq[query_index]
                    target_pos = target_tRNA.seq[target_index]
                    # 匹配得+1，不匹配得-1
                    if query_pos == target_pos and query_pos != "-" and target_pos != "-":
                        individual_scores.append(1)
                    elif query_pos != target_pos and query_pos != "-" and target_pos != "-":
                        individual_scores.append(-1)
                    else:
                        continue

                # 计算该tRNA的平均得分
                if individual_scores:
                    avg_score = sum(individual_scores) / len(individual_scores)
                    isoacceptor_scores.append(avg_score)

            # 计算该氨基酸类型的平均得分
            if isoacceptor_scores:
                aa_score = sum(isoacceptor_scores) / len(isoacceptor_scores)

                # 添加到结果中
                result_data.append({
                    "query_id": query_tRNA.tRNA_id,
                    "query_amino_acid_type": query_tRNA.amino_acid_type,
                    "target_amino_acid_type": aa_type,
                    "orthogonal_score": aa_score
                })

    # 创建结果DataFrame
    orthogonal_score = pd.DataFrame(result_data)
    return orthogonal_score

# def find_positions(source: list[int], target: list[int]) -> list[int]:
#     pos_map = {v: i for i, v in enumerate(target)}
#     return [pos_map[x] for x in source if x in pos_map]

def filter_orthogonal_scores(orthogonal_scores_wide: pd.DataFrame) -> pd.DataFrame:
    """筛选出自身氨基酸类型得分小于0的tRNA"""
    filtered_df = pd.DataFrame()
    for index, row in orthogonal_scores_wide.iterrows():
        tRNA_type = index.split("-")[1]
        if tRNA_type == "iMet":
            tRNA_type = "Met"
        if row[tRNA_type] < 0:
            filtered_df = pd.concat([filtered_df, row.to_frame().T])
    return filtered_df

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="计算tRNA的正交得分")
    parser.add_argument("-q", "--query", required=True, help="查询tRNA的Stockholm文件路径")
    parser.add_argument("-t", "--target", required=True, help="目标物种的tRNA的Stockholm文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出文件路径")
    parser.add_argument("-e", "--identity_elements", required=True, help="identity_elements.txt文件路径")
    args = parser.parse_args()

    query_stk_file = args.query
    target_stk_file = args.target
    identity_elements_file = args.identity_elements
    output_dir = args.output

    # 读取identity_elements.txt文件
    identity_elements = phase_identity_elements(identity_elements_file)

    # 读取Stockholm文件
    query_alignments = Align.parse(query_stk_file, "stockholm")
    target_alignments = Align.parse(target_stk_file, "stockholm")
    query_alignment = next(query_alignments)
    target_alignment = next(target_alignments)

    # 解析alignment
    query_tRNARecords = parse_alignment(query_alignment, identity_elements)
    target_tRNARecords = parse_alignment(target_alignment, identity_elements)

    # 计算正交得分
    orthogonal_scores = calculate_orthogonal_score(query_tRNARecords, target_tRNARecords)
    orthogonal_scores_wide = orthogonal_scores.pivot(index='query_id', columns='target_amino_acid_type', values='orthogonal_score')

    # 输出正交得分
    orthogonal_scores_wide.to_csv(output_dir + "/orthogonal_scores.csv")

    # 筛选候选tRNA - 筛选出tRNA所属氨基酸类型的正交得分小于0的tRNA。
    candidate_tRNAs = filter_orthogonal_scores(orthogonal_scores_wide)

    # 输出结果
    candidate_tRNAs.to_csv(output_dir + "/candidate_tRNAs.csv")

if __name__ == "__main__":
    main()