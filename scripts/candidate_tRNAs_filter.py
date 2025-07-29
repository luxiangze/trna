#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
根据 Rosetta RNA-蛋白质复合物结构预测的 score 结果，
筛选出与所有 aaRs 对接后平均 score 最高的 tRNA
"""

import os
import argparse
import pandas as pd


def extract_trna_id(sample_id):
    """
    从 sample_id 中提取 tRNA 标识符
    例如：A0A8R1WPS3_tRNA-Asn-GTT-2 -> tRNA-Asn-GTT-2
    sample_id 格式为：aaRs_ID_tRNA_ID，用下划线分割
    """
    parts = sample_id.split('_')
    if len(parts) >= 2:
        # 取从第二个下划线开始的所有部分，重新用下划线连接
        # 这样可以处理 tRNA ID 中可能包含下划线的情况
        trna_id = '_'.join(parts[1:])
        return trna_id
    else:
        # 如果格式不符合预期，返回原始 sample_id
        return sample_id


def filter_candidate_trnas(scores_file, out_file, block_list=None):
    """
    筛选候选 tRNA，计算每个 tRNA 与所有蛋白质对接的平均 score
    
    Args:
        scores_file: 输入的 scores CSV 文件
        out_file: 输出的 tRNA scores CSV 文件
        block_list: 需要排除的蛋白 ID 列表文件路径
    """
    print(f"读取 scores 文件: {scores_file}")

    # 读取 scores 文件
    try:
        df = pd.read_csv(scores_file)
        print(f"成功读取 {len(df)} 条记录")
    except Exception as e:
        print(f"读取文件失败: {e}")
        return

    # 检查必要的列是否存在
    required_columns = ['total_score', 'sample_id']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"缺少必要的列: {missing_columns}")
        return
    
    # 如果提供了 block_list，读取需要排除的蛋白 ID
    blocked_proteins = set()
    if block_list and os.path.exists(block_list):
        try:
            with open(block_list, 'r', encoding='utf-8') as f:
                blocked_proteins = set(line.strip() for line in f if line.strip())
            print(f"从 {block_list} 读取到 {len(blocked_proteins)} 个需要排除的蛋白 ID")
        except Exception as e:
            print(f"读取 block_list 文件失败: {e}")
            return
    elif block_list:
        print(f"警告: block_list 文件不存在: {block_list}")
    
    # 过滤掉被阻止的蛋白质样本
    if blocked_proteins:
        original_count = len(df)
        # 提取每个 sample_id 的蛋白质 ID（第一个下划线前的部分）
        df['protein_id'] = df['sample_id'].apply(lambda x: x.split('_')[0])
        # 过滤掉在 block_list 中的蛋白质
        df = df[~df['protein_id'].isin(blocked_proteins)]
        filtered_count = len(df)
        print(f"排除了 {original_count - filtered_count} 条记录，剩余 {filtered_count} 条记录")
        # 删除临时列
        df = df.drop('protein_id', axis=1)

    # 提取 tRNA ID
    df['trna_id'] = df['sample_id'].apply(extract_trna_id)

    print(f"发现 {df['trna_id'].nunique()} 个不同的 tRNA")

    # 按 tRNA 分组，计算平均 score
    trna_stats = df.groupby('trna_id').agg({
        'total_score': ['mean', 'std', 'count', 'min', 'max']
    }).round(3)

    # 展平列名
    trna_stats.columns = ['mean_score', 'std_score', 'count', 'min_score', 'max_score']
    trna_stats = trna_stats.reset_index()

    # 按平均 score 排序（score 越低越好，所以升序排列）
    trna_stats = trna_stats.sort_values('mean_score')

    print(f"\n前10个平均 score 最低（最好）的 tRNA:")
    print(trna_stats.head(10)[['trna_id', 'mean_score', 'count']].to_string(index=False))

    print(f"\n后10个平均 score 最高（最差）的 tRNA:")
    print(trna_stats.tail(10)[['trna_id', 'mean_score', 'count']].to_string(index=False))

    # 创建输出目录（如果不存在）
    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    # 保存结果
    trna_stats.to_csv(out_file, index=False)
    print(f"\n结果已保存到: {out_file}")

    # 输出统计信息
    print(f"\n统计信息:")
    print(f"  总 tRNA 数量: {len(trna_stats)}")
    print(f"  平均对接次数: {trna_stats['count'].mean():.1f}")
    print(f"  最佳平均 score: {trna_stats['mean_score'].min():.3f}")
    print(f"  最差平均 score: {trna_stats['mean_score'].max():.3f}")
    print(f"  平均 score 标准差: {trna_stats['mean_score'].std():.3f}")

    return trna_stats


def main():
    parser = argparse.ArgumentParser(description='筛选候选 tRNA，基于与 aaRs 对接的平均 score')
    parser.add_argument('--scores_file', type=str, required=True, help='输入的 scores CSV 文件')
    parser.add_argument('--out_file', type=str, required=True, help='输出的 tRNA scores CSV 文件')
    parser.add_argument('--block_list', type=str, help='包含需要排除的蛋白 ID 的文件路径')
    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.exists(args.scores_file):
        print(f"错误: 输入文件不存在: {args.scores_file}")
        return

    # 筛选候选 tRNA
    filter_candidate_trnas(args.scores_file, args.out_file, args.block_list)


if __name__ == '__main__':
    main()
