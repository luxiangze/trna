#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
tRNA名称映射文件生成器

基于tRNAscan-SE输出的.ss文件，生成符合命名规则的映射文件。
命名规则：tRNA-氨基酸-反密码子-数字1-数字2
其中数字1和数字2的规则：
- 在相同氨基酸和反密码子基础上：
  - 如果序列相同，第一个数字相同，第二个数字按顺序编号
  - 如果序列不同，第一个数字按序列类型编号
"""

import re
import argparse
from collections import defaultdict

def parse_trnascan_output(input_file):
    """
    解析tRNAscan-SE输出的ss文件，提取tRNA信息

    参数:
        input_file: 输入文件路径

    返回:
        [{"id": "tRNAdb-id", "aa_type": "氨基酸", "anticodon": "反密码子", "sequence": "序列"}]
    """
    trnas = []
    current_trna = None

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()

            # 新的tRNA记录开始
            if re.match(r'^[a-zA-Z0-9]+\.trna\d+', line):
                if current_trna:
                    trnas.append(current_trna)

                trna_id = line.split()[0]
                current_trna = {'id': trna_id}

            # 解析类型和反密码子
            elif "Type:" in line and "Anticodon:" in line:
                parts = line.split()
                aa_type = parts[1]  # 氨基酸类型
                anticodon = parts[3]  # 反密码子
                current_trna['aa_type'] = aa_type
                current_trna['anticodon'] = anticodon

            # 解析序列
            elif line.startswith("Seq:"):
                seq = line[5:].strip()
                current_trna['sequence'] = seq

            # 解析结构
            elif line.startswith("Str:"):
                structure = line[5:].strip()
                current_trna['structure'] = structure

    # 添加最后一个tRNA
    if current_trna:
        trnas.append(current_trna)

    return trnas

def generate_name_map(trnas):
    """
    根据规则生成tRNA名称映射

    参数:
        trnas: tRNA信息列表

    返回:
        映射字典 {原始id: 新id}
    """
    # 按照氨基酸和反密码子分组
    groups = defaultdict(list)
    for trna in trnas:
        key = (trna['aa_type'], trna['anticodon'])
        groups[key].append(trna)

    # 生成映射
    name_map = {}

    for (aa_type, anticodon), group in groups.items():
        # 按照序列进一步分组
        sequence_groups = defaultdict(list)
        for trna in group:
            sequence_groups[trna['sequence']].append(trna)

        # 分配第一个数字（按照不同序列）
        seq_number = 1
        for sequence, seq_group in sequence_groups.items():
            # 对每个序列组，分配第二个数字
            for i, trna in enumerate(seq_group, 1):
                new_id = f"tRNA-{aa_type}-{anticodon}-{seq_number}-{i}"
                name_map[trna['id']] = new_id

            seq_number += 1

    return name_map

def write_name_map(name_map, output_file):
    """
    将名称映射写入文件

    参数:
        name_map: 映射字典
        output_file: 输出文件路径
    """
    with open(output_file, 'w') as f:
        # 写入表头
        f.write("tRNAscan-SE_id\tGtRNAdb_id\n")

        # 自定义排序函数，提取染色体编号和tRNA编号
        def sort_key(trna_id):
            # 提取染色体编号
            chr_match = re.search(r'([a-zA-Z]+)(\d+)', trna_id)
            if chr_match:
                chr_prefix = chr_match.group(1)  # 如 'chr'
                chr_num = int(chr_match.group(2)) if chr_match.group(2) else 0
            else:
                chr_prefix = ''
                chr_num = 0

            # 提取tRNA编号
            trna_match = re.search(r'trna(\d+)', trna_id)
            trna_num = int(trna_match.group(1)) if trna_match else 0

            return (chr_prefix, chr_num, trna_num)

        # 按原始id排序并写入映射
        for trna_id in sorted(name_map.keys(), key=sort_key):
            f.write(f"{trna_id}\t{name_map[trna_id]}\n")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='生成tRNA名称映射文件')
    parser.add_argument('-i', '--input', required=True, help='tRNAscan-SE输出的.ss文件')
    parser.add_argument('-o', '--output', required=True, help='输出映射文件路径')

    args = parser.parse_args()

    # 解析tRNAscan-SE输出
    trnas = parse_trnascan_output(args.input)

    # 生成名称映射
    name_map = generate_name_map(trnas)

    # 写入映射文件
    write_name_map(name_map, args.output)

    print(f"成功处理 {len(trnas)} 条tRNA记录")
    print(f"映射文件已保存至: {args.output}")

if __name__ == "__main__":
    main()
