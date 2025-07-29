#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
设计 tRNA 反密码子的突变体

用法：
python design_mutant_library.py --input_ID <tRNA ID> --input_ID_map <tRNA ID 映射文件> --input_structure <tRNA 结构文件> --out_file <输出FASTA文件路径>
'''

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass, field
from generate_trna_name_map import parse_trnascan_output
from logger_utils import setup_logger, get_logger

@dataclass
class tRNARecord:
    tRNA_name: str
    tRNA_id: str
    seq: str
    anticodon: str
    amino_acid_type: str
    structure: str
    mutants: list[str] = field(default_factory=list)

def parse_arguments():
    '''
    解析命令行参数
    '''
    parser = argparse.ArgumentParser(description='设计tRNA反密码子的突变体')
    parser.add_argument('--input_IDs_file', required=True, help='输入 tRNA ID 文件')
    parser.add_argument('--input_ID_map', required=True, help='输入 tRNA ID 映射文件')
    parser.add_argument('--input_structure', required=True, help='输入 tRNA 结构文件')
    parser.add_argument('--out_file', required=True, help='输出FASTA文件路径')
    return parser.parse_args()

def read_tRNA_id(input_IDs_file, input_ID_map):
    logger = get_logger('design_mutant_library')
    logger.info(f"开始读取tRNA ID，候选文件: {input_IDs_file}，映射文件: {input_ID_map}")

    tRNA_ids = []
    tRNA_ids_dict = {}

    # 读取候选tRNA ID
    logger.info("读取候选tRNA ID...")
    with open(input_IDs_file, 'r') as f:
        for line in f:
            line = line.strip()  # 移除换行符和空格
            if line.startswith('tRNA'):
                tRNA_id = line.split(',')[0]  # 支持逗号分隔格式
                tRNA_ids.append(tRNA_id)
                logger.debug(f"找到候选tRNA: {tRNA_id}")

    logger.info(f"共读取到 {len(tRNA_ids)} 个候选tRNA ID")

    # 读取映射文件并进行匹配
    logger.info("读取映射文件并进行ID匹配...")
    with open(input_ID_map, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                trnascan_id = parts[0]
                gtrnadb_id = parts[1]

                for tRNA_id in tRNA_ids:
                    target_pattern = tRNA_id + "-1"
                    if str(target_pattern) == str(gtrnadb_id):
                        tRNA_ids_dict[tRNA_id] = trnascan_id
                        logger.debug(f"ID匹配成功: {tRNA_id} -> {trnascan_id}")
                        break

    logger.info(f"ID映射完成: {len(tRNA_ids_dict)}/{len(tRNA_ids)} 个tRNA找到了映射")

    # 记录未找到映射的tRNA
    missing_ids = [tid for tid in tRNA_ids if tid not in tRNA_ids_dict]
    if missing_ids:
        logger.warning(f"以下 {len(missing_ids)} 个tRNA未找到映射:")
        for mid in missing_ids:
            logger.warning(f"  未映射: {mid}")

    return tRNA_ids_dict

def tRNA_prepare(input_IDs_file, input_ID_map, input_structure) -> list[tRNARecord]:
    logger = get_logger('design_mutant_library')
    logger.info(f"开始准备tRNA数据，结构文件: {input_structure}")

    tRNA_ids_dict = read_tRNA_id(input_IDs_file, input_ID_map)

    logger.info("解析tRNAscan-SE结构文件...")
    trnas = parse_trnascan_output(input_structure)
    logger.info(f"从结构文件中解析到 {len(trnas)} 个tRNA条目")

    tRNA_records = []
    found_count = 0

    for tRNA_id, trna_name in tRNA_ids_dict.items():
        logger.debug(f"查找tRNA {tRNA_id} (tRNAscan ID: {trna_name})的结构信息")

        structure_found = False
        for trna in trnas:
            if trna_name == trna['id']:
                tRNA_records.append(tRNARecord(trna_name, tRNA_id, trna['sequence'], trna['anticodon'], trna['aa_type'], trna['structure'], []))
                found_count += 1
                structure_found = True
                logger.info(f"找到结构: {tRNA_id} -> 反密码子: {trna['anticodon']}, 氨基酸: {trna['aa_type']}, 序列长度: {len(trna['sequence'])}")
                break

        if not structure_found:
            logger.warning(f"未找到tRNA {tRNA_id} (tRNAscan ID: {trna_name})的结构信息")

    logger.info(f"tRNA数据准备完成: {found_count}/{len(tRNA_ids_dict)} 个tRNA找到了完整信息")

    # 记录未找到结构的tRNA
    missing_structures = [tid for tid in tRNA_ids_dict.keys() if tid not in [tr.tRNA_id for tr in tRNA_records]]
    if missing_structures:
        logger.warning(f"以下 {len(missing_structures)} 个tRNA未找到结构信息:")
        for mid in missing_structures:
            logger.warning(f"  缺失结构: {mid}")

    return tRNA_records

# 更准确的反密码子到氨基酸映射（基于标准遗传密码）
STANDARD_ANTICODON_TO_AA = {
    # Phenylalanine (F) - UUU, UUC
    'AAA': 'Phe', 'GAA': 'Phe',
    # Leucine (L) - UUA, UUG, CUU, CUC, CUA, CUG
    'UAA': 'Leu','UAG': 'Leu', 'CAA': 'Leu', 'AAG': 'Leu', 'GAG': 'Leu', 'CAG': 'Leu',
    # Serine (S) - UCU, UCC, UCA, UCG
    'AGA': 'Ser', 'GGA': 'Ser', 'UGA': 'Ser', 'CGA': 'Ser',
    # Tyrosine (Y) - UAU, UAC
    'AUA': 'Tyr', 'GUA': 'Tyr',
    # Stop - UAA, UAG, UGA
    'UUA': '*', 'CUA': '*', 'UCA': '*',
    # Cysteine (C) - UGU, UGC
    'ACA': 'Cys', 'GCA': 'Cys',
    # Tryptophan (W) - UGG
    'CCA': 'Trp',
    # Proline (P) - CCU, CCC, CCA, CCG
    'AGG': 'Pro', 'GGG': 'Pro', 'UGG': 'Pro', 'CGG': 'Pro',
    # Histidine (H) - CAU, CAC
    'AUG': 'His', 'GUG': 'His',
    # Glutamine (Q) - CAA, CAG
    'UUG': 'Gln', 'CUG': 'Gln',
    # Arginine (R) - CGU, CGC, CGA, CGG, AGA, AGG
    'ACG': 'Arg', 'GCG': 'Arg', 'UCG': 'Arg', 'CCG': 'Arg', 'UCU': 'Arg', 'CCU': 'Arg',
    # Isoleucine (I) - AUU, AUC, AUA
    'AAU': 'Ile', 'GAU': 'Ile', 'UAU': 'Ile',
    # Methionine (M) - AUG
    'CAU': 'Met',
    # Threonine (T) - ACU, ACC, ACA, ACG
    'AGU': 'Thr', 'GGU': 'Thr', 'UGU': 'Thr', 'CGU': 'Thr',
    # Asparagine (N) - AAU, AAC
    'AUU': 'Asn', 'GUU': 'Asn',
    # Lysine (K) - AAA, AAG
    'UUU': 'Lys', 'CUU': 'Lys',
    # Serine (S) - AGU, AGC
    'ACU': 'Ser', 'GCU': 'Ser',
    # Valine (V) - GUU, GUC, GUA, GUG
    'AAC': 'Val', 'GAC': 'Val', 'UAC': 'Val', 'CAC': 'Val',
    # Alanine (A) - GCU, GCC, GCA, GCG
    'AGC': 'Ala', 'GGC': 'Ala', 'UGC': 'Ala', 'CGC': 'Ala',
    # Aspartic acid (D) - GAU, GAC
    'AUC': 'Asp', 'GUC': 'Asp',
    # Glutamic acid (E) - GAA, GAG
    'UUC': 'Glu', 'CUC': 'Glu',
    # Glycine (G) - GGU, GGC, GGA, GGG
    'ACC': 'Gly', 'GCC': 'Gly', 'UCC': 'Gly', 'CCC': 'Gly'
}

def get_amino_acid_from_anticodon(anticodon):
    """根据反密码子获取对应的氨基酸"""
    # 清理反密码子格式
    clean_anticodon = anticodon.upper().replace('T', 'U')
    result = STANDARD_ANTICODON_TO_AA.get(clean_anticodon, 'X')  # X表示未知
    return result

def generate_alternative_anticodons(original_anticodon, original_aa):
    """生成不同氨基酸的反密码子突变体"""
    logger = get_logger('design_mutant_library')
    logger.debug(f"为反密码子 {original_anticodon} (氨基酸: {original_aa}) 生成突变体")

    bases = ['A', 'U', 'G', 'C']
    alternatives = []

    # 生成所有可能的单碱基突变
    logger.debug("尝试单碱基突变...")
    for i in range(3):  # 反密码子有3个位置
        for base in bases:
            if base != original_anticodon[i]:
                mutant = list(original_anticodon)
                mutant[i] = base
                mutant_anticodon = ''.join(mutant)
                mutant_aa = get_amino_acid_from_anticodon(mutant_anticodon)

                # 只保留编码不同氨基酸的突变体
                if mutant_aa != original_aa and mutant_aa != 'X':
                    alternatives.append((mutant_anticodon, mutant_aa))
                    logger.debug(f"  单碱基突变: 位置{i+1} {original_anticodon[i]}->{base}: {mutant_anticodon} -> {mutant_aa}")

    logger.debug(f"单碱基突变找到 {len(alternatives)} 个候选")
    return alternatives

def mutate_anticodon_in_sequence(sequence, structure, original_anticodon, new_anticodon):
    """根据二级结构符号精确定位并替换tRNA序列中的反密码子"""
    logger = get_logger('design_mutant_library')
    logger.debug(f"根据结构符号替换反密码子: {original_anticodon} -> {new_anticodon}")

    # structure参数已经是纯净的结构符号字符串
    structure_str = structure.strip()

    if not structure_str:
        logger.warning("结构符号为空")
        return fallback_simple_replace(sequence, original_anticodon, new_anticodon, logger)

    logger.debug(f"结构符号: {structure_str}")

    # 使用结构符号定位反密码子
    start_pos, end_pos = find_anticodon_position_by_structure(structure_str)

    if start_pos is None or end_pos is None:
        logger.warning("无法通过结构符号定位反密码子")
        return fallback_simple_replace(sequence, original_anticodon, new_anticodon, logger)

    # 验证序列中指定位置的反密码子
    original_anticodon_dna = original_anticodon.replace('U', 'T')
    actual_anticodon = sequence[start_pos:end_pos]

    if actual_anticodon != original_anticodon_dna:
        logger.warning(f"结构定位的位置 {start_pos}-{end_pos-1} 反密码子 '{actual_anticodon}' 与预期 '{original_anticodon_dna}' 不匹配")
        return fallback_simple_replace(sequence, original_anticodon, new_anticodon, logger)

    # 精确替换指定位置的反密码子
    new_anticodon_dna = new_anticodon.replace('U', 'T')
    mutated_sequence = sequence[:start_pos] + new_anticodon_dna + sequence[end_pos:]

    logger.info(f"成功在位置 {start_pos}-{end_pos-1} 替换反密码子: {actual_anticodon} -> {new_anticodon_dna}")
    logger.debug(f"原序列片段: ...{sequence[max(0,start_pos-5):end_pos+5]}...")
    logger.debug(f"新序列片段: ...{mutated_sequence[max(0,start_pos-5):end_pos+5]}...")

    return mutated_sequence

def find_anticodon_position_by_structure(structure_str):
    """根据tRNA二级结构符号定位反密码子位置"""

    logger = get_logger('design_mutant_library')

    # 查找 "<.>" 的位置
    stem_loop_index = structure_str.find("<.>")
    logger.debug(f"找到 <.> 中的“<”的位置: {stem_loop_index}")
    # 计算反密码子起始位置：找到 "<.>" 后第一个 "." 的位置，然后加2
    i = stem_loop_index + 2
    while i < len(structure_str) and structure_str[i] == ">":
        i += 1
    anticodon_start = i + 3
    anticodon_end = anticodon_start + 3
    logger.debug(f"找到反密码子起始位置: {anticodon_start}, 反密码子结束位置: {anticodon_end}")
    return anticodon_start, anticodon_end

def fallback_simple_replace(sequence, original_anticodon, new_anticodon, logger):
    """简单替换序列中的反密码子（不使用结构符号）"""
    original_anticodon_dna = original_anticodon.replace('U', 'T')
    new_anticodon_dna = new_anticodon.replace('U', 'T')

    if original_anticodon_dna in sequence:
        mutated_sequence = sequence.replace(original_anticodon_dna, new_anticodon_dna)
        logger.debug(f"成功替换反密码子")
        return mutated_sequence

    logger.warning(f"无法在序列中找到反密码子 {original_anticodon}")
    return sequence

def generate_mutant_library(tRNA_records, out_file):
    """生成tRNA突变体库"""
    logger = get_logger('design_mutant_library')
    logger.info(f"开始为 {len(tRNA_records)} 个tRNA生成突变体...")

    output_records = []

    for tRNA_record in tRNA_records:
        logger.info(f"处理 {tRNA_record.tRNA_id} ({tRNA_record.amino_acid_type})...")

        # 清理反密码子格式
        original_anticodon = tRNA_record.anticodon.upper().replace('T', 'U')
        original_aa = tRNA_record.amino_acid_type

        logger.debug(f"原始反密码子: {original_anticodon} -> {original_aa}")

        # 添加原始tRNA序列
        original_record = SeqRecord(
            Seq(tRNA_record.seq),
            id=f"{tRNA_record.tRNA_id}_original",
            description=f"aa:{original_aa}"
        )
        output_records.append(original_record)

        # 生成突变体
        mutant_anticodons = generate_alternative_anticodons(original_anticodon, original_aa)

        logger.info(f"为 {tRNA_record.tRNA_id} 生成了 {len(mutant_anticodons)} 个突变体")
        if len(mutant_anticodons) == 0:
            logger.warning(f"未能为 {tRNA_record.tRNA_id} 生成任何突变体")
            continue

        for i, (mutant_anticodon, mutant_aa) in enumerate(mutant_anticodons, 1):
            logger.info(f"  突变体{i}: {mutant_anticodon} -> {mutant_aa}")

            # 在序列中替换反密码子
            mutated_sequence = mutate_anticodon_in_sequence(
                tRNA_record.seq,
                tRNA_record.structure,
                original_anticodon,
                mutant_anticodon
            )

            # 创建突变体记录
            mutant_record = SeqRecord(
                Seq(mutated_sequence),
                id=f"{tRNA_record.tRNA_id}_mutant{i}",
                description=f" aa:{mutant_aa}"
            )
            output_records.append(mutant_record)
            logger.debug(f"创建突变体记录: {mutant_record.id}, 序列长度: {len(mutated_sequence)}")

    # 写入FASTA文件
    logger.info(f"将 {len(output_records)} 个序列写入 {out_file}...")
    try:
        SeqIO.write(output_records, out_file, "fasta-2line")
        logger.info(f"突变体库生成完成！输出文件: {out_file}")
        logger.info(f"总计生成: {len([r for r in output_records if 'original' in r.id])} 个原始tRNA + {len([r for r in output_records if 'mutant' in r.id])} 个突变体")
    except Exception as e:
        logger.error(f"写入FASTA文件失败: {e}")
        raise

def main():
    # 初始化日志系统
    logger = setup_logger(__file__)
    logger.info("="*60)
    logger.info("tRNA突变体设计脚本开始运行")
    logger.info("="*60)

    try:
        args = parse_arguments()
        logger.info(f"命令行参数:")
        logger.info(f"  候选tRNA文件: {args.input_IDs_file}")
        logger.info(f"  ID映射文件: {args.input_ID_map}")
        logger.info(f"  结构文件: {args.input_structure}")
        logger.info(f"  输出文件: {args.out_file}")

        # 检查输入文件是否存在
        for file_path, file_desc in [(args.input_IDs_file, "候选tRNA文件"),
                                     (args.input_ID_map, "ID映射文件"),
                                     (args.input_structure, "结构文件")]:
            if not os.path.exists(file_path):
                logger.error(f"{file_desc}不存在: {file_path}")
                raise FileNotFoundError(f"{file_desc}不存在: {file_path}")
            logger.info(f"{file_desc}存在: {file_path}")

        # 获取 tRNA ID 对应的序列和二级结构
        logger.info("\n步骤1: 准备tRNA数据")
        tRNA_records = tRNA_prepare(args.input_IDs_file, args.input_ID_map, args.input_structure)

        if not tRNA_records:
            logger.error("未找到任何有效的tRNA记录")
            return

        logger.info(f"成功准备了 {len(tRNA_records)} 个tRNA记录")

        # 生成突变体
        logger.info("\n步骤2: 生成突变体库")
        generate_mutant_library(tRNA_records, args.out_file)

        logger.info("\n" + "="*60)
        logger.info("tRNA突变体设计脚本运行完成")
        logger.info("="*60)

    except Exception as e:
        logger.error(f"脚本运行出错: {e}")
        logger.exception("详细错误信息:")
        raise


if __name__ == "__main__":
    main()
