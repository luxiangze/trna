#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
批量生成tRNA和蛋白质复合物结构预测所需的输入文件

该脚本用于生成Rosetta RNA-蛋白质复合物结构预测所需的序列(fasta.txt)和结构(secstruct.txt)文件。
脚本会读取tRNA结构文件，以及蛋白质FASTA文件，
将它们一一配对，并生成相应的输入文件。

使用方法：
    fasta_file_prepare.py \
    --tRNA_ids results/Sf_in_Bm/candidate_tRNAs.csv \
    --tRNA_ids_map work/tRNAscan-SE/sfr-tRNAs-name_map.txt \
    --tRNA_structure work/tRNAscan-SE/sfr-tRNAs-confidence.ss \
    --aaRSs_fasta work/colabfold/Bm_aaRSs.fasta \
    --output_dir work/rosetta
"""

import argparse
import os
from Bio import SeqIO
from generate_trna_name_map import parse_trnascan_output
from logger_utils import setup_logger, get_logger

# 提取 tRNA ID
def extract_tRNA_id(tRNA_file, tRNA_ids_map):
    logger = get_logger('fasta_file_prepare')
    logger.info(f"开始提取tRNA ID，输入文件: {tRNA_file}")
    logger.info(f"映射文件: {tRNA_ids_map}")
    
    tRNA_ids = []
    tRNA_ids_dict = {}
    
    # 读取候选tRNA ID
    with open(tRNA_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('tRNA'):
                tRNA_id = line.split(',')[0]
                tRNA_ids.append(tRNA_id)
                logger.debug(f"第{line_num}行找到tRNA ID: {tRNA_id}")
    
    logger.info(f"从{tRNA_file}中提取到{len(tRNA_ids)}个tRNA ID")
    logger.info(f"tRNA ID列表: {tRNA_ids}")
    
    # 读取映射文件并进行匹配
    mapped_count = 0
    with open(tRNA_ids_map, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line_num == 1:  # 跳过标题行
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                trnascan_id = parts[0]
                gtrnadb_id = parts[1]
                
                for tRNA_id in tRNA_ids:
                    target_pattern = tRNA_id + "-1"
                    if gtrnadb_id == target_pattern:
                        tRNA_ids_dict[tRNA_id] = trnascan_id
                        mapped_count += 1
                        logger.info(f"成功映射: {tRNA_id} -> {trnascan_id}")
                        break
    
    # 记录映射结果
    logger.info(f"映射完成: {mapped_count}/{len(tRNA_ids)} 个tRNA ID成功映射")
    
    # 记录未映射的tRNA ID
    unmapped_ids = [tid for tid in tRNA_ids if tid not in tRNA_ids_dict]
    if unmapped_ids:
        logger.warning(f"以下{len(unmapped_ids)}个tRNA ID未能成功映射:")
        for uid in unmapped_ids:
            logger.warning(f"  未映射: {uid} (查找模式: {uid}-1)")
    
    logger.info(f"最终映射字典: {tRNA_ids_dict}")
    return tRNA_ids_dict

def extract_tRNA_structure(tRNA_ids_dict, tRNA_structure):
    logger = get_logger('fasta_file_prepare')
    logger.info(f"开始提取tRNA结构，结构文件: {tRNA_structure}")
    logger.info(f"需要查找{len(tRNA_ids_dict)}个tRNA的结构信息")
    
    # 解析tRNAscan-SE输出文件
    logger.info("解析tRNAscan-SE输出文件...")
    tRNAs = parse_trnascan_output(tRNA_structure)
    logger.info(f"从结构文件中解析到{len(tRNAs)}个tRNA条目")
    
    tRNA_structures = {}
    found_count = 0
    
    # 为每个映射的tRNA ID查找对应的结构
    for tRNA_id, tRNA_name in tRNA_ids_dict.items():
        logger.debug(f"查找tRNA {tRNA_id} (对应tRNAscan ID: {tRNA_name})的结构")
        
        structure_found = False
        for tRNA in tRNAs:
            if tRNA['id'] == tRNA_name:
                tRNA_structures[tRNA_id] = (tRNA['sequence'], tRNA['structure'])
                found_count += 1
                structure_found = True
                logger.info(f"找到结构: {tRNA_id} -> 序列长度: {len(tRNA['sequence'])}, 结构长度: {len(tRNA['structure'])}")
                logger.debug(f"  序列: {tRNA['sequence']}")
                logger.debug(f"  结构: {tRNA['structure']}")
                break
        
        if not structure_found:
            logger.warning(f"未找到tRNA {tRNA_id} (tRNAscan ID: {tRNA_name})的结构信息")
    
    logger.info(f"结构提取完成: {found_count}/{len(tRNA_ids_dict)} 个tRNA找到了结构信息")
    
    # 记录未找到结构的tRNA
    missing_structures = [tid for tid in tRNA_ids_dict.keys() if tid not in tRNA_structures]
    if missing_structures:
        logger.warning(f"以下{len(missing_structures)}个tRNA未找到结构信息:")
        for mid in missing_structures:
            logger.warning(f"  缺失结构: {mid} (tRNAscan ID: {tRNA_ids_dict[mid]})")
    
    return tRNA_structures

def extract_aaRSs(aaRSs_fasta):
    logger = get_logger('fasta_file_prepare')
    logger.info(f"开始提取aaRS蛋白质序列，FASTA文件: {aaRSs_fasta}")
    
    aaRSs = {}
    count = 0
    for record in SeqIO.parse(aaRSs_fasta, "fasta"):
        aaRSs[record.id] = record.seq
        count += 1
        logger.info(f"提取aaRS: {record.id}, 序列长度: {len(record.seq)}")
        logger.debug(f"  序列: {str(record.seq)[:50]}{'...' if len(record.seq) > 50 else ''}")
    
    logger.info(f"aaRS提取完成，共提取到{count}个蛋白质序列")
    logger.info(f"aaRS ID列表: {list(aaRSs.keys())}")
    return aaRSs

def pair_tRNA_and_aaRS(tRNA_structures, aaRSs):
    logger = get_logger('fasta_file_prepare')
    logger.info(f"开始配对tRNA和aaRS")
    logger.info(f"可用tRNA数量: {len(tRNA_structures)}")
    logger.info(f"可用aaRS数量: {len(aaRSs)}")
    
    pairs = []
    pair_count = 0
    
    for tRNA_id, tRNA_structure in tRNA_structures.items():
        logger.debug(f"为tRNA {tRNA_id} 配对aaRS...")
        for aaRS_id, aaRS_sequence in aaRSs.items():
            pairs.append((tRNA_id, tRNA_structure, aaRS_id, aaRS_sequence))
            pair_count += 1
            logger.debug(f"  创建配对: {tRNA_id} + {aaRS_id}")
    
    expected_pairs = len(tRNA_structures) * len(aaRSs)
    logger.info(f"配对完成: 创建了{pair_count}个配对 (预期: {expected_pairs})")
    
    if pair_count != expected_pairs:
        logger.warning(f"配对数量不匹配! 实际: {pair_count}, 预期: {expected_pairs}")
    
    return pairs

def generate_input_files(pairs, output_dir):
    logger = get_logger('fasta_file_prepare')
    logger.info(f"开始生成输入文件，输出目录: {output_dir}")
    logger.info(f"需要处理{len(pairs)}个配对")
    
    rnp_ids = []
    generated_count = 0
    
    for i, pair in enumerate(pairs, 1):
        tRNA_id, tRNA_structure, aaRS_id, aaRS_sequence = pair
        
        # 处理序列和结构
        aaRS_structure = '.' * len(aaRS_sequence)
        tRNA_sequence = tRNA_structure[0].lower().replace("t", "u")
        tRNA_structure_processed = tRNA_structure[1].replace("<", ")").replace(">", "(")
        
        # 生成RNP ID
        rnp_id = f"{aaRS_id}_{tRNA_id}"
        rnp_ids.append(rnp_id)
        
        # 组合序列和结构
        rnp_sequence = aaRS_sequence + tRNA_sequence
        rnp_structure = aaRS_structure + tRNA_structure_processed
        
        logger.debug(f"处理配对 {i}/{len(pairs)}: {rnp_id}")
        logger.debug(f"  aaRS长度: {len(aaRS_sequence)}, tRNA长度: {len(tRNA_sequence)}")
        logger.debug(f"  总长度: {len(rnp_sequence)}")
        
        # 生成FASTA文件
        fasta_file = os.path.join(output_dir, f"{rnp_id}.fasta")
        with open(fasta_file, "w") as f:
            f.write(f">{rnp_id}\n{rnp_sequence}\n")
        logger.debug(f"  生成FASTA文件: {fasta_file}")
        
        # 生成结构文件
        struct_file = os.path.join(output_dir, f"{rnp_id}.txt")
        with open(struct_file, "w") as f:
            f.write(f"{rnp_structure}\n{rnp_sequence}\n")
        logger.debug(f"  生成结构文件: {struct_file}")
        
        generated_count += 1
        
        if i % 10 == 0 or i == len(pairs):
            logger.info(f"已处理 {i}/{len(pairs)} 个配对")
    
    # 生成RNP ID列表文件
    rnp_ids_file = os.path.join(output_dir, "rnp_ids.txt")
    with open(rnp_ids_file, "w") as f:
        for rnp_id in rnp_ids:
            f.write(f"{rnp_id}\n")
    
    logger.info(f"文件生成完成!")
    logger.info(f"  生成了{generated_count}个FASTA文件")
    logger.info(f"  生成了{generated_count}个结构文件")
    logger.info(f"  生成了RNP ID列表文件: {rnp_ids_file}")
    logger.info(f"  所有文件已保存到: {output_dir}")


def main():
    # 初始化日志系统
    logger = setup_logger(__file__)
    
    parser = argparse.ArgumentParser(description='批量生成tRNA和蛋白质复合物结构预测所需的输入文件')
    parser.add_argument('--tRNA_ids', required=True, help='tRNA候选文件')
    parser.add_argument('--tRNA_ids_map', required=True, help='tRNA ID 映射文件')
    parser.add_argument('--tRNA_structure', required=True, help='tRNAscan-SE输出的.ss文件')
    parser.add_argument('--aaRSs_fasta', required=True, help='aaRSs FASTA文件')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("开始执行fasta_file_prepare.py脚本")
    logger.info("=" * 80)
    logger.info(f"输入参数:")
    logger.info(f"  tRNA候选文件: {args.tRNA_ids}")
    logger.info(f"  tRNA ID映射文件: {args.tRNA_ids_map}")
    logger.info(f"  tRNA结构文件: {args.tRNA_structure}")
    logger.info(f"  aaRS FASTA文件: {args.aaRSs_fasta}")
    logger.info(f"  输出目录: {args.output_dir}")
    
    try:
        # 确保输出目录存在
        logger.info(f"创建输出目录: {args.output_dir}")
        os.makedirs(args.output_dir, exist_ok=True)
        
        # 步骤1: 提取 tRNA ID
        logger.info("\n" + "="*50)
        logger.info("步骤1: 提取tRNA ID")
        logger.info("="*50)
        tRNA_ids_dict = extract_tRNA_id(args.tRNA_ids, args.tRNA_ids_map)
        logger.info(f"步骤1完成: 成功映射{len(tRNA_ids_dict)}个tRNA ID")
        
        # 步骤2: 提取 tRNA 结构
        logger.info("\n" + "="*50)
        logger.info("步骤2: 提取tRNA结构")
        logger.info("="*50)
        tRNA_structures = extract_tRNA_structure(tRNA_ids_dict, args.tRNA_structure)
        logger.info(f"步骤2完成: 获得{len(tRNA_structures)}个tRNA的结构信息")
        
        # 步骤3: 提取蛋白质序列
        logger.info("\n" + "="*50)
        logger.info("步骤3: 提取aaRS蛋白质序列")
        logger.info("="*50)
        aaRSs = extract_aaRSs(args.aaRSs_fasta)
        logger.info(f"步骤3完成: 获得{len(aaRSs)}个aaRS蛋白质序列")
        
        # 步骤4: 配对 tRNA 和蛋白质
        logger.info("\n" + "="*50)
        logger.info("步骤4: 配对tRNA和aaRS")
        logger.info("="*50)
        pairs = pair_tRNA_and_aaRS(tRNA_structures, aaRSs)
        logger.info(f"步骤4完成: 创建了{len(pairs)}个tRNA-aaRS配对")
        
        # 步骤5: 生成输入文件
        logger.info("\n" + "="*50)
        logger.info("步骤5: 生成输入文件")
        logger.info("="*50)
        generate_input_files(pairs, args.output_dir)
        logger.info("步骤5完成: 所有输入文件已生成")
        
        # 执行总结
        logger.info("\n" + "="*80)
        logger.info("执行总结")
        logger.info("="*80)
        logger.info(f"输入tRNA候选数量: 从{args.tRNA_ids}读取")
        logger.info(f"成功映射的tRNA ID: {len(tRNA_ids_dict)}个")
        logger.info(f"获得结构信息的tRNA: {len(tRNA_structures)}个")
        logger.info(f"可用的aaRS蛋白质: {len(aaRSs)}个")
        logger.info(f"最终生成的配对: {len(pairs)}个")
        logger.info(f"输出文件保存位置: {args.output_dir}")
        
        if len(tRNA_structures) < len(tRNA_ids_dict):
            missing_structures = len(tRNA_ids_dict) - len(tRNA_structures)
            logger.warning(f"注意: 有{missing_structures}个tRNA ID映射成功但未找到结构信息")
        
        logger.info("脚本执行完成!")
        
    except Exception as e:
        logger.error(f"脚本执行过程中发生错误: {str(e)}")
        logger.error(f"错误类型: {type(e).__name__}")
        import traceback
        logger.error(f"详细错误信息:\n{traceback.format_exc()}")
        raise

if __name__ == '__main__':
    main()
