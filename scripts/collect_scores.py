#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
收集 Rosetta RNA-蛋白质复合物结构预测的 score
"""

import os
import argparse
import pandas as pd
import subprocess
import io
from concurrent.futures import ThreadPoolExecutor, as_completed


def read_task_list(input_dir):
    sample_ids = []
    for entry in os.listdir(input_dir):
        entry_path = os.path.join(input_dir, entry)
        if os.path.isdir(entry_path):
            sample_id = entry.split("/")[-1].split(".")[0]
            sample_ids.append(sample_id)
    return sample_ids

def run_score_jd2_for_sample(input_dir, sample_id):
    """为单个样本运行score_jd2命令"""
    sample_dir = input_dir + "/" + sample_id
    score_jd2_cmd = "score_jd2 -in:file:silent default.out -out:file:scorefile scores.sc"

    if os.path.exists(sample_dir + "/scores.sc"):
        print(f"{sample_id} 的 score 已存在，跳过")
        return sample_id, "skipped"

    try:
        subprocess.run(score_jd2_cmd, shell=True, check=True, cwd=sample_dir)
        print(f"{sample_id} 的 score 计算完成")
        return sample_id, "completed"
    except subprocess.CalledProcessError as e:
        print(f"{sample_id} 的 score 计算失败: {e}")
        return sample_id, "failed"

def calculate_score(input_dir, sample_ids, max_workers=75):
    """使用多线程并行计算score"""
    print(f"开始使用 {max_workers} 个线程并行计算 {len(sample_ids)} 个样本的 score...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_sample = {executor.submit(run_score_jd2_for_sample, input_dir, sample_id): sample_id
                           for sample_id in sample_ids}

        # 收集结果
        completed_count = 0
        skipped_count = 0
        failed_count = 0

        for future in as_completed(future_to_sample):
            sample_id = future_to_sample[future]
            try:
                sample_id, status = future.result()
                if status == "completed":
                    completed_count += 1
                elif status == "skipped":
                    skipped_count += 1
                elif status == "failed":
                    failed_count += 1
            except Exception as exc:
                print(f"{sample_id} 产生异常: {exc}")
                failed_count += 1

    print(f"\n任务完成统计:")
    print(f"  完成: {completed_count}")
    print(f"  跳过: {skipped_count}")
    print(f"  失败: {failed_count}")
    print(f"  总计: {len(sample_ids)}")

def sum_scores(input_dir, sample_ids):
    # 汇总 score 到输出文件
    all_scores = []
    for sample_id in sample_ids:
        sample_dir = os.path.join(input_dir, sample_id)
        score_path = os.path.join(sample_dir, "scores.sc")
        if not os.path.exists(score_path):
            print(f"{sample_id} 的 score 不存在，跳过")
            continue
        with open(score_path) as f:
            lines = [line for line in f if line.startswith('SCORE:')]
        if len(lines) < 2:
            print(f"{sample_id} 的 scores.sc 文件内容不足，跳过")
            continue
        try:
            df = pd.read_csv(io.StringIO(''.join(lines)), sep=r'\s+', engine='python')

            # 删除第一列 "SCORE:"
            if "SCORE:" in df.columns:
                df.drop(columns=["SCORE:"], inplace=True)

            # 选 total_score 最小的一行
            best_row = df.loc[df['total_score'].idxmin()].to_frame().T  # 转回 DataFrame
            best_row["sample_id"] = sample_id
            all_scores.append(best_row)
        except Exception as e:
            print(f"{sample_id} 处理失败: {e}")
            continue
    if not all_scores:
        return pd.DataFrame()
    scores = pd.concat(all_scores, ignore_index=True)
    return scores

def main():
    parser = argparse.ArgumentParser(description='收集 Rosetta RNA-蛋白质复合物结构预测的 score')
    parser.add_argument('--input_dir', type=str, required=True, help='输入目录')
    parser.add_argument('--out_file', type=str, required=True, help='输出文件')
    args = parser.parse_args()

    # 读取输入文件夹下面的文件名, 识别 sample_id
    sample_ids = read_task_list(args.input_dir)

    # 利用 score_jd2 计算 score
    calculate_score(args.input_dir, sample_ids)

    # 汇总 score 到输出文件
    scores = sum_scores(args.input_dir, sample_ids)
    scores.to_csv(args.out_file, index=False)

if __name__ == '__main__':
    main()
