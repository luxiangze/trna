#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
绘制 Rosetta RNA-蛋白质复合物结构预测的 score 散点图

用法：
    plot_scores.py --scores_file scores.csv --out_file scatter.pdf
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description='绘制 Rosetta RNA-蛋白质复合物结构预测的 score 散点图')
    parser.add_argument('--scores_file', type=str, required=True, help='score 文件')
    parser.add_argument('--out_file', type=str, required=True, help='输出文件')
    parser.add_argument('--x_col', type=str, default='total_score', help='X轴列名')
    parser.add_argument('--y_col', type=str, default='fa_rep', help='Y轴列名')
    parser.add_argument('--hue_col', type=str, default='N_WC', help='颜色映射列名')
    parser.add_argument('--size_col', type=str, default='hbond_sc', help='大小映射列名')
    parser.add_argument('--log_scale', action='store_true', help='使用对数坐标轴')
    args = parser.parse_args()

    # 检查文件是否存在
    if not os.path.exists(args.scores_file):
        print(f"错误: 文件 {args.scores_file} 不存在")
        return

    # 创建输出目录（如果不存在）
    os.makedirs(os.path.dirname(args.out_file), exist_ok=True)

    # 读取数据
    try:
        scores = pd.read_csv(args.scores_file)
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return

    # 检查列是否存在
    required_cols = [args.x_col, args.y_col, args.hue_col, args.size_col]
    missing_cols = [col for col in required_cols if col not in scores.columns]
    if missing_cols:
        print(f"错误: 缺少以下列: {', '.join(missing_cols)}")
        print(f"可用列: {', '.join(scores.columns)}")
        return

    # 设置绘图风格
    sns.set_theme(style="whitegrid", font_scale=1.2)
    plt.figure(figsize=(10, 8))

    # 创建散点图
    cmap = sns.cubehelix_palette(reverse=True, as_cmap=True)
    g = sns.relplot(
        data=scores,
        x=args.x_col, y=args.y_col,
        hue=args.size_col, size=args.hue_col,
        palette=cmap, sizes=(20, 200),legend = "auto",
        alpha=1, height=7, aspect=1.2
    )

    # 设置坐标轴
    if args.log_scale:
        g.set(xscale="log", yscale="log")
    g.ax.xaxis.grid(True, "minor", linewidth=.25)
    g.ax.yaxis.grid(True, "minor", linewidth=.25)

    # 添加标题和标签
    # g.fig.suptitle(f"RNA-蛋白质复合物结构预测 Score 分布", fontsize=16)
    g.set_xlabels(args.x_col, fontsize=14)
    g.set_ylabels(args.y_col, fontsize=14)

    # 调整图例
    # g._legend.set_title("参数")

    # 保存图像
    try:
        g.savefig(args.out_file, dpi=300, bbox_inches='tight')
        print(f"图像已保存至 {args.out_file}")
    except Exception as e:
        print(f"保存图像时出错: {e}")

    plt.close()


if __name__ == "__main__":
    main()