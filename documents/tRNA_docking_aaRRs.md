# tRNA docking aaRs

分析思路如下：

1. 构建 aaRs 的三维结构
2. 构建 tRNA 和 aaRs 复合物的三维结构
3. 计算复合物的亲合能

## 构建 aaRs 的三维结构

使用 ColabFold 预测 aaRs 的三维结构，详情请见：[ColabFold](https://github.com/luxiangze/alphafold).

把预测结果链接到 `work/colabfold` 目录下。

```bash
ln -s /home/gyk/alphafold/work work/colabfold
```

## 预测 tRNA 和 aaRs 复合物的三维结构

使用 Rosetta 预测 tRNA 和 aaRs 复合物的三维结构，详情请见：[rnp-modeling](https://docs.rosettacommons.org/docs/latest/application_documentation/rna/rnp-modeling)

### 环境配置&软件安装

<!-- ```bash
cd /path/to/project/trna
# # 使用 apptainer 拉取 rosetta 镜像
# apptainer pull sifs/rosetta.sif docker://rosettacommons/rosetta
``` -->

1. 下载 rosetta 预编译的二进制文件

```bash
cd & cd software
wget https://downloads.rosettacommons.org/downloads/academic/3.14/rosetta_bin_ubuntu_3.14_bundle.tar.bz2
extract rosetta_bin_ubuntu_3.14_bundle.tar.bz2
```

2. 安装RNA tools
按照 [RNA-tools](https://docs.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools) 中的说明编辑 .bashrc 文件。

```bash
vim ~/.bashrc # 编辑 RNA_TOOLS 路径
source ~/.bashrc
python $RNA_TOOLS/sym_link.py
```

3. 安装ViennaRNA

```bash
cd & cd software
wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_22_10/viennarna_2.7.0-1_amd64.deb
sudo dpkg -i viennarna_2.7.0-1_amd64.deb

# 如果遇到依赖性问题，可以使用apt-get install -f
```

### run demo

1. 在运行之前，需要配置好 rna_denovo 和其他额外的环境变量。

```bash
# 进入 rosetta 的 bin 目录
cd & cd software/rosetta.binary.ubuntu.release-371/main/source/bin
# 创建 rna_denovo 的软链接
ln -s software/rosetta.binary.ubuntu.release-371/main/source/bin/rna_denovo.static.linuxgccrelease /home/gyk/.local/bin/rna_denovo
# 如果遇到错误：assert( exists( MINI_EXE ) ) AssertionError 则运行下面的命令
ln -s extract_pdbs.static.linuxgccrelease extract_pdbs.linuxgccrelease
```

2. 运行 demo
```bash
# 复制 demo 到当前目录
cp /home/gyk/software/rosetta.binary.ubuntu.release-371/main/demos/public/rnp_structure_prediction demos/ -r

# 进入 demo 目录
cd demos/rnp_structure_prediction

# 运行 demo
rna_denovo @flags

# 提取 PDB 文件
extract_lowscore_decoys.py 2qux_fold_and_dock.out 5

# 如果遇到错误：assert( exists( MINI_EXE ) ) AssertionError 则在 rosetta 的 bin 目录下运行下面的命令，然后再尝试提取 PDB 文件
ln -s extract_pdbs.static.linuxgccrelease extract_pdbs.linuxgccrelease
```

### 开始预测 tRNA 和蛋白质复合物的三维结构

1. 蛋白质 PDB 文件

见构建 aaRs 的三维结构部分内容，PDB 文件和 aaRSs 序列文件位于 `work/colabfold` 目录下。

2. 包含 tRNA 和 蛋白质序列的 FASTA 文件

```bash
# Sf_in_Bm
python scripts/fasta_file_prepare.py \
--tRNA_ids results/Sf_in_Bm/candidate_tRNAs.csv \
--tRNA_ids_map work/tRNAscan-SE/sfr-tRNAs-name_map.txt \
--tRNA_structure work/tRNAscan-SE/sfr-tRNAs-confidence.ss \
--aaRSs_fasta work/colabfold/Bm_aaRSs.fasta \
--output_dir work/rosetta/Sf_in_Bm

# Bm_in_Sf
python scripts/fasta_file_prepare.py \
--tRNA_ids results/Bm_in_Sf/candidate_tRNAs.csv \
--tRNA_ids_map work/tRNAscan-SE/bom-tRNAs-name_map.txt \
--tRNA_structure work/tRNAscan-SE/bmo-tRNAs-confidence.ss \
--aaRSs_fasta work/colabfold/Sf_aaRSs.fasta \
--output_dir work/rosetta/Bm_in_Sf
```

3. 预测 tRNA 和蛋白质复合物的三维结构

```bash
# 使用单个测试数据预测 tRNA 和蛋白质复合物的三维结构
cd demos/rnp_structure_prediction_test
rna_denovo -fasta A0A1Q1NKM3_tRNA-Asn-ATT-1.fasta \
-secstruct_file A0A1Q1NKM3_tRNA-Asn-ATT-1.txt \
-s A0A1Q1NKM3.pdb -minimize_rna false -nstruct 10
# extract_lowscore_decoys.py default.out 5
# fields.py default.out
score_jd2 -in:file:silent default.out 2qux_fold_and_dock.out -out:file:scorefile scores.sc

# 计算接口亲和能
# extract_lowscore_decoys.py default.out 1
# InterfaceAnalyzer \
#  -s default.out.1.pdb \
#  -interface A_B \
#  -pack_separated \
#  -out:file:scorefile default.out.1.sc

# 构建批量预测的任务列表
bash scripts/slurm.sh # 设置工作目录为/home/gyk/project/trna/work/rosetta/Sf_in_Bm
bash scripts/slurm.sh # 设置工作目录为/home/gyk/project/trna/work/rosetta/Bm_in_Sf

# 在计算集群上批量预测
sbatch scripts/slurm.sh # Sf_in_Bm
sbatch scripts/slurm.sh # Bm_in_Sf
```

4. 整理所有预测结果的 score

```bash
# Sf_in_Bm
python scripts/collect_scores.py \
--input_dir work/rosetta/Sf_in_Bm/results \
--out_file work/rosetta/Sf_in_Bm/results/scores.csv

# Bm_in_Sf
python scripts/collect_scores.py \
--input_dir work/rosetta/Bm_in_Sf/results \
--out_file work/rosetta/Bm_in_Sf/results/scores.csv
```

5. 统计 tRNA 和所有蛋白质的平均亲和能

```bash
# Sf_in_Bm
python scripts/candidate_tRNAs_filter.py \
--scores_file work/rosetta/Sf_in_Bm/results/scores.csv \
--out_file results/Sf_in_Bm/tRNAs_score.csv \
--block_list work/rosetta/Sf_in_Bm/block_list.txt

# Bm_in_Sf
python scripts/candidate_tRNAs_filter.py \
--scores_file work/rosetta/Bm_in_Sf/results/scores.csv \
--out_file results/Bm_in_Sf/tRNAs_score.csv \
--block_list work/rosetta/Bm_in_Sf/block_list.txt
```

6. 分别计算家蚕和草地贪夜蛾中 tRNA 和对应的 aaR 的平均亲和能，以作参考。

```bash
# Bm_in_Bm
bash scripts/rna_denovo.sh /home/gyk/project/trna/work/rosetta/Bm_in_Bm /home/gyk/project/trna/work/colabfold/Bm_aaRSs_output_pdb
python scripts/collect_scores.py \
--input_dir /home/gyk/project/trna/work/rosetta/Bm_in_Bm/results \
--out_file /home/gyk/project/trna/work/rosetta/Bm_in_Bm/results/scores.csv
python scripts/candidate_tRNAs_filter.py \
--scores_file work/rosetta/Bm_in_Bm/results/scores.csv \
--out_file results/Bm_in_Bm/tRNAs_score.csv

# Sf_in_Sf
bash scripts/rna_denovo.sh /home/gyk/project/trna/work/rosetta/Sf_in_Sf /home/gyk/project/trna/work/colabfold/Sf_aaRSs_output_pdb
python scripts/collect_scores.py \
--input_dir /home/gyk/project/trna/work/rosetta/Sf_in_Sf/results \
--out_file /home/gyk/project/trna/work/rosetta/Sf_in_Sf/results/scores.csv
python scripts/candidate_tRNAs_filter.py \
--scores_file work/rosetta/Sf_in_Sf/results/scores.csv \
--out_file results/Sf_in_Sf/tRNAs_score.csv
```

### 设计候选 tRNA 突变体库

```bash
# Sf_in_Bm
python scripts/design_mutant_library.py \
--input_IDs_file results/Sf_in_Bm/candidate_tRNAs_id.txt \
--input_ID_map work/tRNAscan-SE/sfr-tRNAs-name_map.txt \
--input_structure work/tRNAscan-SE/sfr-tRNAs-confidence.ss \
--out_file results/Sf_in_Bm/Sf_in_Bm-candidate_tRNAs_mutant_library.fasta

# Bm_in_Sf
python scripts/design_mutant_library.py \
--input_IDs_file results/Bm_in_Sf/candidate_tRNAs_id.txt \
--input_ID_map work/tRNAscan-SE/bom-tRNAs-name_map.txt \
--input_structure work/tRNAscan-SE/bmo-tRNAs-confidence.ss \
--out_file results/Bm_in_Sf/Bm_in_Sf-candidate_tRNAs_mutant_library.fasta
```
