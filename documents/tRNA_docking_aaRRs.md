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

### 准备输入文件

1.

