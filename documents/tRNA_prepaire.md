# trna
本文档用于生成分析报告，包含详细的分析过程和结果。

## 环境配置&软件安装

使用pixi管理环境，初始化项目。
```bash
cd /path/to/project/trna
pixi init
```

安装tRNAscan-SE
```bash
git clone https://github.com/UCSC-LoweLab/tRNAscan-SE.git
cd tRNAscan-SE
./configure
make
sudo make install
```

## tRNA预测

1. 使用tRNAscan-SE预测家蚕和草地贪夜蛾的tRNA。
```bash
# 家蚕
cd /home/gyk/reference/bm_ncbi/GCF_030269925.1
tRNAscan-SE -o work/tRNAscan-SE/bmo-tRNAs.out \
-f work/tRNAscan-SE/bmo-tRNAs.ss \
-s work/tRNAscan-SE/bmo-tRNAs.iso \
-m work/tRNAscan-SE/bmo-tRNAs.stats \
-b work/tRNAscan-SE/bmo-tRNAs.bed \
-j work/tRNAscan-SE/bmo-tRNAs.gff \
-a work/tRNAscan-SE/bmo-tRNAs.fa \
-l work/tRNAscan-SE/bmo-tRNAs.log \
-d -D -H --detail --thread 70 /home/gyk/reference/bm_ncbi/GCF_030269925.1/GCF_030269925.1_ASM3026992v2_genomic_chr.fa

# 草地贪夜蛾
cd /home/gyk/reference/sf_ncbi
tRNAscan-SE -o work/tRNAscan-SE/sfr-tRNAs.out \
-f work/tRNAscan-SE/sfr-tRNAs.ss \
-s work/tRNAscan-SE/sfr-tRNAs.iso \
-m work/tRNAscan-SE/sfr-tRNAs.stats \
-b work/tRNAscan-SE/sfr-tRNAs.bed \
-j work/tRNAscan-SE/sfr-tRNAs.gff \
-a work/tRNAscan-SE/sfr-tRNAs.fa \
-l work/tRNAscan-SE/sfr-tRNAs.log \
-d -D -H --detail --thread 70 /home/gyk/reference/sf_ncbi/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic_chr.fa
```

2. 生成tRNA名称映射文件
```bash
python scripts/generate_trna_name_map.py -i work/tRNAscan-SE/bmo-tRNAs.ss -o work/tRNAscan-SE/bmo-tRNAs-name_map.txt
python scripts/generate_trna_name_map.py -i work/tRNAscan-SE/sfr-tRNAs.ss -o work/tRNAscan-SE/sfr-tRNAs-name_map.txt
```
3. 使用 tDRnamer 制作数据库，以获得二级结构信息。

```bash
# 在项目根目录下创建sifs文件夹
mkdir sifs & cd sifs

# 使用 Apptainer 拉取 tDRnamer 镜像
apptainer pull tdrnamer.sif docker://ucsclowelab/tdrnamer

# 使用 Apptainer 运行 tDRnamer
cd ../work/tRNAscan-SE
# 家蚕
apptainer exec /home/gyk/project/trna/sifs/tdrnamer.sif \
create_tDRnamer_db --db bmo --genome /home/gyk/reference/bm_ncbi/GCF_030269925.1/GCF_030269925.1_ASM3026992v2_genomic_chr.fa \
--trna bmo-tRNAs.out \
--ss bmo-tRNAs.ss \
--namemap bom-tRNAs-name_map.txt --force

# 草地贪夜蛾
apptainer exec /home/gyk/project/trna/sifs/tdrnamer.sif \
create_tDRnamer_db --db sfr --genome /home/gyk/reference/sf_ncbi/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic_chr.fa \
--trna sfr-tRNAs.out \
--ss sfr-tRNAs.ss \
--namemap sfr-tRNAs-name_map.txt --force
```

4. 计算正交得分

以草地贪夜蛾为查询物种，家蚕为参考物种，找到能在草地贪夜蛾中被识别而在家蚕中未被识别的tRNA。即tRNA所属氨基酸类型的正交得分小于0的tRNA。

```bash
python scripts/trna_orthogonal_score.py \
-q work/tRNAscan-SE/sfr-tDRnamer_db/sfr-trnaalign.stk \
-t work/tRNAscan-SE/bmo-tDRnamer_db/bmo-trnaalign.stk \
-o results/Sf_in_Bm \
-e data/identity_elements.txt
```

以家蚕为查询物种，草地贪夜蛾为参考物种，找到能在家蚕中被识别而在草地贪夜蛾中未被识别的tRNA。即tRNA所属氨基酸类型的正交得分小于0的tRNA。

```bash
python scripts/trna_orthogonal_score.py \
-q work/tRNAscan-SE/bmo-tDRnamer_db/bmo-trnaalign.stk \
-t work/tRNAscan-SE/sfr-tDRnamer_db/sfr-trnaalign.stk \
-o results/Bm_in_Sf \
-e data/identity_elements.txt
```
