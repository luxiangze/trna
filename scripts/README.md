# 此文件夹包含用于tRNA分析的脚本。

## generate_trna_name_map.py

用于生成tRNA名称映射文件。

用法：
```bash
python generate_trna_name_map.py -i <tRNAscan-SE输出的.ss文件> -o <输出文件>
```

## trna_orthogonal_score.py

用于计算tRNA的正交得分。

用法：
```bash
python trna_orthogonal_score.py -q <查询物种的Stockholm文件> -t <目标物种的Stockholm文件> -o <输出文件夹> -e <identity_elements.txt文件>
```
