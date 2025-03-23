Codon Optimization Script

简介

本项目旨在对外源基因进行密码子优化，以提高外源蛋白在特定宿主（如细菌或植物）中的表达效率。通过加载宿主密码子使用频率表（codon_usage.txt）和外源基因序列（cas9_ori.txt），脚本会优化基因序列使其更符合宿主的密码子使用偏好，并计算优化前后的 Codon Adaptation Index (CAI) 值，还可以生成相关图表用于结果可视化。

功能

加载密码子使用频率表：从指定文件中解析密码子与氨基酸的映射及其使用频率。
基因序列优化：根据宿主密码子偏好表优化外源基因序列。
CAI计算：对优化前后基因进行 CAI 值的计算。
结果输出：将优化后的基因序列与 CAI 值保存到文件。
结果可视化：生成并保存包括优化前后 CAI 值对比的图片（SVG、PDF 和 PNG 格式）。
目录结构

假设以下所有文件均在 /public/home/zhangxiaolin/lyr/taxus/codon_opt_250323 目录下：

/public/home/zhangxiaolin/lyr/taxus/codon_opt_250323/
│── codon_usage.txt         # 密码子使用频率文件
│── cas9_ori.txt            # 外源基因原始序列文件
│── optimized_seq.txt       # 优化后的基因序列（脚本生成）
│── cai_values.txt          # 优化前后的 CAI 值对比（脚本生成）
│── cai_plot.{svg/pdf/png}  # CAI 可视化图表（脚本生成）
└── script.py               # 此仓库的主要脚本
依赖

此脚本依赖以下 Python 库：

Biopython：用于处理 DNA 序列。
Collections：用于字典结构。
Matplotlib：用于数据可视化。
安装必要的库：

pip install biopython matplotlib
使用方法

1. 准备输入文件

密码子使用频率表 （codon_usage.txt）
文件格式如下，每个氨基酸及其对应密码子的使用频率需正确分组：
Amino Acid:Isoleucine
ATA: 0.15
ATC: 0.45
ATT: 0.40
Amino Acid:Leucine
CTA: 0.10
CTC: 0.20
... （其余氨基酸）
外源基因序列（cas9_ori.txt）
文件内为连续的 DNA 序列：
ATGGACATTGACAAGAAGTTGACTTT...
2. 运行脚本

python script.py
3. 输出文件

优化后的基因序列（optimized_seq.txt）
保存优化后的 DNA 序列，用于后续实验。
CAI 值文件（cai_values.txt）
保存优化前后的密码子 CAI 值对比，每个密码子对应一行。
图片文件
cai_plot.svg：矢量图格式，适合编辑。
cai_plot.pdf：便于打印和分享。
cai_plot.png：便于快速查看。
脚本结构说明

主要分为以下几个模块：

密码子使用表加载（load_codon_usage）
将 codon_usage.txt 文件解析为嵌套的字典结构 {amino_acid: {codon: freq}}。
基因序列优化（optimize_gene）
根据宿主密码子使用频率，选择频率最高的同义密码子替换原始密码子。
CAI计算（calculate_single_cai）
逐密码子比对频率，计算优化前后的适应性 CAI 值。
可视化（plot_cai）
绘制 CAI 值对比图表，便于观察优化效果。
关键输入与输出文件名

文件名	描述
codon_usage.txt	宿主密码子使用频率表
cas9_ori.txt	外源基因的原始DNA序列
optimized_seq.txt	脚本生成：优化后的DNA序列文件
cai_values.txt	脚本生成：每个密码子的 CAI 值对比表
cai_plot.svg/pdf/png	脚本生成：绘制的优化前后 CAI 值对比图表
示例

假设：

codon_usage.txt 内容如下：
Amino Acid:Isoleucine
ATA: 0.15
ATC: 0.45
ATT: 0.40
...
cas9_ori.txt 含有如下序列：
ATGGACATTGACAAGAAGTTGACTTT...
运行脚本后：

optimized_seq.txt 优化后的序列：
ATGACCATTGACCAGAAAATTGACCC...
cai_values.txt CAI 数据：
Codon    Original CAI    Optimized CAI
1          0.40           0.85
2          0.30           0.80
...
生成图表文件 cai_plot.svg/pdf/png，如下所示： 示例图表
注意事项

确保 codon_usage.txt 和 cas9_ori.txt 的路径和内容格式正确。
文件目录需要写权限，脚本生成文件会保存到指定路径。
如果要适配不同物种，请提供对应物种的 codon_usage.txt。
联系方式

如有问题，请联系：liyirou@caas.cn
或者提交 Issue 提交至本项目的 GitHub 页面。
