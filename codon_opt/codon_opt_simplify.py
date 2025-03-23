import os
from Bio.Seq import Seq
from collections import defaultdict
import matplotlib.pyplot as plt

# 定义工作目录
work_dir = "/public/home/zhangxiaolin/lyr/taxus/codon_opt_250323"
os.makedirs(work_dir, exist_ok=True)  # 如果目录不存在则创建

# 密码子到氨基酸的映射字典
codon_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

# 函数：加载密码子使用表
def load_codon_usage(codon_usage_file):
    categorized_usage = defaultdict(dict)
    with open(codon_usage_file, "r") as f:
        lines = f.readlines()
        amino_acid = None
        for line in lines:
            if line.startswith("Amino Acid:"):
                amino_acid = line.split(":")[1].strip()
            elif line.strip() and amino_acid:
                codon, freq = line.strip().split(":")
                categorized_usage[amino_acid][codon.strip()] = float(freq.strip())
    return categorized_usage

# 函数：计算单个密码子的CAI值
def calculate_single_cai(codon, categorized_usage):
    if codon in codon_to_aa:
        amino_acid = codon_to_aa[codon]
        return categorized_usage[amino_acid].get(codon, 0)
    return 0

# 函数：外源基因密码子优化
def optimize_gene(gene_sequence, categorized_usage):
    optimized_sequence = ""
    # 按密码子分割外源基因序列
    for i in range(0, len(gene_sequence), 3):
        codon = gene_sequence[i:i + 3]
        if len(codon) == 3:
            # 将密码子转换为 Seq 对象
            codon_seq = Seq(codon)
            # 找到使用频率最高的同义密码子
            synonyms = [c for c in categorized_usage[codon_seq.translate()].keys()]
            best_codon = max(synonyms, key=lambda x: categorized_usage[codon_seq.translate()].get(x, 0))
            optimized_sequence += best_codon
    return optimized_sequence

# 绘图函数
def plot_cai(original_cai_values, optimized_cai_values):
    # 创建图形对象
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 绘制优化前后的CAI值
    ax.plot(original_cai_values, label='Original Sequence CAI', marker='o', color='skyblue')
    ax.plot(optimized_cai_values, label='Optimized Sequence CAI', marker='x', color='lightpink')
    
    # 计算平均CAI值
    avg_original_cai = sum(original_cai_values) / len(original_cai_values)
    avg_optimized_cai = sum(optimized_cai_values) / len(optimized_cai_values)

    # 添加平均CAI值到图表
    ax.text(0.5, 0.95, f'Original Sequence CAI: {avg_original_cai:.2f}\nOptimized Sequence CAI: {avg_optimized_cai:.2f}',
        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
        fontsize=14, bbox=dict(facecolor='white', alpha=0.5))

    # 设置标题和标签
    ax.set_title('CAI Values for Original and Optimized Sequences', fontsize=14)
    ax.set_xlabel('Codon Number', fontsize=12)
    ax.set_ylabel('Relative Adaptiveness (CAI)', fontsize=12)
    
    # 添加图例和网格
    ax.legend(fontsize=10)
    ax.grid(True)
    
    # 返回图形对象
    return fig
# 主程序
if __name__ == "__main__":
    # 输入文件和外源基因
    codon_usage_file = "/public/home/zhangxiaolin/lyr/taxus/codonoptimization/codon_usage.txt"
    foreign_gene_file = "/public/home/zhangxiaolin/lyr/taxus/codon_opt_250323/cas9_ori.txt"  # 替换为你的外源基因序列
    with open(foreign_gene_file, "r") as f:
        foreign_gene = f.read().strip()

    # 加载密码子使用表
    categorized_usage = load_codon_usage(codon_usage_file)

    # 对外源基因进行密码子优化
    optimized_gene = optimize_gene(foreign_gene, categorized_usage)

    # 保存优化后的序列
    optimized_seq_file = os.path.join(work_dir, "optimized_seq.txt")
    with open(optimized_seq_file, "w") as f:
        f.write(optimized_gene)
    print(f"优化后的序列已保存到: {optimized_seq_file}")

    # CAI计算
    original_cai_values = []
    optimized_cai_values = []
    
    # 计算优化前CAI值
    for i in range(0, len(foreign_gene), 3):
        codon = foreign_gene[i:i + 3]
        if len(codon) == 3:
            ca_value = calculate_single_cai(codon, categorized_usage)
            original_cai_values.append(ca_value)

    # 计算优化后CAI值
    for i in range(0, len(optimized_gene), 3):
        codon = optimized_gene[i:i + 3]
        if len(codon) == 3:
            ca_value = calculate_single_cai(codon, categorized_usage)
            optimized_cai_values.append(ca_value)

    # 保存优化前后的CAI值
    cai_values_file = os.path.join(work_dir, "cai_values.txt")
    with open(cai_values_file, "w") as f:
        f.write("Codon\tOriginal CAI\tOptimized CAI\n")
        for i in range(len(original_cai_values)):
            f.write(f"{i + 1}\t{original_cai_values[i]}\t{optimized_cai_values[i]}\n")
    print(f"优化前后的CAI值已保存到: {cai_values_file}")

    # 生成图表并保存在变量 p1 中
    p1 = plot_cai(original_cai_values, optimized_cai_values)

    # 显示图表
    plt.show()

    # 保存图表为矢量图（支持AI编辑）
    plot_file_svg = os.path.join(work_dir, "cai_plot.svg")
    plot_file_pdf = os.path.join(work_dir, "cai_plot.pdf")
    plot_file_png = os.path.join(work_dir, "cai_plot.png")
    p1.savefig(plot_file_svg, format="svg", bbox_inches="tight")
    p1.savefig(plot_file_pdf, format="pdf", bbox_inches="tight")
    p1.savefig(plot_file_png, format="png", dpi=300, bbox_inches="tight")
    print(f"图表已保存为: {plot_file_svg}, {plot_file_pdf}, {plot_file_png}")
