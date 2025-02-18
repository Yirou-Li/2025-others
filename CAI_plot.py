import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import matplotlib.pyplot as plt

# 定义工作目录
work_dir = "/public/home/zhangxiaolin/lyr/taxus/codonoptimization"
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

# 函数：计算密码子使用表
def calculate_codon_usage(cds_fasta):
    codon_count = defaultdict(int)  # 存储密码子计数
    total_codons = 0                # 总密码子数

    # 遍历CDS序列
    for record in SeqIO.parse(cds_fasta, "fasta"):
        sequence = str(record.seq).upper()
        # 按密码子分割序列
        for i in range(0, len(sequence) - 2, 3):  # 确保完整的密码子
            codon = sequence[i:i + 3]
            if len(codon) == 3:  # 确保是完整的密码子
                codon_count[codon] += 1
                total_codons += 1

    # 计算密码子使用频率
    codon_usage = {codon: count / total_codons for codon, count in codon_count.items()}
    
    # 验证总频率
    total_freq = sum(codon_usage.values())
    if not abs(total_freq - 1.0) < 1e-6: 
        raise ValueError(f"密码子使用频率总和不为 1，当前总和为: {total_freq}")
    
    return codon_usage

# 函数：将密码子归类到氨基酸
def categorize_codon_usage(codon_usage):
    categorized_usage = defaultdict(dict)
    
    # 遍历密码子使用表
    for codon, freq in codon_usage.items():
        if codon in codon_to_aa:  # 只处理已知密码子
            amino_acid = codon_to_aa[codon]
            categorized_usage[amino_acid][codon] = freq
    
    # 计算每个氨基酸的密码子使用频率之和
    for amino_acid, codons in categorized_usage.items():
        total = sum(codons.values())
        if total > 0:
            for codon in codons:
                codons[codon] /= total  # 归一化，确保总和为 1
    
    return categorized_usage

# 函数：保存密码子使用表到文件
def save_codon_usage(categorized_usage, output_file):
    with open(output_file, "w") as f:
        for amino_acid, codons in categorized_usage.items():
            f.write(f"Amino Acid: {amino_acid}\n")
            for codon, freq in codons.items():
                f.write(f"    {codon}: {freq:.4f}\n")
            f.write("\n")

# 函数：绘制密码子使用频率图表
def plot_codon_usage(categorized_usage, output_dir):
    for amino_acid, codons in categorized_usage.items():
        # 创建图表
        plt.figure(figsize=(10, 6))
        plt.bar(codons.keys(), codons.values(), color='skyblue')
        plt.title(f"Codon Usage for Amino Acid: {amino_acid}")
        plt.xlabel("Codon")
        plt.ylabel("Frequency")
        plt.xticks(rotation=45)
        plt.tight_layout()

        # 保存图表
        plot_file = os.path.join(output_dir, f"codon_usage_{amino_acid}.png")
        plt.savefig(plot_file, dpi=300)
        plt.close()

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

    #返回图形对象
    return fig
# 主程序
if __name__ == "__main__":
    # 输入文件和外源基因
    cds_fasta = "/public/home/zhangxiaolin/lyr/taxus/taxus_T2T/Taxus_chinensis_var_mairei_hap1_T2T_20230807.cds"
    foreign_gene = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"  # 替换为你的外源基因序列

    # 计算密码子使用表
    codon_usage = calculate_codon_usage(cds_fasta)

    # 将密码子使用情况按氨基酸分类
    categorized_usage = categorize_codon_usage(codon_usage)

    # 保存密码子使用表到文件
    codon_usage_file = os.path.join(work_dir, "codon_usage.txt")
    save_codon_usage(categorized_usage, codon_usage_file)
    print(f"密码子使用表已保存到: {codon_usage_file}")

    # 绘制密码子使用频率图表
    plot_codon_usage(categorized_usage, work_dir)
    print(f"密码子使用图表已保存到: {work_dir}")

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
