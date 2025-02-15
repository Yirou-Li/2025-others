# 2025-others
#Codon Optimization Tool

#This repository contains a Python script for codon optimization of foreign genes based on the codon usage of a given species. The tool calculates the codon usage frequency from a CDS (coding sequence) FASTA file, optimizes a foreign gene sequence, and generates visualizations and tables for analysis.

#Features

#Codon Usage Calculation: Computes the codon usage frequency from a CDS FASTA file.
#Codon Optimization: Optimizes a foreign gene sequence by replacing codons with the most frequently used synonymous codons.
#Visualization: Generates bar charts for codon usage frequency and line plots for CAI (Codon Adaptation Index) values.
#Output Files: Saves the codon usage table, optimized sequence, and CAI values in text files.
#Requirements

#Python 3.7+
#Biopython (pip install biopython)
#Matplotlib (pip install matplotlib)

#Usage

Usage

Clone the Repository:
git clone https://github.com/your-username/codon-optimization.git
cd codon-optimization
Run the Script: Modify the cds_fasta and foreign_gene variables in the script to point to your input files. Then run:
python codon_optimization.py
Output Files:
codon_usage.txt: Codon usage frequency table.
optimized_seq.txt: Optimized foreign gene sequence.
cai_values.txt: CAI values for the original and optimized sequences.
codon_usage_*.png: Bar charts for codon usage frequency of each amino acid.
cai_plot.png: Line plot comparing CAI values of the original and optimized sequences.
Functions

calculate_codon_usage(cds_fasta)

Description: Calculates the codon usage frequency from a CDS FASTA file.
Input: Path to the CDS FASTA file.
Output: Dictionary of codon usage frequencies.
categorize_codon_usage(codon_usage)

Description: Groups codon usage frequencies by amino acid.
Input: Codon usage frequency dictionary.
Output: Dictionary of codon usage frequencies categorized by amino acid.
optimize_gene(gene_sequence, categorized_usage)

Description: Optimizes a foreign gene sequence based on codon usage.
Input: Foreign gene sequence and categorized codon usage dictionary.
Output: Optimized gene sequence.
plot_cai(original_cai_values, optimized_cai_values)

Description: Plots the CAI values of the original and optimized sequences.
Input: Lists of CAI values for the original and optimized sequences.
Output: Matplotlib figure object.
Example

Input

cds_fasta: Path to the CDS FASTA file.
foreign_gene: Foreign gene sequence to optimize (e.g., "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG").
Output

codon_usage.txt:
Amino Acid: I
    ATA: 0.1700
    ATC: 0.4700
    ATT: 0.3600
Amino Acid: M
    ATG: 1.0000
...
optimized_seq.txt: Optimized gene sequence.
cai_values.txt:
Codon   Original CAI   Optimized CAI
1       0.025789       0.035524
2       0.022347       0.028766
...
codon_usage_*.png: Bar charts for each amino acid.
cai_plot.png: Line plot comparing CAI values.
License

This project is licensed under the MIT License. See the LICENSE file for details.

Author

Your Name

Acknowledgments

Biopython: https://biopython.org/
Matplotlib: https://matplotlib.org/


---

### How to Use the README

1. Replace `your-username` with your GitHub username.
2. Update the `Author` section with your name and GitHub profile link.
3. Add any additional acknowledgments or references if needed.
4. Save the file as `README.md` in the root of your repository.

This README provides a clear and concise overview of your project, making it easy for others to understand and use your code. Let me know if you need further adjustments! 😊
