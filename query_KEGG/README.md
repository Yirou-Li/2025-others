# KEGG ID 查询脚本

这是一个 Python 脚本，用于从 KEGG 在线数据库中查询基因的 KEGG ID，并将结果保存到 TSV 文件中。脚本支持以下功能：

1. **根据 COG 或 EC_number 查询 KEGG ID**：
   - 如果 `COG` 存在，通过 `COG` 查询 KEGG ID。
   - 如果 `EC_number` 存在，通过 `EC_number` 查询 KEGG ID。

2. **断点续传**：
   - 每次查询后，将结果保存到临时文件中。
   - 如果脚本中断，下次运行时可以从临时文件中读取已查询的结果，继续查询未处理的行。

3. **自动重试机制**：
   - 如果查询失败（例如网络中断），脚本会等待 1 分钟后重试，最多重试 3 次。

4. **并行处理**：
   - 使用 `ThreadPoolExecutor` 并行处理查询，提高检索速度。

5. **命令行支持**：
   - 脚本支持通过命令行参数传入 `.tsv` 文件路径，自动生成输出文件和临时文件。

---

## 使用说明

### 1. 安装依赖

确保已安装以下 Python 库：

```bash
pip install requests pandas beautifulsoup4
```

### 2. 运行脚本

将脚本保存为 `query_kegg_c.py`，然后在终端中运行以下命令：

```bash
python3 query_kegg_c.py /path/to/your_file.tsv
```

### 3. 输入文件格式

输入文件应为 TSV 格式，包含以下列：

| locus_tag       | ftype | length_bp | gene   | EC_number | COG      | product                              |
|-----------------|-------|-----------|--------|-----------|----------|--------------------------------------|
| GGAONODO_00001  | CDS   | 1353      | dnaA_1 |           | COG0593  | Chromosomal replication initiator protein DnaA |
| GGAONODO_00002  | CDS   | 1137      | dnaN   |           | COG0592  | Beta sliding clamp                   |
| GGAONODO_00038  | CDS   | 750       | yfiC   | 2.1.1.223 |          | tRNA1(Val) (adenine(37)-N6)-methyltransferase |

### 4. 输出文件格式

输出文件为 TSV 格式，包含以下列：

| locus_tag       | ftype | length_bp | gene   | EC_number | COG      | product                              | KEGG                  |
|-----------------|-------|-----------|--------|-----------|----------|--------------------------------------|-----------------------|
| GGAONODO_00001  | CDS   | 1353      | dnaA_1 |           | COG0593  | Chromosomal replication initiator protein DnaA | K02313                |
| GGAONODO_00002  | CDS   | 1137      | dnaN   |           | COG0592  | Beta sliding clamp                   | K02311                |
| GGAONODO_00038  | CDS   | 750       | yfiC   | 2.1.1.223 |          | tRNA1(Val) (adenine(37)-N6)-methyltransferase | K02469;K02470;K02621  |

---

## 批量处理脚本

如果你需要批量处理多个 `.tsv` 文件，可以使用以下 Shell 脚本：

### **Shell 脚本：`query_all_kegg.sh`**

```bash
#!/bin/bash

# 设置输入路径
input_base="/home/lyr/Microbial_fermentation/CGMCC/Geobacillus_fasta/prokka_output"
script_path="/home/lyr/Microbial_fermentation/CGMCC/Geobacillus_fasta/prokka_output/KEGG/query_kegg_c.py"

# 遍历每个子文件夹
for dir in "$input_base"/*; do
  if [ -d "$dir" ]; then
    # 在子文件夹中查找 .tsv 文件
    for tsv_file in "$dir"/*.tsv; do
      if [ -f "$tsv_file" ]; then
        # 输出当前处理的文件
        echo "Processing: $tsv_file"
        # 调用 Python 脚本进行 KEGG 查询
        python3 "$script_path" "$tsv_file"
      fi
    done
  fi
done

echo "All KEGG queries have been processed."
```

### **运行步骤**
1. 将 Shell 脚本保存为 `query_all_kegg.sh`。
2. 为脚本添加执行权限：

```bash
chmod +x query_all_kegg.sh
```

3. 运行脚本：

```bash
./query_all_kegg.sh
```

---

## 脚本参数

脚本中的以下参数可以根据需要调整：

1. **并行线程数**：
   ```python
   max_workers = 5  # 设置并行线程数
   ```

2. **重试次数和延迟**：
   ```python
   retries = 3  # 最大重试次数
   delay = 60   # 重试延迟时间（秒）
   ```

---

## 示例

### 输入文件

`input.tsv`：

| locus_tag       | ftype | length_bp | gene   | EC_number | COG      | product                              |
|-----------------|-------|-----------|--------|-----------|----------|--------------------------------------|
| GGAONODO_00001  | CDS   | 1353      | dnaA_1 |           | COG0593  | Chromosomal replication initiator protein DnaA |
| GGAONODO_00002  | CDS   | 1137      | dnaN   |           | COG0592  | Beta sliding clamp                   |
| GGAONODO_00038  | CDS   | 750       | yfiC   | 2.1.1.223 |          | tRNA1(Val) (adenine(37)-N6)-methyltransferase |

### 输出文件

`output.tsv`：

| locus_tag       | ftype | length_bp | gene   | EC_number | COG      | product                              | KEGG                  |
|-----------------|-------|-----------|--------|-----------|----------|--------------------------------------|-----------------------|
| GGAONODO_00001  | CDS   | 1353      | dnaA_1 |           | COG0593  | Chromosomal replication initiator protein DnaA | K02313                |
| GGAONODO_00002  | CDS   | 1137      | dnaN   |           | COG0592  | Beta sliding clamp                   | K02311                |
| GGAONODO_00038  | CDS   | 750       | yfiC   | 2.1.1.223 |          | tRNA1(Val) (adenine(37)-N6)-methyltransferase | K02469;K02470;K02621  |

---

## 注意事项

1. **网络连接**：
   - 确保你的网络可以访问 KEGG 在线数据库。

2. **查询频率**：
   - KEGG 对查询频率有限制，建议根据实际情况调整并行线程数和延迟时间。

3. **临时文件**：
   - 如果脚本中断，临时文件会保存已查询的结果。下次运行时，脚本会从临时文件中读取已查询的结果，继续查询未处理的行。

---

## 许可证

本项目采用 [MIT 许可证](LICENSE)。
