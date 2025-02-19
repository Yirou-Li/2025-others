import requests
import pandas as pd
from bs4 import BeautifulSoup
import time
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

# 定义函数：通过 COG 或 EC_number 查询 KEGG ID
def get_kegg_ids(query, query_type, retries=3, delay=60):
    if pd.isna(query) or query == "":
        return []
    if query_type == "COG":
        # KEGG 的 COG 查询需要将 "COG" 前缀去掉
        cog_id = query.replace("COG", "")
        url = f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=kegg&keywords={cog_id}"
    elif query_type == "EC":
        url = f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=kegg&keywords={query}&page=1"
    else:
        return []
    
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=10)  # 设置超时时间
            soup = BeautifulSoup(response.text, "html.parser")
            kegg_ids = []
            # 查找包含 KEGG ID 的标签
            for a in soup.find_all("a", href=True):
                href = a["href"]
                if "/entry/" in href:
                    kegg_id = href.split("/entry/")[1]
                    if kegg_id.startswith("K"):
                        kegg_ids.append(kegg_id)
            return kegg_ids
        except (requests.RequestException, ConnectionError) as e:
            print(f"Attempt {attempt + 1} failed for {query_type}={query}: {e}")
            if attempt < retries - 1:
                time.sleep(delay)  # 等待一段时间后重试
            else:
                print(f"Max retries reached for {query_type}={query}. Skipping.")
                return []
    return []

# 主函数
def main(input_file):
    # 自动生成输出文件和临时文件路径
    base_dir = os.path.dirname(input_file)
    file_name = os.path.basename(input_file).replace(".tsv", "")
    output_file = os.path.join(base_dir, f"{file_name}_ko.tsv")
    temp_file = os.path.join(base_dir, f"{file_name}_ko_temp.tsv")

    # 读取输入文件
    df = pd.read_csv(input_file, sep="\t")

    # 检查临时文件是否存在
    if os.path.exists(temp_file):
        # 如果临时文件存在，读取已查询的结果
        temp_df = pd.read_csv(temp_file, sep="\t")
        # 合并已查询的结果到原始 DataFrame
        df["KEGG"] = temp_df["KEGG"]
    else:
        # 如果临时文件不存在，初始化 KEGG 列
        df["KEGG"] = ""

    # 定义函数：处理单行查询
    def process_row(index, row):
        # 如果 KEGG 列已有值，跳过该行
        if pd.notna(row["KEGG"]) and row["KEGG"] != "":
            return index, row["KEGG"]
        # 查询 COG 对应的 KEGG ID
        cog_kegg_ids = get_kegg_ids(row["COG"], "COG")
        # 查询 EC_number 对应的 KEGG ID
        ec_kegg_ids = get_kegg_ids(row["EC_number"], "EC")
        # 合并 KEGG ID
        kegg_ids = cog_kegg_ids + ec_kegg_ids
        return index, ";".join(kegg_ids)

    # 使用并行处理查询剩余行
    with ThreadPoolExecutor(max_workers=5) as executor:  # 设置并行线程数
        futures = {executor.submit(process_row, index, row): index for index, row in df.iterrows()}
        for future in as_completed(futures):
            index, kegg_ids = future.result()
            df.at[index, "KEGG"] = kegg_ids
            # 保存中间结果到临时文件
            df.to_csv(temp_file, sep="\t", index=False)
            print(f"Processed row {index + 1}/{len(df)}")

    # 保存最终结果到输出文件
    df.to_csv(output_file, sep="\t", index=False)
    # 删除临时文件
    if os.path.exists(temp_file):
        os.remove(temp_file)

    print(f"KEGG IDs have been added to {output_file}")

# 脚本入口
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python query_kegg_c.py <input_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    main(input_file)
