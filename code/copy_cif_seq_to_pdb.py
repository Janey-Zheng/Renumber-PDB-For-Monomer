import os
import shutil

# === 路径配置 ===
mapping_file = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\INDEX_general_PL_name.txt"   # 你的映射表文件
uniprot_dir = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\str_seq"
pdbbind_seq_dir = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\pdbbind_seq"

# 读取目标目录下已存在的 PDB 子文件夹（小写集合，做大小写无关匹配）
existing_pdb_dirs = {
    d.lower() for d in os.listdir(pdbbind_seq_dir)
    if os.path.isdir(os.path.join(pdbbind_seq_dir, d))
}

# 解析映射表：第1列是 PDB ID，第3列是 UniProt ID
# 只保留 PDB ID 在 existing_pdb_dirs 中的行
# 一个 PDB 只保留最后一条映射（若同一 PDB 出现多次，可改成列表）
filtered_map = {}  # pdb_id(lower) -> uniprot_id
with open(mapping_file, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()  # 空格/Tab 都能分
        if len(parts) < 3:
            continue
        pdb_id, uniprot_id = parts[0].lower(), parts[2]
        if pdb_id in existing_pdb_dirs:
            filtered_map[pdb_id] = uniprot_id

print(f"将在 {len(filtered_map)} 个已存在的 PDB 子文件夹中复制序列。")

copied, missing_txt = 0, 0
for pdb_id, uniprot_id in filtered_map.items():
    dst_subfolder = os.path.join(pdbbind_seq_dir, pdb_id)
    src_txt = os.path.join(uniprot_dir, f"{uniprot_id}.txt")
    if os.path.exists(src_txt):
        dst_txt = os.path.join(dst_subfolder, f"{uniprot_id}.txt")
        # 如果不想覆盖已存在文件，可改为：if not os.path.exists(dst_txt): shutil.copy2(...)
        shutil.copy2(src_txt, dst_txt)
        copied += 1
        print(f"✔ 复制 {src_txt} -> {dst_txt}")
    else:
        missing_txt += 1
        print(f"⚠ 未找到 UniProt 文本：{src_txt}（PDB: {pdb_id}）")

print(f"完成：复制 {copied} 个；缺失文本 {missing_txt} 个；跳过（因无对应子文件夹）的已在筛选阶段处理。")
