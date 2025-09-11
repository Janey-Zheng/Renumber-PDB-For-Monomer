# -*- coding: utf-8 -*-
from pathlib import Path

ROOT = Path(r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\pdbbind_seq")
PATTERNS = ("*_renumbered_vs_ref.pdb", "*_renumbered_vs_ref.txt")

DRY_RUN = False

deleted = 0
targets = []

for pat in PATTERNS:
    for p in ROOT.rglob(pat):
        # 只要“子文件夹中的”文件；如需包含根目录本身的文件，把下面这一行条件去掉
        if p.is_file() and p.parent != ROOT:
            targets.append(p)

# 预览
print(f"将处理 {len(targets)} 个文件（DRY_RUN={DRY_RUN}）：")
for p in targets:
    print(" -", p)

# 真删
if not DRY_RUN:
    for p in targets:
        try:
            p.unlink()
            deleted += 1
        except Exception as e:
            print(f"[失败] {p} -> {e}")
    print(f"已删除 {deleted} 个文件。")
