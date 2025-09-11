# -*- coding:utf-8 -*-

import os
import re
import csv
import re
from collections import namedtuple
from typing import List, Tuple, Dict
from collections import namedtuple, deque

three_to_one = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
    'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'HSD':'H','HSE':'H','HSP':'H','MSE':'M','SEC':'U','PYL':'O'
}
# 定义一个不可变的轻量级类Residue，包含chain/resseq/icode/resname四个字段
Residue = namedtuple("Residue", "chain resseq icode resname")

# 读取参考序列
def load_reference_sequence_and_numbers(ref_txt_path):
    with open(ref_txt_path, "r") as f:
        text = f.read()
        # 提取参考序列
        m = re.search(r"Sequence:\s*([A-Z]+)", text)
        if m:
            seq = m.group(1).strip()
        else: # 处理FASTA格式的参考序列文件
            lines = [ln.strip() for ln in text.splitlines() if not ln.strip().startswith(">")]
            letters = "".join(re.findall(r"[A-Z]+", "".join(lines)))
            if letters:
                seq = letters
        # 提取参考序列对应的编号
        n = re.search(r"Auth\s+Seq\s+IDs:\s*([0-9,\-\s]+)", text)
        if n:
            raw = n.group(1)
            parts = re.split(r"[,\-\s]+", raw.strip())
            nums = [int(x) for x in parts if x.isdigit()]
        return seq, nums
    
# 解析PDB结构中的残基序列
def parse_pdb_residues(str_PDB_path):
    residues = [] # 保存残基信息
    seen = set() # 去重相同残基
    seq_letters = [] # 存储一字母AA序列
    
    with open(str_PDB_path, "r") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            # 取残基三字母名
            resname = line[17:20].strip()
            # 跳过小分子(MSE/SEC/PYL除外)
            if line.startswith("HETATM") and resname not in ("MSE", "SEC", "PYL"):
                continue
            # 仅取A构象
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            chain = line[21].strip() or " " # 链ID
            resseq = int(line[22:26]) # 残基序号
            icode = line[26].strip() # 插入码
            resid = (chain, resseq, icode) # 唯一标识残基
            if resid in seen:
                continue
            seen.add(resid)
            residues.append(Residue(chain, resseq, icode, resname))
            seq_letters.append(three_to_one.get(resname.upper(), "X"))
    return residues, "".join(seq_letters)

# 计算未改编号前的匹配率
def match_rate_before_renumbering(residues, ref_seq, ref_nums, renumber_start=None):
    pos_map = {}
    for idx, num in enumerate(ref_nums):
        pos_map.setdefault(num, deque()).append(idx)
        # print("pos_map:", pos_map)
        # pos_map: {4: deque([0]), 5: deque([1]), 6: deque([2]), 7: deque([3]), ...
    matches = 0
    for k, r in enumerate(residues):
        # r: Residue(chain='A', resseq=243, icode='', resname='ASP')
        # k: 239
        # 找到结构残基需要对齐到的参考编号
        num = (renumber_start + k) if renumber_start is not None else r.resseq
        dp = pos_map.get(num)
        if not dp: # 如果不存在该编号
            continue
        ref_idx = dp.popleft()
        a = three_to_one.get(r.resname.upper(), "X") # 结构序列残基
        b = ref_seq[ref_idx] # 参考序列残基
        if a != "X" and b != "X" and a == b:
                matches += 1
    match_rate_raw = matches / len(ref_seq)
    return match_rate_raw
                

# 在ref_seq上滑动pdb_seq
def best_sliding_alignment(pdb_seq, ref_seq):
    best = (-1, -1, -1)
    Lp, Lr = len(pdb_seq), len(ref_seq)
    for ref_start in range(-Lp+1, Lr):
        matches = 0
        overlap_len = 0
        for i in range(Lp):
            j = ref_start + i
            if 0 <= j < Lr:
                overlap_len += 1
                a = pdb_seq[i]
                b = ref_seq[j]
                if a != 'X' and b != 'X' and a == b:
                    matches += 1
        if overlap_len > 0:
            if best[1] < matches or (best[1] == matches and best[2] < overlap_len):
                best = (ref_start, matches, overlap_len)
    matched_rate_new = best[1] / Lp
    print(matched_rate_new)
    return best, matched_rate_new

# 计算重编号起点(new_start_num)
def cal_new_start_num(best, ref_nums):
    ref_start, matches, overlap_len = best
    if ref_nums and ref_start is not None and ref_start >= 0 and ref_start < len(ref_nums):
        new_start_num = ref_nums[ref_start]
    elif ref_start is not None and ref_start >= 0:
        new_start_num = ref_start + 1  # 无编号时，默认参考序列从 1 开始
    else:
        if ref_nums:
            offset = -ref_start if ref_start is not None else 0
            new_start_num = max(1, ref_nums[0] - offset)
        else:
            new_start_num = 1
    return new_start_num

# 实现PDB重编号和Pocket重编号
def renumber_pdb(str_PDB_path, str_PDB_outpath, new_start_num, PDB_pocket_path, PDB_pocket_outpath):
    residues, pdb_seq = parse_pdb_residues(str_PDB_path)
    mapping = {}
    current = new_start_num
    for r in residues:
        key = (r.chain, r.resseq, r.icode)
        if key not in mapping:
            mapping[key] = current
            current += 1
    with open(str_PDB_path, "r") as fin, open(str_PDB_outpath, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                altloc = line[16].strip()
                if altloc not in ("", "A"):  # 跳过除主要构象之外的原子
                    continue
                chain = line[21].strip() or " "
                resseq = int(line[22:26])
                icode = line[26].strip()
                key = (chain, resseq, icode)
                if key in mapping:
                    new_resseq = mapping[key]
                    line = line[:22] + f"{new_resseq:>4}" + line[26:] # 把列22-25替换成新编号
                    line = line[:26] + " " + line[27:] # 清空iCode
            fout.write(line)
    with open(PDB_pocket_path, "r") as fin, open(PDB_pocket_outpath, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                altloc = line[16].strip()
                if altloc not in ("", "A"):  # 跳过除主要构象之外的原子
                    continue
                chain = line[21].strip() or " "
                resseq = int(line[22:26])
                icode = line[26].strip()
                key = (chain, resseq, icode)
                if key in mapping:
                    new_resseq = mapping[key]
                    line = line[:22] + f"{new_resseq:>4}" + line[26:] # 把列22-25替换成新编号
                    line = line[:26] + " " + line[27:] # 清空iCode
            fout.write(line)    
    return mapping


# 扫描子文件夹
def get_inf_in_subdir(subdir):
    pdb_id = os.path.basename(subdir.rstrip("\\/")) # 从路径中提取最后一级目录名
    str_PDB_path, PDB_pocket_path, ref_txt_path, uniprot = None, None, None, None
    for name in os.listdir(subdir):
        lower = name.lower()
        path = os.path.join(subdir, name)
        if os.path.isfile(path):
            if str_PDB_path is None and lower.endswith("_protein_clean_singlechain.pdb"):
                str_PDB_path = path
            elif PDB_pocket_path is None and lower.endswith("_pocket_clean_singlechain.pdb"):
                PDB_pocket_path = path
            elif lower.endswith(".txt") and ref_txt_path is None:
                ref_txt_path = path
                uniprot = os.path.splitext(name)[0]
            
    return str_PDB_path, PDB_pocket_path, ref_txt_path, uniprot, pdb_id


def main():
    # PATH
    ROOT = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\pdbbind_seq"
    SUMMARY_CSV = os.path.join(ROOT, "Renumber_Table.csv")
    # 记录残基正确匹配率
    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PDB", "Uniprot", "Matched_Rate_Raw", "Matched_Rate_New"])
        for entry in os.listdir(ROOT):
            print("PDB:", entry)
            subdir = os.path.join(ROOT, entry)
            if not os.path.isdir(subdir):   # 只处理子文件夹
               continue
            str_PDB_path, PDB_pocket_path, ref_txt_path, uniprot, pdb_id = get_inf_in_subdir(subdir)
            ref_seq, ref_nums = load_reference_sequence_and_numbers(ref_txt_path)
            # print("ref_nums:", ref_nums)
            # ref_nums: [4, 5, 6, 7, 8, 9, 10, ...]
            residues, pdb_seq = parse_pdb_residues(str_PDB_path)
            # print("residues:", residues)
            best, matched_rate_new = best_sliding_alignment(pdb_seq, ref_seq)
            ref_start, matches, overlap = best
            new_start_num = cal_new_start_num(best, ref_nums)
            str_PDB_outpath = os.path.join(subdir, f"{pdb_id}_renumbered.pdb")
            PDB_pocket_outpath = os.path.join(subdir, f"{pdb_id}_pocket_renumbered.pdb")
            renumber_pdb(str_PDB_path, str_PDB_outpath, new_start_num, PDB_pocket_path, PDB_pocket_outpath)
            match_rate_raw = match_rate_before_renumbering(residues, ref_seq, ref_nums, renumber_start=None)
            writer.writerow([pdb_id, uniprot, f"{match_rate_raw:.6f}", f"{matched_rate_new:.6f}"])
            
            
            
            
            

if __name__ == "__main__":
    main()
            
            
    