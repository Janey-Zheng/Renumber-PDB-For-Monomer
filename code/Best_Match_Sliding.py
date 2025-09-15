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
    residues_real = []            # 真实存在的残基
    seen = set()                  # 去重
    residue_by_num = {}           # {resseq: Residue}
    letter_by_num  = {}           # {resseq: one-letter or 'X'}
    chains = set()

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
            resseq = int(line[22:26].strip()) # 残基序号
            icode = line[26].strip() # 插入码
            resid = (chain, resseq, icode) # 唯一标识残基
            if resid in seen:
                continue
            seen.add(resid)
            chains.add(chain)
            r = Residue(chain, resseq, icode, resname)
            residues_real.append(r)
            residue_by_num[resseq] = r
            letter_by_num[resseq]  = three_to_one.get(resname.upper(), "X")

    if not residues_real:
        return [], ""

    # 用出现频率最高的链作为占位的链（单链文件就是那条链）
    from collections import Counter
    default_chain = Counter(chains).most_common(1)[0][0] if chains else " "

    # 连续编号轴，缺口填 '-'
    min_num = min(r.resseq for r in residues_real)
    max_num = max(r.resseq for r in residues_real)

    residues_filled = []
    seq_letters = []
    for n in range(min_num, max_num + 1):
        if n in residue_by_num:
            r = residue_by_num[n]
            residues_filled.append(r)
            seq_letters.append(letter_by_num[n])
        else:
            # 缺号占位：resname='-', icode=''（空），链用 default_chain
            residues_filled.append(Residue(default_chain, n, "", "-"))
            seq_letters.append("-")

    return residues_filled, "".join(seq_letters)

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
    # print("matches:", matches)
    # print("len_ref_seq:", len(ref_seq))
    match_rate_raw = matches / len(residues)
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
    return best

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
def renumber_pdb(str_PDB_path, str_PDB_outpath, ref_nums, ref_start, PDB_pocket_path, PDB_pocket_outpath):
    residues, pdb_seq = parse_pdb_residues(str_PDB_path)
    # print("residues:", residues)
    # residues: [Residue(chain='A', resseq=185, icode='', resname='GLU'),...add()
    # print("ref_nums:", ref_nums)
    # ref_nums: [10, 11, 12, 13, 14, 15, 16, 17, ...
    Lp = len(residues)
    Lr = len(ref_nums)
    # 计算结构序列与参考序列的overlap窗口, k∈[k0, k1)
    k0 = max(0, -ref_start)
    k1 = min(Lp, Lr - ref_start)
    # 计算整体的偏移量
    if k1 <= k0:
        delta_all = (ref_nums[0] if Lr else 1) - residues[0].resseq
        delta_left = delta_right = delta_all
    else:
        delta_left  = ref_nums[ref_start + k0] - residues[k0].resseq
        delta_right = ref_nums[ref_start + (k1-1)] - residues[k1-1].resseq
    mapping = {}
    assigned = set() # 避免出现重复编号
    for k, r in enumerate(residues):
        # print("k:", k)
        # print("r:", r)
        # k: 0
        # r: Residue(chain='A', resseq=185, icode='', resname='GLU')
        key = (r.chain, r.resseq, r.icode)
        if k0 <= k < k1:
            new_resseq_num = ref_nums[ref_start + k] # 索引查找编号
            # print("new_resseq:", new_resseq)
            # new_resseq: 189 ...
        elif k < k0:
            new_resseq_num = r.resseq + delta_left
        elif k >= k1:
            new_resseq_num = r.resseq + delta_right
        assigned.add(new_resseq_num)
        mapping[key] = new_resseq_num
    # print("assigned:", assigned)
    # print("mapping:", mapping)
    # 写入结构PDB文件
    with open(str_PDB_path, "r") as fin, open(str_PDB_outpath, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                altloc = line[16].strip()
                if altloc not in ("", "A"):  # 跳过除主要构象之外的原子
                    continue
                # chain = line[21].strip() or " "
                chain = "A"
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

# 计算修改过编号后的正确匹配率
def match_rate_by_number(residues, ref_seq, ref_nums):
    num_to_idx = {num: i for i, num in enumerate(ref_nums)}
    matches = 0
    valid_residues = 0  # 分母：仅统计 resseq > 0 的残基数量
    for r in residues:
        if r.resseq <= 0:
            continue
        idx = num_to_idx.get(r.resseq)
        if idx is None or idx >= len(ref_seq):
            continue
        a = three_to_one.get(r.resname.upper(), "X")
        b = ref_seq[idx]
        if a != "X" and b != "X" and a == b:
            matches += 1
        valid_residues += 1
    return matches / valid_residues if valid_residues else 0.0


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
            best = best_sliding_alignment(pdb_seq, ref_seq)
            ref_start, matches, overlap = best
            # new_start_num = cal_new_start_num(best, ref_nums)
            str_PDB_outpath = os.path.join(subdir, f"{pdb_id}_renumbered.pdb")
            PDB_pocket_outpath = os.path.join(subdir, f"{pdb_id}_pocket_renumbered.pdb")
            renumber_pdb(str_PDB_path, str_PDB_outpath, ref_nums, ref_start, PDB_pocket_path, PDB_pocket_outpath)
            match_rate_raw = match_rate_before_renumbering(residues, ref_seq, ref_nums, renumber_start=None)
            residues_new, pdb_seq_new = parse_pdb_residues(str_PDB_outpath)
            match_rate_new = match_rate_by_number(residues_new, ref_seq, ref_nums)
            writer.writerow([pdb_id, uniprot, f"{match_rate_raw:.6f}", f"{match_rate_new:.6f}"])
            

if __name__ == "__main__":
    main()
            
            
    