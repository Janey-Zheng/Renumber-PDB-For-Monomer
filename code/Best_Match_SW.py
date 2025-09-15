# -*- coding:utf-8 -*-

import os
import re
import csv
from collections import namedtuple
from typing import List, Tuple

three_to_one = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
    'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'HSD':'H','HSE':'H','HSP':'H','MSE':'M','SEC':'U','PYL':'O'
}

# 定义一个不可变的轻量级类Residue，包含chain/resseq/icode/resname四个字段
Residue = namedtuple("Residue", "chain resseq icode resname")

# 读取参考序列及编号
# 支持两种格式(纯文本包含Sequence:/Auth Seq IDs:1-2-3-...和FASTA:>...)
def load_reference_sequence_and_numbers(ref_txt_path):
    seq, nums = "", []
    with open(ref_txt_path, "r", encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"Sequence:\s*([A-Z]+)", text)
    if m:
        seq = m.group(1).strip()
    else:    # FASTA
        lines = [ln.strip() for ln in text.splitlines() if not ln.strip().startswith(">")]
        letters = "".join(re.findall(r"[A-Z]+", "".join(lines)))
        if letters:
            seq = letters
    # Auth 序列编号
    n = re.search(r"Auth\s+Seq\s+IDs:\s*([0-9,\-\s]+)", text)
    if n:
        raw = n.group(1)
        parts = re.split(r"[,\-\s]+", raw.strip())
        nums = [int(x) for x in parts if x.isdigit()]
    else:
        nums = list(range(1, len(seq) + 1))
    return seq, nums

# 解析结构：补齐缺口，占位'-'，链统一为 'A'(针对单体)
def parse_pdb_residues(str_PDB_path):
    residues_real = []           # 真实存在的残基（用于确定范围）
    seen = set()
    residue_by_num = {}          # {resseq: Residue}
    letter_by_num  = {}          # {resseq: 'A'..'Z' or 'X'}
    with open(str_PDB_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            resname = line[17:20].strip()
            if line.startswith("HETATM") and resname not in ("MSE", "SEC", "PYL"):
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            chain = "A"  # 统一为 A 链
            try:
                resseq = int(line[22:26].strip())
            except ValueError:
                continue
            icode  = line[26].strip()
            resid = (chain, resseq, icode)
            if resid in seen:
                continue
            seen.add(resid)
            r = Residue(chain, resseq, icode, resname)
            residues_real.append(r)
            residue_by_num[resseq] = r
            letter_by_num[resseq]  = three_to_one.get(resname.upper(), "X")
    if not residues_real:
        return [], ""
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
            # 缺号占位
            residues_filled.append(Residue("A", n, "", "-"))
            seq_letters.append("-")
    return residues_filled, "".join(seq_letters)

# 根据给定的参考序列和对应的编号列表，生成一个带gap的序列片段，并返回该片段的范围
def build_ref_letters_and_axis(ref_seq, ref_nums):
    if not ref_seq or not ref_nums:
        return "", 0, -1
    lo, hi = min(ref_nums), max(ref_nums)
    slot = {n: ref_seq[i] for i, n in enumerate(ref_nums)}
    letters = "".join(slot.get(n, "-") for n in range(lo, hi + 1))
    return letters, lo, hi

# 局部配对(Smith–Waterman),用于将结构序列和参考序列对齐，输出匹配位置映射关系
def _sw_local_align_map(struct_letters, ref_letters):
    # 压缩结构序列
    s_compact = []
    s_map = []   # 原索引 -> 压缩后索引（若原位是 '-'，置 -1）
    cnt = 0
    for i, ch in enumerate(struct_letters):
        if ch == '-':
            s_map.append(-1)
        else:
            s_map.append(cnt)
            s_compact.append(ch)
            cnt += 1
    s = "".join(s_compact)
    t = ref_letters
    # struct_letters = "A-CG"
    # s_compact = ["A", "C", "G"] → "ACG"
    # s_map = [0, -1, 1, 2]
    # 定义打分规则
    def score(a, b):
        if a == 'X' or b == 'X':
            return 0
        return 2 if a == b else -1
    gap = -2
    # 创建动态规划矩阵
    m, n = len(s), len(t)
    H = [[0]*(n+1) for _ in range(m+1)]
    P = [[0]*(n+1) for _ in range(m+1)]  # 0 stop, 1 diag, 2 up, 3 left
    best, bi, bj = 0, 0, 0
    # H[i][j]: s[0..i) 与 t[0..j) 的最优局部对齐分数
    # P[i][j]: 记录方向（0 停止，1 对角，2 上，3 左）
    # 进行填表
    for i in range(1, m+1):
        si = s[i-1]
        for j in range(1, n+1):
            tj = t[j-1]
            diag = H[i-1][j-1] + score(si, tj)
            up   = H[i-1][j] + gap
            left = H[i][j-1] + gap
            h = max(0, diag, up, left)
            H[i][j] = h
            if h == 0:
                P[i][j] = 0
            elif h == diag:
                P[i][j] = 1
            elif h == up:
                P[i][j] = 2
            else:
                P[i][j] = 3
            if h > best:
                best, bi, bj = h, i, j
    # best, bi, bj: 记录最大分数和对应位置
    # 回溯,找到最优路径
    pairs = []
    i, j = bi, bj
    while i > 0 and j > 0 and P[i][j] != 0:
        if P[i][j] == 1:      # diag
            pairs.append((i-1, j-1))
            i -= 1; j -= 1
        elif P[i][j] == 2:    # up: 结构插入缺口（跳过一个结构字符）
            i -= 1
        else:                 # left: 参考插入缺口
            j -= 1
    pairs.reverse()
    # 从最佳位置 (bi, bj) 开始往回走
    # diag: 对齐成功 → 记录配对 (i-1, j-1)
    # up: 结构多了一个字符（gap 对齐），跳过结构字符
    # left: 参考多了一个字符，跳过参考字符
    # pairs中存储压缩序列下标对
    # 映射回原始序列
    compact2orig = []
    for orig_idx, comp_idx in enumerate(s_map):
        if comp_idx != -1:
            if comp_idx == len(compact2orig):
                compact2orig.append(orig_idx)
    pairs_orig = [(compact2orig[i_c], j) for (i_c, j) in pairs]
    # 示例:[(0, 0), (2, 1), (3, 2)]
    # 结构的第0个残基对应参考的第0个, 结构的第2个残基对应参考的第1个
    return pairs_orig

# 根据局部比对得到的配对结果，把结构中的残基编号重新映射到参考序列的编号轴上，生成逐位编号映射表
def build_number_mapping_from_pairs(residues, lo_ref, pairs_orig):
    if not residues or not pairs_orig:
        return {}
    # 填配对区锚点
    idx2num = {}
    for i_struct, j_ref_axis in pairs_orig:
        target_num = lo_ref + j_ref_axis
        idx2num[i_struct] = target_num
    # 2) 左/右外推
    i_list = sorted(idx2num.keys())
    iL, iR = i_list[0], i_list[-1]
    # 左端 Δ
    delta_L = idx2num[iL] - residues[iL].resseq
    for i in range(0, iL):
        idx2num[i] = residues[i].resseq + delta_L
    # 右端 Δ
    delta_R = idx2num[iR] - residues[iR].resseq
    for i in range(iR + 1, len(residues)):
        idx2num[i] = residues[i].resseq + delta_R
    # 3) 转成 PDB key -> new_resseq
    mapping = {}
    for i, r in enumerate(residues):
        if r.resname == '-':
            continue
        newnum = idx2num.get(i)
        if newnum is not None:
            key = (r.chain, r.resseq, r.icode)
            mapping[key] = int(newnum)
    return mapping

# PDB重编号
def renumber_pdb_with_mapping(str_PDB_path, str_PDB_outpath, PDB_pocket_path, PDB_pocket_outpath, mapping):
    def _rewrite(in_path, out_path):
        with open(in_path, "r", encoding="utf-8", errors="ignore") as fin, open(out_path, "w", encoding="utf-8") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    altloc = line[16].strip()
                    if altloc not in ("", "A"):
                        continue
                    chain = "A"  # 设定为默认A链
                    try:
                        resseq = int(line[22:26])
                    except ValueError:
                        fout.write(line)
                        continue
                    icode = line[26].strip()
                    key = (chain, resseq, icode)
                    if key in mapping:
                        new_resseq = mapping[key]
                        # 写新链 & 新编号 & 清 iCode
                        line = line[:21] + chain + line[22:]
                        line = line[:22] + f"{new_resseq:>4}" + line[26:]
                        line = line[:26] + " " + line[27:]
                fout.write(line)
    _rewrite(str_PDB_path, str_PDB_outpath)
    _rewrite(PDB_pocket_path, PDB_pocket_outpath)

# 计算重编号前的正确匹配率
def match_rate_before_renumbering(residues: List[Residue], ref_seq: str, ref_nums: List[int]) -> float:
    if not residues or not ref_seq or not ref_nums:
        return 0.0
    num_to_idx = {n: i for i, n in enumerate(ref_nums)}
    matches = compares = 0
    for r in residues:
        if r.resname in ('-',):
            continue
        j = num_to_idx.get(r.resseq)
        if j is None or j < 0 or j >= len(ref_seq):
            continue
        a = three_to_one.get(r.resname.upper(), "X")
        b = ref_seq[j]
        if a != "X" and b != "X":
            compares += 1
            if a == b:
                matches += 1
    return (matches / compares) if compares else 0.0

# 计算重编号后的正确匹配率
def match_rate_via_mapping(residues: List[Residue], ref_seq: str, ref_nums: List[int], mapping) -> float:
    if not residues or not ref_seq or not ref_nums or not mapping:
        return 0.0
    num_to_idx = {n: i for i, n in enumerate(ref_nums)}
    matches = compares = 0
    for r in residues:
        if r.resname in ('-',):
            continue
        jnum = mapping.get((r.chain, r.resseq, r.icode))
        if jnum is None:
            continue
        j = num_to_idx.get(jnum)
        if j is None or j < 0 or j >= len(ref_seq):
            continue
        a = three_to_one.get(r.resname.upper(), "X")
        b = ref_seq[j]
        if a != "X" and b != "X":
            compares += 1
            if a == b:
                matches += 1
    return (matches / compares) if compares else 0.0

# 目录与文件的扫描
def get_inf_in_subdir(subdir):
    pdb_id = os.path.basename(subdir.rstrip("\\/"))
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
    ROOT = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\renum_try"
    SUMMARY_CSV = os.path.join(ROOT, "Renumber_Table.csv")
    
    with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(["PDB", "Uniprot", "Matched_Rate_Raw", "Matched_Rate_New"])
        for entry in os.listdir(ROOT):
            print("entry:", entry)
            subdir = os.path.join(ROOT, entry)
            if not os.path.isdir(subdir):
                continue
            str_PDB_path, PDB_pocket_path, ref_txt_path, uniprot, pdb_id = get_inf_in_subdir(subdir)
            if not (str_PDB_path and PDB_pocket_path and ref_txt_path):
                continue
            # 参考序列以及编号
            ref_seq, ref_nums = load_reference_sequence_and_numbers(ref_txt_path)
            # entry: 5t1m
            # ref_seq: RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
            # ref_nums: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,...]
            # 结构序列以及编号
            residues, struct_letters = parse_pdb_residues(str_PDB_path)
            # print("residues:", residues)
            # print("struct_letters:", struct_letters)
            # residues: [...Residue(chain='A', resseq=277, icode='', resname='ARG'), 
            # Residue(chain='A', resseq=278, icode='', resname='-'), ...]
            # struct_letters: ERPREEFTLCRKLGSGYFGEVFEGLWKDR
            # 重编号前的匹配率
            match_rate_raw = match_rate_before_renumbering(residues, ref_seq, ref_nums)
            # print("match_rate_raw:", match_rate_raw)
            # match_rate_raw: 0.06542056074766354
            # # 参考轴
            ref_letters_axis, lo, hi = build_ref_letters_and_axis(ref_seq, ref_nums)
            # ref_letters_axis: RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
            # lo: 1
            # hi: 107
            pairs_orig = _sw_local_align_map(struct_letters, ref_letters_axis)
            # print("pairs_orig:", pairs_orig)
            # pairs_orig: [(107, 0), (108, 1), (109, 2), (110, 3),...
            # # 由配对生成逐位映射（目标编号 = lo + j_ref_axis）
            mapping = build_number_mapping_from_pairs(residues, lo, pairs_orig)
            # print("mapping:", mapping)
            # mapping: {('A', 1, ''): -106, ('A', 2, ''): -105, ('A', 3, ''): -104, ('A', 4, ''): -103,...
            #  ('A', 276, ''): 276, ('A', 277, ''): 277, ('A', 279, ''): 279, ...
            str_PDB_outpath = os.path.join(subdir, f"{pdb_id}_renumbered.pdb")
            PDB_pocket_outpath = os.path.join(subdir, f"{pdb_id}_pocket_renumbered.pdb")
            # 重新编号PDB文件中的残基
            renumber_pdb_with_mapping(str_PDB_path, str_PDB_outpath, PDB_pocket_path, PDB_pocket_outpath, mapping)
            match_rate_new = match_rate_via_mapping(residues, ref_seq, ref_nums, mapping)
            writer.writerow([pdb_id, uniprot or "", f"{match_rate_raw:.6f}", f"{match_rate_new:.6f}"])
            # print(f"[{pdb_id}] raw={match_rate_raw:.4f}, new={match_rate_new:.4f}")

if __name__ == "__main__":
    main()
