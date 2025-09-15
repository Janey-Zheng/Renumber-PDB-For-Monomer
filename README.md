# Renumber-PDB-For-Monomer
Given the reference sequence, renumber the PDB structure file and output the matching rate before and after. The renumbering of the pocket can be optionally performed or not.

Use sliding alignment or local alignment (SW) for residue matching. Reference sequences can be structure files or FASTA sequences, preferably without gaps (i.e., missing residues).

# Note
Before alignment, small molecules, ions, etc., in PDB structure files must be removed beforehand (refer to the provided script remove_HOH_HETATM.py). This code is only applicable for modifying single-chain residue numbering. Structures downloaded from databases like PDB require the above processing before subsequent amino acid numbering corrections can be performed.

# Gap
Due to residue deletions or redundant gaps present in some protein structures, the algorithm includes specialized handling for these issues. Sliding alignment may encounter problems when renumbering such structures, so the local alignment algorithm is strongly recommended for renumbering in these cases.

# iCode
For structures containing insertion codes (icodes), first remove the icode residues during structure reading. Then perform local pairing on the sequence without icodes. Finally, renumber the removed icode residues according to optimal renumbering rules (e.g., 62A becomes 101A, 99B becomes 42B).

---

使用滑动配对或者局部配对（SW）进行残基的匹配。参考序列可以是结构文件或者FASTA序列，最好不包含gap，即空缺的残基。

# 注意
在进行配对之前，pdb结构文件中的小分子、离子等需要提前删除（可参考给定的脚本remove_HOH_HETATM.py），并且本代码仅适用于单链残基编号的修改。从PDB等数据库中下载得到的结构需要经过上述处理才能够进行后续的氨基酸编号修正。

# 关于Gap
由于部分蛋白质结构中存在残基缺失或者存在多余gap的情况，因此在算法中特别添加了针对这类问题的处理。滑动配对对存在该类情况的结构进行重编号时可能会出现问题，因此更推荐局部配对算法进行重编号。

# 关于插入码
针对插入码（icode）存在的结构，首先在读入结构的时候删除掉icode残基，然后对无插入码的序列进行局部配对，最后根据最优编号规则对删除的icode残基进行重编号（例如62A改为101A， 99B改为42B）。
