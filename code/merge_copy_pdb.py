import os
import shutil

# 源目录和目标目录
src_dir = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\PDBbind\PDBbind_structure"
dst_dir = r"F:\Monomer\Cluster_Tables\supplement_inf\pocket\sequences\pdbbind_seq"

# 遍历源目录下的子文件夹
for folder in os.listdir(src_dir):
    subfolder_path = os.path.join(src_dir, folder)
    
    if os.path.isdir(subfolder_path):
        # 在目标目录下新建一个同名子文件夹
        dst_subfolder = os.path.join(dst_dir, folder)
        os.makedirs(dst_subfolder, exist_ok=True)

        # 遍历子文件夹中的文件
        for file in os.listdir(subfolder_path):
            if file.endswith("_clean_singlechain.pdb"):
                src_file = os.path.join(subfolder_path, file)
                dst_file = os.path.join(dst_subfolder, file)
                
                # 复制文件
                shutil.copy2(src_file, dst_file)
                print(f"已复制: {src_file} -> {dst_file}")

print("所有文件复制完成！")
