import os
import shutil

# 定义文件夹路径
poscar_dir = "poscars"
input_dir = "inputs"
output_base_dir = "."

# 遍历 i 和 j 创建文件夹并复制文件
for i in range(1, 9):
    for j in range(1, 6):
        # 创建"ij"文件夹
        folder_name = f"{i}{j}"
        folder_path = os.path.join(output_base_dir, folder_name)
        os.makedirs(folder_path, exist_ok=True)

        # 找到对应的POSCAR文件
        poscar_file = f"{i}{j}.vasp"
        poscar_src_path = os.path.join(poscar_dir, poscar_file)

        # 复制POSCAR文件到"ij"文件夹
        if os.path.exists(poscar_src_path):
            shutil.copy(poscar_src_path, os.path.join(folder_path, "POSCAR"))
        else:
            print(f"POSCAR文件 {poscar_file} 不存在！")

        # 复制inputs文件夹的内容到"ij"文件夹
        if os.path.exists(input_dir):
            for item in os.listdir(input_dir):
                item_path = os.path.join(input_dir, item)
                shutil.copy(item_path, folder_path)
        else:
            print(f"inputs文件夹不存在！")
