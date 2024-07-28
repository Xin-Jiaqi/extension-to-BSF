import os
import numpy as np
import pandas as pd
import shutil


def read_poscar_lattice_vectors(poscar_file):
    """
    从POSCAR文件中读取晶格常数信息，并计算非周期方向及周期方向的信息。

    Parameters:
        poscar_file (str): POSCAR文件路径。

    Returns:
        dict: 包含晶格信息及计算结果的字典。
    """
    with open(poscar_file, 'r') as file:
        lines = file.readlines()

    # 读取晶格常数信息
    lattice_vectors = [list(map(float, line.split())) for line in lines[2:5]]

    # 计算每个晶格向量的模长
    lengths = [np.linalg.norm(vector) for vector in lattice_vectors]

    # 找到模长最大的方向（非周期方向）
    non_periodic_direction = lengths.index(max(lengths))

    # 找到两个周期方向的索引
    periodic_directions = [i for i in range(3) if i != non_periodic_direction]

    # 提取两个周期方向的向量
    vector1 = lattice_vectors[periodic_directions[0]]
    vector2 = lattice_vectors[periodic_directions[1]]

    # 计算两个周期方向的模长
    length1 = lengths[periodic_directions[0]]
    length2 = lengths[periodic_directions[1]]

    # 计算两个周期方向的夹角（角度制）
    dot_product = np.dot(vector1, vector2)
    angle_radians = np.arccos(dot_product / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
    angle_degrees = np.degrees(angle_radians)  # 将弧度制转换为角度制

    return {
        "length1": length1,
        "length2": length2,
        "angle_degrees": angle_degrees
    }


def process_vasp_files(base_directory):
    """
    遍历给定目录中的所有.vasp文件，提取周期方向信息并存储到Excel文件。
    同时，将满足条件的.vasp及.cif文件复制到新建的文件夹中。

    Parameters:
        base_directory (str): 要遍历的目录路径。

    Returns:
        None
    """
    # 初始化数据列表
    data = []

    # 定义centered lattice的lgnum值
    centered_lattice_lgnums = {10, 13, 18, 22, 26, 35, 36, 47, 48}

    # 遍历 base_directory 下的所有子文件夹
    for folder in os.listdir(base_directory):
        if folder.startswith("lgnum_"):  # 确保文件夹名称以 'lgnum_' 开头
            lgnum = int(folder.split('_')[1])  # 提取并转换为整数
            directory = os.path.join(base_directory, folder)

            # 新建一个文件夹用于存放筛选后的.vasp和.cif文件
            filtered_dir = os.path.join(directory, 'filtered_vasp_files')
            os.makedirs(filtered_dir, exist_ok=True)

            # 遍历文件夹中的所有.vasp文件
            for filename in os.listdir(directory):
                if filename.endswith('.vasp'):
                    file_path = os.path.join(directory, filename)
                    print(f"Processing {filename} in {folder}...")

                    # 读取并处理每个.vasp文件
                    result = read_poscar_lattice_vectors(file_path)

                    # 提取化学式
                    chemical_formula = filename.split('.')[0]

                    # 追加到数据列表中
                    data.append(
                        [chemical_formula, result["length1"], result["length2"], result["angle_degrees"], lgnum])

    # 创建DataFrame
    df = pd.DataFrame(data, columns=["Chemical Formula", "Length1", "Length2", "Angle (degrees)", "lgnum"])

    # 确保lgnum列是整数类型
    df["lgnum"] = df["lgnum"].astype(int)

    # 使用 ExcelWriter 将数据写入到多个 Sheet
    output_file = os.path.join(base_directory, 'vasp_results.xlsx')
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # 将所有数据写入 Sheet1
        df.to_excel(writer, sheet_name='Sheet1', index=False)

        # 筛选出满足条件的数据
        filtered_df = df[
            ((df["Length1"] - df["Length2"]).abs() < 0.2) &
            (
                    ((df["Angle (degrees)"] >= 85) & (df["Angle (degrees)"] <= 95)) |
                    ((df["Angle (degrees)"] >= 115) & (df["Angle (degrees)"] <= 125))
            )
            ]

        # 将满足条件的.vasp和.cif文件复制到新文件夹中
        for _, row in filtered_df.iterrows():
            chemical_formula = row["Chemical Formula"]
            lgnum = row["lgnum"]
            folder_name = f'lgnum_{lgnum}'
            directory = os.path.join(base_directory, folder_name)
            filtered_dir = os.path.join(directory, 'filtered_vasp_files')

            # 复制.vasp文件
            source_vasp_file = os.path.join(directory, f"{chemical_formula}.vasp")
            destination_vasp_file = os.path.join(filtered_dir, f"{chemical_formula}.vasp")
            shutil.copy2(source_vasp_file, destination_vasp_file)
            print(f"Copied {chemical_formula}.vasp to {filtered_dir}")

            # 复制.cif文件（如果存在）
            source_cif_file = os.path.join(directory, f"{chemical_formula}.cif")
            if os.path.exists(source_cif_file):
                destination_cif_file = os.path.join(filtered_dir, f"{chemical_formula}.cif")
                shutil.copy2(source_cif_file, destination_cif_file)
                print(f"Copied {chemical_formula}.cif to {filtered_dir}")

        # 按照lgnum值分割数据为simple lattice和centered lattice
        simple_lattice_df = filtered_df[~filtered_df["lgnum"].isin(centered_lattice_lgnums)]
        centered_lattice_df = filtered_df[filtered_df["lgnum"].isin(centered_lattice_lgnums)]

        # 分离simple lattice的不同角度范围
        simple_lattice_85_95 = simple_lattice_df[
            (simple_lattice_df["Angle (degrees)"] >= 85) & (simple_lattice_df["Angle (degrees)"] <= 95)
            ]

        simple_lattice_115_125 = simple_lattice_df[
            (simple_lattice_df["Angle (degrees)"] >= 115) & (simple_lattice_df["Angle (degrees)"] <= 125)
            ]

        # 分离centered lattice的不同角度范围
        centered_lattice_85_95 = centered_lattice_df[
            (centered_lattice_df["Angle (degrees)"] >= 85) & (centered_lattice_df["Angle (degrees)"] <= 95)
            ]

        centered_lattice_115_125 = centered_lattice_df[
            (centered_lattice_df["Angle (degrees)"] >= 115) & (centered_lattice_df["Angle (degrees)"] <= 125)
            ]

        # 将不同范围的数据写入到不同的Sheet
        simple_lattice_85_95.to_excel(writer, sheet_name='simple lattice (85-95)', index=False)
        simple_lattice_115_125.to_excel(writer, sheet_name='simple lattice (115-125)', index=False)

        centered_lattice_85_95.to_excel(writer, sheet_name='centered lattice (85-95)', index=False)
        centered_lattice_115_125.to_excel(writer, sheet_name='centered lattice (115-125)', index=False)

    print(f"Results have been saved to {output_file}")


# 指定文件夹路径
vasp_directory = './'  # 请根据实际情况修改文件夹路径

# 处理所有.vasp文件并输出到Excel
process_vasp_files(vasp_directory)