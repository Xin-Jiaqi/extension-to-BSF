import numpy as np
import pandas as pd
import os
import json
import sqlite3
import ase.db
from ase.io import write

def prehandle_c2db():
    # Step1： connect database and get the features（功能）
    con = sqlite3.connect("c2db-2022-11-30.db")
    df = pd.read_sql_query("select * from systems", con)

    # Step2: convert database to csv
    key_value_pairs = []
    for line in df.key_value_pairs:
        key_value_pairs.append(json.loads(line))
    key_value_pairs = pd.DataFrame(key_value_pairs)

    # Step3: rename formula（fromula 指化学式）
    formula_values = list(key_value_pairs['folder'].str.split('/'))
    flag = []
    for line in formula_values:
        flag.append(line[-2])
    key_value_pairs['formula'] = flag

    # Step4：修改 lgnum 范围条件，并将结果保存到不同文件夹
    for lgnum in range(8, 49):  # 遍历空间群 8 到 48
        structures = key_value_pairs[key_value_pairs['lgnum'] == lgnum]
        if not structures.empty:
            folder_name = f'lgnum_{lgnum}'
            os.makedirs(folder_name, exist_ok=True)  # 创建文件夹
            structures.to_csv(f'{folder_name}/sorted_by_condition.csv', index=False)

def output_structure():
    dbpath = os.path.join(os.getcwd(), 'c2db-2022-11-30.db')
    db = ase.db.connect(dbpath)

    def float2line(mart):
        string = ''
        for line in mart:
            for s in line:
                string += '%21.10f' % float(s)
            string += '\n'
        return string

    def count_element(mart):
        string = ['', '']
        elements = sorted(set(mart), key=mart.index)
        for key in elements:
            string[0] += key + '    '
            string[1] += str(mart.count(key)) + '    '
        string[0] += '\n'
        string[1] += '\n'
        return string

    # Step5：遍历每个 lgnum 文件夹并输出结构
    for lgnum in range(8, 49):  # 遍历所有 lgnum 文件夹
        folder_name = f'lgnum_{lgnum}'
        csv_path = f'{folder_name}/sorted_by_condition.csv'
        if os.path.exists(csv_path):
            csv = pd.read_csv(csv_path)
            for i in range(len(csv)):
                atom = db.get_atoms(uid=csv.iloc[i]['uid'], add_additional_information=True)
                formula = csv.iloc[i]['formula']
                # 输出 CIF 文件
                write(f'{folder_name}/{formula}.cif', atom, format='cif')
                # 输出 VASP 文件
                with open(f'{folder_name}/{formula}.vasp', 'w') as f:
                    f.write(csv.iloc[i]['uid'] + '\n')
                    f.write('1.0\n')
                    f.write(float2line(atom.get_cell()))
                    atomInfo = count_element(atom.get_chemical_symbols())
                    f.write(atomInfo[0])
                    f.write(atomInfo[1])
                    f.write('Cartesian\n')
                    f.write(float2line(atom.get_positions()))

if __name__ == "__main__":
    prehandle_c2db()
    output_structure()