#!/bin/bash

# 创建一个空的矩阵
declare -A energy_matrix

# 循环读取文件夹中的能量信息
for i in {1..8}; do
    for j in {1..5}; do
        folder="${i}${j}"
        if [ -d "$folder" ]; then
            # 从OUTCAR中提取TOTEN的最后一行
            energy=$(grep TOTEN "$folder/OUTCAR" | tail -1 | awk '{print $5}')
            energy_matrix[$j,$i]=$energy
        else
            echo "文件夹 $folder 不存在。"
        fi
    done
done

# 输出矩阵，j行i列
for j in {1..5}; do
    for i in {1..8}; do
        printf "%s\t" "${energy_matrix[$j,$i]}"
    done
    echo
done
