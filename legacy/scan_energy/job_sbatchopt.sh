#!/bin/bash
#SBATCH -p <partition>
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 60
#SBATCH -o vasp.out
#SBATCH -e vasp.err
#SBATCH -J vasp

# 脱敏模板：由本地环境显式提供可执行文件，仓库不记录集群路径。
: "${VASP_EXEC:?Set VASP_EXEC to the local VASP executable before submission}"

# 定义SBATCH脚本内容
SBATCH_SCRIPT="#!/bin/bash
#SBATCH -p <partition>
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 60
#SBATCH -o vasp.out
#SBATCH -e vasp.err
#SBATCH -J vasp

mpirun -n 32 \"$VASP_EXEC\""

# 循环遍历所有子目录
for dir in */; do
    cd "$dir"
    
    # 判断是否已经收敛
    if grep -q "reached required accuracy" OUTCAR; then
        echo "目录 $dir: 计算已收敛，跳过作业提交。"
    else
        # 如果没有CONTCAR文件
        if [ ! -f CONTCAR ]; then
            echo "目录 $dir: 未收敛且没有CONTCAR，提交作业。"
        else
            # 如果已经有CONTCAR文件，重命名并提交作业
            echo "目录 $dir: 未收敛，但已有CONTCAR，将其复制为POSCAR并重新提交作业。"
            cp CONTCAR POSCAR
        fi
        sbatch <<< "$SBATCH_SCRIPT"
    fi

    cd ..
done
