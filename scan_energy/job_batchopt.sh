#!/bin/bash
#SBATCH -o vasp.out
#SBATCH -e vasp.err
#SBATCH -A hpc1106163703
#SBATCH --partition=C032M0128G
#SBATCH -J vasp
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-node=32

# 定义VASP执行路径
VASP_EXEC="/gpfs/share/home/1301111484/apps/vasp-wannier/vasp_std"

# 循环遍历40个文件夹
for dir in */; do
    cd "$dir"
    
    # 判断是否已经收敛
    if grep -q "reached required accuracy" OUTCAR; then
        echo "目录 $dir: 计算已收敛，跳过作业提交。"
    else
        # 如果没有CONTCAR文件
        if [ ! -f CONTCAR ]; then
            echo "目录 $dir: 未收敛且没有CONTCAR，提交作业。"
            sbatch <<EOF
#!/bin/bash
#SBATCH -o vasp.out
#SBATCH -e vasp.err
#SBATCH -A hpc1106163703
#SBATCH --partition=C032M0128G
#SBATCH -J vasp
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-node=32

mpirun -n 32 $VASP_EXEC
EOF
        else
            # 如果已经有CONTCAR文件，重命名并提交作业
            echo "目录 $dir: 未收敛，但已有CONTCAR，将其复制为POSCAR并重新提交作业。"
            cp CONTCAR POSCAR
            sbatch <<EOF
#!/bin/bash
#SBATCH -o vasp.out
#SBATCH -e vasp.err
#SBATCH -A hpc1106163703
#SBATCH --partition=C032M0128G
#SBATCH -J vasp
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-node=32

mpirun -n 32 $VASP_EXEC
EOF
        fi
    fi

    cd ..
done
