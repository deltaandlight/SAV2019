#!/bin/bash
#SBATCH -J SAV-3B           # 任务名
#SBATCH -p cpu                   # 提交到 cpu 分区
#SBATCH -t 18:00:00              # 
#SBATCH -N 1                     # 申请 1 个节点
#SBATCH --ntasks-per-node=4    # 每个节点进程
#SBATCH --cpus-per-task=1        # 每个进程占用 1 个 core

module add mpich/3.2.1           # 添加 mpich/3.2.1 模块，注意，不带 -pmi 后缀
export PETSC_DIR=/home/yangchao/zl/petsc-3.6.4
export PETSC_ARCH=linux-gnu-c-debug
# 运行程序
mpiexec -n 4 ./main -ksp_monitor -ksp_converged_reason -ksp_atol 1.e-13 -ksp_rtol 1.e-6\
  -ksp_type gmres -ksp_gmres_restart 30 -ksp_pc_side right -pc_type asm -pc_asm_type restrict -pc_asm_overlap 2 -sub_ksp_type preonly -sub_pc_type lu\
  -da_grid_x 1024 -da_grid_y 1024 -BDFi 3 -BDFAorB 1 -beta 0 -eps 0.1 -endT 0.032 -Tmin 1e-3 -Tmax 1e-3 -interP 1

