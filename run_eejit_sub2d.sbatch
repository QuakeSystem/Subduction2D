#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --gres=gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH -o job.%N.%j.out  # STDOUT
#SBATCH -e job.%N.%j.err  # STDERR
#SBATCH --job-name sub2d_eejit
#SBATCH -t 30-00:00:00


# Eejit:
# There is 4GB of memory per core. Only if this is not submittingcient
# ntasks-per-node should be reduced.

# -- How to run other things --
# - execute script, once *c and *t3c are ok -
# source modload_eejit.sh
# make 
# sbatch run_lab_eejit.sbatch
#
#
# - or get an interactive job: OLD ETH -
# salloc -p fichtner_compute -t 10:00:00 -N 1 --mem=30GB
# than to run on allocated compute node (i.s.o. login node) use:
# srun -n 1 gdb ./i2stm       (-T 2)
# Note if you get comments about time limit than you are not in the
# sallocated terminal


module purge
module load development
module load CUDA/12.5 # on gpu nodes
# module load dev2024 # on cpu nodes


export JULIA_SSL_NO_VERIFY_HOSTS="*"
export OPENSSL_DIR=""
export JULIA_PKG_PRECOMPILE_AUTO=0

echo "Job Configuration:"
echo "==================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "CPUs per Task: $SLURM_CPUS_PER_TASK"
echo "Threads per Core: $SLURM_THREADS_PER_CORE"
echo "Tasks per Node: $SLURM_NTASKS_PER_NODE"
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "==================="
echo ""

echo "Running on hosts: $SLURM_NODELIST"
echo "Current working directory is `pwd`"

echo "  submitting the first job: `date` "
julia --project ./Subduction2D_SZU2019/Subduction2D_DYREL.jl
echo "  first job finished: `date` "
