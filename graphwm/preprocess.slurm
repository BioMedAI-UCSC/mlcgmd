#!/bin/bash
#SBATCH --job-name="preProcessData"
#SBATCH --output="output.txt"
#SBATCH --partition=gpuA40x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # 4 because each node has 4 gpus
#SBATCH --cpus-per-task=16  # spread out to use 1 core per numa, set to 64 if tasks is 1
#SBATCH --constraint="scratch"
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=closest   # select a cpu close to gpu on pci bus topology
#SBATCH --account=bbpa-delta-gpu
#SBATCH --exclusive  # dedicated node for this job
#SBATCH --no-requeue
#SBATCH -t 48:00:00
#SBATCH --array=0

export OMP_NUM_THREADS=2  # if code is not multithreaded, otherwise set to 8 or 16
echo "TASK_ID", $SLURM_ARRAY_TASK_ID
# srun -N 1 -n 1 python generate_data.py $SLURM_ARRAY_TASK_ID

source /projects/bbpa/coarseGrained/miniconda3/bin/activate kbcoarse
srun python preprocess/preprocess.py split_protein datasets/proteins datasets/protein_split
# py-torch example, --ntasks-per-node=1 --cpus-per-task=64
# srun python3 multiple_gpu.py

# nodes mean number of nodes NODES NODELIST(REASON) ex  4 gpub[020,030,035,038]
# -N means number of nodes to run on 
# -n means number of tasks to run