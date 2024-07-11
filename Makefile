interactive:
	salloc --account=bbpa-delta-gpu --partition=gpuA40x4-interactive \
     --nodes=2 --gpus-per-node=1 \
     --tasks-per-node=4 --cpus-per-task=8\

# after allocating an interactive session, run on the same terminal with srun <your_command>
srun:
	srun python generate_data.py

#submit a batch job
sbatch-kb:
	sbatch launchkb.slurm
sbatch-ds:
	sbatch ../launchds.slurm
preprocess:
	python graphwm/preprocess/preprocess.py split_protein ./graphwm/datasets/proteins/ ./graphwm/datasets/protein_train_ready

eval-env:
	srun --mem=64g --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --partition=gpuA100x4-interactive,gpuA40x4-interactive --account=bbpa-delta-gpu --gpus-per-node=1 --time=01:00:00 --pty /bin/bash
eval-protein:
	python eval.py --config-name eval_protein

sbatch-preprocess:
	sbatch graphwm/preprocess.slurm

# TRAIN
# MAKE SURE TO CHECK THE TRAJECTORY LENGTHS (10k or 2k??) what is idx?
setup:
	srun --mem=16g --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --partition=gpuA100x4-interactive,gpuA40x4-interactive --account=bbpa-delta-gpu --gpus-per-node=2 --time=01:00:00 --pty /bin/bash

train:
	python mlcgmd/train.py