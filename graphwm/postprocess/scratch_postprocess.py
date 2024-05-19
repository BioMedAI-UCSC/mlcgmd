import torch as pt
from mdtraj import load_frame, Trajectory, load_hdf5
import numpy as np
# import pandas as pd
# from Bio.PDB.PDBParser import PDBParser
# import utils

# eval_protein_dir = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/splits/protein/eval.txt"
eval_protein_h5 = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/proteins/2JLE/result/output_2JLE.h5"

pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/seed_10.pt', map_location=pt.device('cpu'))

original_topology = load_frame(eval_protein_h5, 1)
# original_topology = load_hdf5(eval_protein_h5)

print("Shape of original file xyz 2JLE: ", original_topology.topology)
evaluated_traj = pickle["position"]
# print("Shape of evaluated positions 2JLE100: ", evaluated_traj.shape)

transposed_traj = np.transpose(evaluated_traj, (1, 0, 2))
transposed_traj = transposed_traj.detach().cpu().numpy()
processed_traj = np.divide(transposed_traj, 10)

print("Shape of transposed evaluated positions 2JLE100: ", processed_traj.shape)


# traj = Trajectory(processed_traj, original_topology.topology)
# print(traj)

# traj.save("test.pdb")