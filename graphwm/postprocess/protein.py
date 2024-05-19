import torch as pt
from mdtraj import load_frame, Trajectory, load_hdf5
import numpy as np
import os

from p_tqdm import p_umap


def load_protein_topology(trajectoryFile):
    """
    load the protein topology from the provided trajectory directory
    """
    traj = load_hdf5(trajectoryFile)
    topology = traj.topology
    
# def load_protein_eval_trajectory(evalTraj):
#     evalTraj = np.tra
    