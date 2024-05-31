import torch as pt
from mdtraj import load_frame, Trajectory, load_hdf5
import numpy as np
import os
# import mdTraj

def protein_postprocess(data_dir, eval_list, model_eval):
    """
    data_dir: string
        directory of the original dataset
    eval_list: string
        path to the file with a list of the evaluated protein names
    model_eval: string
        path to the evaluated protein file
    return: protein pdb file(s)
    """
    # gets list of protein names/splits that were evaluated
    f = open(eval_list, "r")
    evaluated_proteins = [protein for protein in f]
    f.close()
    
    # load pickle file
    pk = pt.load(model_eval)
    evaluated_trajectories = pk["position"]
    
    # process and split evaluated trajectories based on the split size of the original protein
    
    
    
    # def process_one_file()
    
    # load evaluated protein

    
    def getProteinTopology(data_dir, eval_file):
        f = open(eval_file, "r")
        protein_topology_list = []
        for protein_split_name in f:
            protein = protein_split_name[:4] # protein name extracted from the protein split to extract the topology
            topology_dir = os.path.join(data_dir, protein, "result", "output_" + protein + ".h5")
            
    
# gets list of protein names/splits that were evaluated
def get_eval_list(eval_file):
    f = open(eval_file, "r")
    eval_list = [protein for protein in f]
    f.close()
    return get_eval_list


    
        
        
        

# eval_protein_dir = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/splits/protein/eval.txt"
# eval_protein_split = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/protein_split/2JLE100"

eval_protein_h5 = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/proteins/2JLE/result/output_2JLE.h5"
# eval_protein_h5 = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/proteins/output_1Y6R/result/output_1Y6R.h5"


# pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/2JLE100_20_2500rollout.pt', map_location=pt.device('cpu'))

pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/2JLE100_20_2000rollout.pt', map_location=pt.device('cpu'))

original_topology = load_frame(eval_protein_h5, 1)
# original_topology = load_hdf5(eval_protein_h5)

print("2JLE Topology: ", original_topology.topology)
print("2JLE xyz shape: ", original_topology.xyz.shape)

evaluated_traj = pickle["rollout_u_pos"]

transposed_traj = np.transpose(evaluated_traj, (1, 0, 2))
transposed_traj = transposed_traj.detach().cpu().numpy()
processed_traj = np.divide(transposed_traj, 10)

print("2JLE100 Evaluated Rollout Position Shape: ", processed_traj.shape)


# traj = Trajectory(processed_traj, original_topology.topology)
# print(traj)

# traj.save("test.pdb")