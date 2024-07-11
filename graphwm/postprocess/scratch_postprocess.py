import torch as pt
# from mdtraj import load_frame, Trajectory, load_hdf5
import mdtraj as md
import numpy as np
# import pandas as pd
# from Bio.PDB.PDBParser import PDBParser
# import utils


# eval_protein_dir = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/splits/protein/eval.txt"
eval_protein_split = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/protein_split_3/2JLE100"

eval_protein_h5 = "/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/proteins/do_not_delete/2JLE/result/output_2JLE.h5"

pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/2JLE100_seed20_length3000.pt', map_location=pt.device('cpu'))

original_topology = md.load_frame(eval_protein_h5,5)

# original_topology = original_topology.atom_slice(original_topology.top.select('not name NA or not name CL'))

# Define the list of ion residue names you want to remove
ion_residue_names = ['NA', 'CL', 'MG', 'CA', 'K']

# Select all atoms that are not part of the ion residues
non_ion_atoms = original_topology.topology.select(' or '.join([f'resname != {res}' for res in ion_residue_names]))

# Slice the trajectory to keep only non-ion atoms
original_topology = original_topology.atom_slice(non_ion_atoms)

# original_topology = load_hdf5(eval_protein_h5)
# print("Shape of original file xyz: ", original_topology.xyz.shape)
print("Shape of original file topology: ", original_topology.topology)
for x in pickle:
    print(x)
    # print(x, pickle[x])
# print("rollout_u_pos", pickle["rollout_u_pos"].shape)
# print("cg_u_pos", pickle["cg_u_pos"].shape)
print("cg_bonds", pickle["cg_bonds"].shape)
# print("cg_weights", pickle["cg_weights"])
# print("n_cg_bond", pickle["n_cg_bond"])
# print("cg_bonds", pickle["cg_bonds"])
print("cluster", pickle["cluster"], pickle["cluster"].shape)


evaluated_traj = pickle["rollout_u_pos"]
cg_bonds = pickle['cg_bonds']

transposed_traj = np.transpose(evaluated_traj, (1, 0, 2))
transposed_traj = transposed_traj.detach().cpu().numpy()
processed_traj = np.divide(transposed_traj, 10)
# processed traj shape should be [rollout count, bead count, position]
# for 2JLE100 w/ 2000 rollouts: 20, 551, 3 

# carbon
top = md.Topology()
chain = top.add_chain()
residue = top.add_residue('ACE', chain)

element = md.element.carbon

frame_number = 5

for i, x in enumerate(processed_traj[frame_number]):
    atom_name = f'C{i}'  # Naming atoms C1, C2, C3, etc.
    top.add_atom(atom_name, element, residue)
    
for i, x in enumerate(cg_bonds):
    # print(x)
    atom1 = top.atom(x[0])
    atom2 = top.atom(x[1])
    top.add_bond(atom1, atom2)

# create a pdb from the first frame
traj = md.Trajectory(processed_traj[frame_number], top)

print(traj)



traj.save("updated_2JLE100.pdb")

print("")

particle_types = pickle["particle_types"]
print(pt.min(particle_types))
print(pt.min(pickle["cluster"]))

# 176 is the cg number of one of the stretched bond beads
# to do: find the atoms that are in that bead
atoms_bead_176 = [] # not too far (has 3 chlorine atoms)
atoms_bead_142 = [] # a little far (has 6 sodium atoms)
atoms_bead_114 = [] # close to main cluster (2 sodium, 5 chlorine)
atoms_bead_391 = [] # within main cluster (no ions)
atoms_bead_187 = []
atoms_bead_10 = []

for i, x in enumerate(pickle["cluster"]):
    # if(x == 176):
    #     atoms_bead_176.append(particle_types[i])
    # if(x == 142):
    #     atoms_bead_142.append(particle_types[i])
    # if(x == 114):
    #     atoms_bead_114.append(particle_types[i])
    # if(x == 391):
    #     atoms_bead_391.append(particle_types[i])
    if(x == 187):
        atoms_bead_187.append(particle_types[i])
    if(x == 10):
        atoms_bead_10.append(particle_types[i])

# print("far but not too far bead: ", atoms_bead_176, "\n")
# print("somewhat far bead: ", atoms_bead_142, "\n")
# print("close to cluster bead: ", atoms_bead_114, "\n")
# print("within range bead: ", atoms_bead_391, "\n")
print(atoms_bead_187)
print(atoms_bead_10)



# take the first rollout for now
# for x in processed_traj[0]:
#     top.add_atom()
    


# print("Shape of transposed evaluated positions 2JLE100: ", processed_traj)


# print("Evaluated topology value", evaluated_traj)
# 2JLE101

# traj = Trajectory(processed_traj, original_topology.topology)
# print(traj)

# traj.save("test.pdb")