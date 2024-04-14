import torch as pt
import mdtraj as md

pdb = md.load('/projects/bbpa/coarseGrained/mlcgmd/graphwm/datasets/proteins/1EBY/simulation/final_state.pdb')
pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/seed_15.pt', map_location=pt.device('cpu'))

# print("pdb keys", pdb)
print("pickle keys", pickle.keys())
# print("bonds", pickle["bonds"].shape)

# position is init positions
# print("position", pickle["position"].shape)

# rollout u pos is predicted trajectory positions
# print("rollout_u_pos", pickle["rollout_u_pos"].shape)
unimportant = ["cg_u_pos", "cg_weights", "cg_bonds", "n_cg_bond", "time_elapsed", "model_params"]

for key in pickle.keys():
    if key not in unimportant:
        print(key, ":", pickle[key])

# coordinates = np.array([[0.0, 0.0, 0.0],
#                         [1.0, 0.0, 0.0],
#                         [0.0, 1.0, 0.0]])

coordinates = pickle["rollout_u_pos"]

# # Atom names and residue numbers
# atom_names = ['CA', 'CA', 'CA']
# residue_numbers = [1, 2, 3]

# # Create a topology
# topology = md.Topology()
# chain = topology.add_chain()
# for res_num in residue_numbers:
#     residue = topology.add_residue('ALA', chain)
#     for atom_name in atom_names:
#         topology.add_atom(atom_name, md.element.carbon, residue)

# # Create a trajectory
# trajectory = md.Trajectory(coordinates, topology)

# # Save trajectory to PDB file
# trajectory.save('output.pdb')