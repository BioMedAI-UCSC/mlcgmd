import mdtraj as md

traj = md.load_hdf5("./../datasets/protein_split/1Y6R0/bond.h5")
print(traj)