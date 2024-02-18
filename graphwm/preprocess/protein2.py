import numpy as np
import os
import shutil
from pathlib import Path
from p_tqdm import p_umap

import mdtraj


def load_protein_traj(data_dir):
    """
    atom_coords: float, (num_traj, num_atom, 3)
    atom_types: int, (num_atom,)
    bond_indices: int, (num_bonds, 2)
    """
    full_path = data_dir + '/result/output_' + data_dir[len(data_dir) - 4:] + '.h5'
    traj = mdtraj.load(full_path)
    # scale units to A
    print('read traj')
    atom_coords = traj.xyz * 10
    lattices = traj.unitcell_lengths * 10.

    # get topology of the graph
    table, bonds = traj.top.to_dataframe()

    atom_types = [atom.element.atomic_number for atom in traj.top.atoms]
    atom_types = np.array(atom_types)

    # shape (num_bonds, 2), denotes indices of atoms connected by bonds
    bond_indices = bonds[:, :2]
    #rgs = 0

    #assert atom_coords.shape[0] == lattices.shape[0] == rgs.shape[0]
    #import pdb; pdb.set_trace()
    print(atom_coords, atom_types, bond_indices)
    return [atom_coords, atom_types, bond_indices]

def split_protein_traj(data_dir, data_save_dir, nsplit=None, traj_len=200):

    full_path = data_dir + '/result/output_' + data_dir[len(data_dir) - 4:] + '.h5'
    traj = mdtraj.load(full_path)
    
    n = traj.n_frames
    max_s = n - traj_len
    if(nsplit == None):
        nsplit = int(n / traj_len + 1) * 2
    print(n, nsplit)
    idx = np.arange(0, max_s, traj_len)
    idx = np.append(idx, [max_s])
    
    
    needed_split = nsplit - len(idx) 
    
    lw = np.arange(0, max_s, int(max_s / needed_split))
    hg = lw[1:]
    hg[-1] = max_s
    lw = lw[:-1]
    idx = np.append(idx, np.random.randint(lw, hg))
    print(idx, n, max_s, lw, hg)
    
    splits = []
    print('1')
    for i in idx:
        splits.append(traj[i:i+traj_len])
    print('2')
    
    dir_name = Path(data_dir).name
    save_path = Path(data_save_dir)
    
    
    def save_one_split(idx, split):
        print(f'start split {idx}')
        p = Path(os.path.join(save_path, dir_name + str(idx)))
        p.mkdir(exist_ok=True)
        split.save_dcd(os.path.join(p, 'trace.dcd'))
        shutil.copy(full_path,
                    os.path.join(p, 'bstate.pdb'))
        print(f'finished split {idx}')
    
    print('3')
    #splits = np.array(splits)
    i = np.arange(len(splits))
    print('4')
    
    
    save_one_split(0, splits[0])
    p_umap(save_one_split, i, splits)

def train_protein_test_split(data_dir, data_save_dir, n_split=0.9):
    traj = mdtraj.load(os.path.join(data_dir, '???.h5'),
                       top=os.path.join(data_dir, 'bstate.pdb'))
    
    n = traj.n_frames

    split_idx = int(n*(1-n_split))
    train_traj = traj[:split_idx]
    test_traj = traj[split_idx:]

    dir_name = Path(data_dir).name
    save_path = Path(data_save_dir)

    p = Path(os.path.join(save_path, dir_name + '_train'))
    p.mkdir(exist_ok=True)
    train_traj.save_dcd(os.path.join(p, 'trace.dcd'))
    shutil.copy(os.path.join(data_dir, 'bstate.pdb'),
                os.path.join(p, 'bstate.pdb'))

    p = Path(os.path.join(save_path, dir_name + '_test'))
    p.mkdir(exist_ok=True)
    test_traj.save_dcd(os.path.join(p, 'trace.dcd'))
    shutil.copy(os.path.join(data_dir, 'bstate.pdb'),
                os.path.join(p, 'bstate.pdb'))


split_protein_traj("./../datasets/proteins/1EBY", "./testing")