"""
clean up and rename to polymer.py in the end.
"""
import os
import sys
import time
from pathlib import Path
import multiprocessing as mp
from p_tqdm import p_umap
import mdtraj as md
import numpy as np

from battery import load_battery_data
from chain import load_polymer_rg
from hbv import load_hbv_traj
from hbv import split_hbv_traj
from hbv import train_test_split
from protein import load_protein_traj
from protein import split_protein_traj
from protein import protein_train_test_split

# from ..data.utils import store_data
from utils import store_data

def polymer_to_h5(data_dir, data_save_dir):
  """
  save whole trajectory data directly from .txt polymer trajectory.
  """
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print('poly_file_dirs', poly_file_dirs)
  print(f"Found {len(poly_file_dirs)} polymer trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def process_one_file(poly_file):
    poly_index = poly_file.parts[-1]
    os.makedirs(os.path.join(data_save_dir, poly_index), exist_ok=True)
    if not Path(str(os.path.join(data_save_dir, poly_index, 'bond.h5'))).exists():
      try:
        data = load_polymer_rg(poly_file)
        store_data(['position'], [data[0]], os.path.join(data_save_dir, poly_index, 'position.h5'))
        store_data(['lattice'], [data[1]], os.path.join(data_save_dir, poly_index, 'lattice.h5'))
        store_data(['rgs'], [data[2]], os.path.join(data_save_dir, poly_index, 'rgs.h5'))
        store_data(['particle_type'], [data[3]], os.path.join(data_save_dir, poly_index, 'ptype.h5'))
        store_data(['bond_indices'], [data[4]], os.path.join(data_save_dir, poly_index, 'bond.h5'))
      except Exception as e:
        print(poly_index)
        print(e)
        pass

  now = time.time()
  process_one_file(poly_file_dirs[0])
  p_umap(process_one_file, poly_file_dirs)
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

def battery_to_h5(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} polymer trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")
  
  def process_one_file(poly_file):
    poly_index = poly_file.parts[-1]
    os.makedirs(os.path.join(data_save_dir, poly_index), exist_ok=True)
    if not Path(str(os.path.join(data_save_dir, poly_index, 'diffusivity.h5'))).exists():
      try:
        data = load_battery_data(poly_file)
        store_data(['wrapped_position'], [data[0]], os.path.join(data_save_dir, poly_index, 'wrapped_position.h5'))
        store_data(['unwrapped_position'], [data[1]], os.path.join(data_save_dir, poly_index, 'unwrapped_position.h5'))
        store_data(['lattice'], [data[2]], os.path.join(data_save_dir, poly_index, 'lattice.h5'))
        store_data(['raw_particle_type'], [data[3]], os.path.join(data_save_dir, poly_index, 'raw_ptype.h5'))
        store_data(['particle_type'], [data[4]], os.path.join(data_save_dir, poly_index, 'ptype.h5'))
        store_data(['bond_indices'], [data[5]], os.path.join(data_save_dir, poly_index, 'bond.h5'))
        store_data(['bond_type'], [data[6]], os.path.join(data_save_dir, poly_index, 'bond_type.h5'))
        store_data(['diffusivity'], [data[7]], os.path.join(data_save_dir, poly_index, 'diffusivity.h5'))
      except OSError:
        pass
    
  now = time.time()
  p_umap(process_one_file, poly_file_dirs)  
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")


def protein_to_h5(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} protein trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def process_one_file(poly_file):
    poly_index = poly_file.parts[-1]
    os.makedirs(os.path.join(data_save_dir, poly_index), exist_ok=True)
    if not Path(str(os.path.join(data_save_dir, poly_index, 'bond.h5'))).exists():
      try:
        data = load_protein_traj(poly_file)
        print('successfully loaded')
        store_data(['position'], [data[0]], os.path.join(data_save_dir, poly_index, 'position.h5'))
        store_data(['particle_type'], [data[1]], os.path.join(data_save_dir, poly_index, 'ptype.h5'))
        store_data(['bond_indices'], [data[2]], os.path.join(data_save_dir, poly_index, 'bond.h5'))
      except Exception as e:
        print(poly_index)
        print(e)
        pass
    
  now = time.time()
  process_one_file(poly_file_dirs[0])
  p_umap(process_one_file, poly_file_dirs)  
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

def split_protein(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} protein trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def process_one_file(poly_file):
    split_protein_traj(poly_file, data_save_dir)

  now = time.time()
  for f in poly_file_dirs:
    process_one_file(f)
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

def protein_train_test_split(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  prot_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(prot_file_dirs)} protein trajectories. First : {prot_file_dirs[0]}")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def processOneFile(protein):
    protein_train_test_split(str(protein), str(data_save_dir))

  now = time.time()
  for f in prot_file_dirs:
    print(f"processing {str(f)[len(str(f)) - 4:]}")
    processOneFile(f)
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(prot_file_dirs)} || Time Elapsed: {elapsed}")


def hbv_to_h5(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} hbv trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def process_one_file(poly_file):
    poly_index = poly_file.parts[-1]
    os.makedirs(os.path.join(data_save_dir, poly_index), exist_ok=True)
    if not Path(str(os.path.join(data_save_dir, poly_index, 'bond.h5'))).exists():
      try:
        data = load_hbv_traj(poly_file)
        print('successfully loaded')
        store_data(['position'], [data[0]], os.path.join(data_save_dir, poly_index, 'position.h5'))
        store_data(['particle_type'], [data[1]], os.path.join(data_save_dir, poly_index, 'ptype.h5'))
        store_data(['bond_indices'], [data[2]], os.path.join(data_save_dir, poly_index, 'bond.h5'))
      except Exception as e:
        print(poly_index)
        print(e)
        pass

  now = time.time()
  process_one_file(poly_file_dirs[0])
  p_umap(process_one_file, poly_file_dirs)  
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

def split_hbv(data_dir, data_save_dir):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} hbv trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  def process_one_file(poly_file):
    split_hbv_traj(poly_file, data_save_dir)

  now = time.time()
  for f in poly_file_dirs:
    process_one_file(f)
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

def make_test(data_dir, data_save_dir, dilation = 5, seq_len = 20):
  input_path = Path(data_dir)
  output_path = Path(data_save_dir)
  traj = md.load_dcd(input_path/'trace.dcd', top=input_path/'bstate.pdb')
  idx = np.arange(0, seq_len*dilation, dilation)
  print(idx)
  traj2 = traj.slice(idx)
  print(input_path)
  output_path = output_path/Path(input_path.name + '_d' + str(dilation) + '_sql' + str(seq_len))
  output_path.mkdir(exist_ok=True)
  print(output_path)

  traj2.save_dcd(output_path / Path('trace.dcd'))

def hbv_test_split(data_dir, data_save_dir, n_split=0.9):
  data_dir = Path(data_dir)
  poly_file_dirs = [d for d in list(data_dir.iterdir()) if os.path.isdir(d)]
  print(f"Found {len(poly_file_dirs)} hbv trajectories.")
  print(f"Use {mp.cpu_count()} cores.")
  print("Start processing...")

  now = time.time()
  train_test_split(poly_file_dirs[0], data_save_dir, n_split)
  p_umap(train_test_split, poly_file_dirs, data_save_dir, n_split)  
  elapsed = time.time() - now
  print(f"Done. Number of rollouts: {len(poly_file_dirs)} || Time Elapsed: {elapsed}")

if __name__ == '__main__':
    import sys
    print(sys.prefix)
    
    if len(sys.argv) < 2:
        print('enter command')
        sys.exit(1)
    command = sys.argv[1]
    if command == 'chain':
      data_dir, data_save_dir = sys.argv[2:]
      polymer_to_h5(data_dir, data_save_dir)
    elif command == 'battery':
      data_dir, data_save_dir = sys.argv[2:]
      battery_to_h5(data_dir, data_save_dir)
    elif command == 'hbv': # hepatitis B virus
      data_dir, data_save_dir = sys.argv[2:]
      hbv_to_h5(data_dir, data_save_dir)
    elif command == 'protein': # protein
      data_dir, data_save_dir = sys.argv[2:]
      protein_to_h5(data_dir, data_save_dir)
    elif command == 'ttprot': #train test split proteins
      data_dir, data_save_dir = sys.argv[2:]
      protein_train_test_split(data_dir, data_save_dir)
    elif command == 'split_protein':
      data_dir, data_save_dir = sys.argv[2:]
      split_protein(data_dir, data_save_dir)
    elif command == 'split_hbv':
      data_dir, data_save_dir = sys.argv[2:]
      split_hbv(data_dir, data_save_dir)
    elif command == 'make_test':
      data_dir, data_save_dir = sys.argv[2:4]
      dilation = 5
      seq_len = 40
      if(len(sys.argv) >= 5):
        dilation = int(sys.argv[4])
      if(len(sys.argv) == 6):
        seq_len = int(sys.argv[5])
      make_test(data_dir,data_save_dir,dilation,seq_len)
    elif command == 'hbv_test_split':
      data_dir, data_save_dir = sys.argv[2:4]
      n_split = 0.9
      if(len(sys.argv) == 5):
        n_split= int(sys.argv[4])
      hbv_test_split(data_dir, data_save_dir, n_split)