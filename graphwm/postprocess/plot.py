import plotly.express as px
import pandas as pd
import torch as pt
import numpy as np

# df = px.data.iris()

pickle = pt.load('/projects/bbpa/coarseGrained/mlcgmd/mlcgmd_models/protein_pnr/nsteps5_stepsize_0.0001/2JLE100_20_2000rollout.pt', map_location=pt.device('cpu'))

evaluated_traj = pickle["rollout_u_pos"]

transposed_traj = np.transpose(evaluated_traj, (1, 0, 2))
transposed_traj = transposed_traj.detach().cpu().numpy()
processed_traj = np.divide(transposed_traj, 10)

graph_traj_format = np.transpose(processed_traj, (0, 2, 1))

df = pd.DataFrame({'x':graph_traj_format[0][0], 'y':graph_traj_format[0][1], 'z':graph_traj_format[0][2]})

fig = px.scatter_3d(df, x='x', y='y', z='z')
fig.show()