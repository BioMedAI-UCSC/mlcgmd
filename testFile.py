import torch as pt
from graphwm.common import log_hyperparameters, PROJECT_ROOT

print(pt.cuda.is_available())
print(PROJECT_ROOT)