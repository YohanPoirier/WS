from __future__ import division
from sys import path
path.insert(1,r'D:\2 - Projects\Weak_Scatterer\Weak Scatterer\Weak Scatterer')
from WSC_Sim import *

# This script runs WSC_Sim.py.

# Current working directory.
cwd = getcwd()

# Run WSC_Sim.py.
WSC_Sim(cwd)