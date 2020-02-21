import numpy as np
import math
from matplotlib import pyplot as plt

# Line parameter definition
C = 100E-6
L = 250E-9
R = 0.05
G = 0.025
Zo = math.sqrt(L/C)
vo = 1.0/math.sqrt(L*C)
d = 0.5
ncells = 100
cfln = 0.99
simtime = 100

# compute the discretization
dx = d / ncells
dt = cfln * dx / vo
nx = ncells + 1
nt = math.floor(simtime/dt)

# initialize line voltages and currents to zero:
V = np.zeros(nx, order='C')
I = np.zeros(nx-1, order='F')

class txLineFDTD():
    def vsupdate(self):
        t = dt*(n-0.5)





