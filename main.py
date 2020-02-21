import numpy as np
import math
from matplotlib import pyplot as plt

# Line parameter definition
C = 100E-6
L = 250E-9
R = 0.05
G = 0.025
Zo = math.sqrt(L / C)
vo = 1.0 / math.sqrt(L * C)
d = 0.5
ncells = 100

cfln = 0.99
simtime = 100


# compute the discretization
dx = d / ncells
dt = cfln * dx / vo
nx = ncells + 1
nt = math.floor(simtime / dt)

# Internal Resistance
Rg = 50
RL = 0
R

# Voltage Source spec
Vamp = 1

# initialize line voltages and currents to zero:
V = np.zeros(nx, order='F')
I = np.zeros(nx - 1, order='F')

## TODO: Create Class for all respective update equations

# Set time step to 1 for initial program
n = 1
def voltage_update():
    """
     Update interior node voltages
     compute the multiplier coefficient:
    """
    cv = dt / (C * dx)

    # recursively update the line voltage at all interior nodes
    for k in range(2, nx-1):
        V[k] = V[k] - cv * (I[k] - I[k-1])
        print(V[k])


def voltage_source_update(n, dt):
    """
    In-line function to update the source voltage V(1)
    """
    Vs = voltage_source(n, dt)

    if Rg > 0:
        # compute the multiplier coefficients (from (2.62))
        b1 = C*dx*0.5/dt
        b2 = 0.5/Rg
        c1 = 1.0/(b1 + b2)
        c2 = b1 - b2
        V[1] = c1 * (c2 * V[1] - I[1] + Vs/Rg)

    else:
        V[1] = Vs

def load_voltage_update():
    if RL == 0:
        V[nx] = 0.
    else:
        b1 = C*dx*0.5/dt
        b2 = 0.5/RL
        c1 = 1.0/(b1 + b2)
        c2 = b1 - b2

        V[nx] = c1 * (c2 * V[nx] + I[nx-1])

def current_update():
    """
    In-line function to compute the update of the line current
    interior to the homogeneous line
    """
    # Compute the multiplier coefficient:
    ci = dt / (L * dx)
    for k in range(1, nx-1):
        I[k] = I[k] - ci * (V[k+1] - V[k])
        print(I[k])




def voltage_source(n, dt, source_type = 'unit_step'):
    if source_type == 'unit_step':
        if n > 0:
            Vs = Vamp


for n in range(0,10):
    voltage_update()
    current_update()







