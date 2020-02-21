import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib

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
print('Running to time step sample {} at spacial sample distance of {}'.format(nt, nx))

# Internal Resistance
Rg = 50
RL = 10

# Voltage Source spec
Vamp = 1

# initialize line voltages and currents to zero:
V = np.zeros(nx, order='F')
I = np.zeros(nx - 1, order='F')

# TODO: Create Class for all respective update equations

# Set time step to 1 for initial program


def voltage_update(n):
    """
     Update interior node voltages
     compute the multiplier coefficient:
    """
    global t
    t = dt*(n-0.5)
    print(t)

    cv = dt / (C * dx)


    # recursively update the line voltage at all interior nodes
    for k in range(0, nx - 1):
        V[k] = V[k] - cv * (I[k] - I[k - 1])

def voltage_source_update(n, dt):
    """
    In-line function to update the source voltage V(1)
    """
    Vs = voltage_source(n)

    if Rg > 0:
        # compute the multiplier coefficients (from (2.62))
        b1 = (C * dx * 0.5) / dt
        b2 = 0.5 / Rg
        c1 = 1.0 / (b1 + b2)
        c2 = b1 - b2
        V[0] = c1 * (c2 * V[0] - I[0] + (Vs / Rg))

    else:
        V[0] = Vs


def voltage_load_update(n, RL):
    if RL == 0:
        V[nx] = 0.
    else:
        b1 = C * dx * 0.5 / dt
        b2 = 0.5 / RL
        c1 = 1.0 / (b1 + b2)
        c2 = b1 - b2

        V[n] = c1 * (c2 * V[n] + I[n-1])


def current_update():
    """
    In-line function to compute the update of the line current
    interior to the homogeneous line
    """
    # Compute the multiplier coefficient:
    ci = dt / (L * dx)
    for k in range(0, nx-1):
        I[k] = I[k] - ci * (V[k + 1] - V[k])
        print(I[k])


def voltage_source(n, source_type='unit_step'):
    Vs = 0
    if source_type == 'unit_step':
        if n >= 0:
            Vs = Vamp
    else:
        Vs = 0

    return Vs



fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
for n in range(1, nt):

    x = np.arange(0, d, nx)
    ax.cla()
    print(V)
    voltage_update(n)
    voltage_source_update(n, dt)
    voltage_load_update(n-1, RL)

    current_update()
    line.set_data(x, V)

    print(V)
    plt.show()
    plt.pause(1)

