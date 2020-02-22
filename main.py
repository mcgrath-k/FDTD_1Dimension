import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import animation

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
simtime = .1

# compute the discretization
dx = d / ncells
dt = cfln * dx / vo
nx = ncells + 1
nt = math.floor((simtime / dt)/1000)
print(nt)
print('Running to time step sample {} at spacial sample distance of {}'.format(nt, nx))

# Internal Resistance
Rg = 50
RL = 0

# Voltage Source spec
Vamp = 1.

# initialize line voltages and currents to zero:
V = np.zeros(nx + 1, order='F')
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

    cv = dt / (C * dx)

    # recursively update the line voltage at all interior nodes
    for k in range(0, nx - 1):
        V[k] = V[k] - cv * (I[k] - I[k - 1])

def voltage_source_update(n):
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


def voltage_load_update(RL):
    if RL == 0:
        V[nx] = 0
    else:
        b1 = (C * dx * 0.5) / dt
        b2 = 0.5 / RL
        c1 = 1.0 / (b1 + b2)
        c2 = b1 - b2

        V[nx-1] = c1 * (c2 * V[nx-1] + I[nx-2])



def current_update():
    """
    In-line function to compute the update of the line current
    interior to the homogeneous line
    """
    # Compute the multiplier coefficient:
    ci = dt / (L * dx)
    for k in range(0, nx-1):
        I[k] = I[k] - (ci * (V[k + 1] - V[k]))


def voltage_source(n, source_type='unit_step'):
    Vs = 0
    if source_type == 'unit_step':
        if n >= 0:
            Vs = Vamp
    else:
        Vs = 0

    return Vs


fig = plt.figure()
ax = plt.axes(xlim=(0,d), ylim=(-2,2))
line, = ax.plot([], [], lw=2)

def animate(n):
    x = np.linspace(0, d, nx+1)
    voltage_update(n)
    voltage_source_update(n)
    voltage_load_update(RL)

    current_update()
    print(I)
    print(V)
    print(len(x), len(V))
    line.set_data(x, V)
    return line,

def init():
    line.set_data([], [])
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames = 200, interval=20, blit=True)

plt.show()