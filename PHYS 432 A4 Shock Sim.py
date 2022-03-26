"""
PHYS 432 A4 Q2

Evolving adiabatic shocks 
(adapted from Prof. Eve Lee "phys432_w2022_sound_wave.py").
"""
import numpy as np
import matplotlib.pyplot as plt

# setting up grid, constants and time/space steps 
Ngrid = 100
Nsteps = 5000
dt = 0.1
dx = 2.0

# create sound speed grid, since sound speed changes across shock
cs2 = np.ones(Ngrid)
gamma = 5/3 # adiabatic index

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.ones(Ngrid) # rho x e_tot
u = np.zeros(Ngrid+1) # advective velocity 

def adv(f, u, dt, dx):
    '''
    Takes in a function f and updates it based on flux information. 
    '''
    # calculating flux terms
    J = np.zeros(len(f)+1)
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) # update

    return f

# Apply initial Gaussian perturbation at the left edge
Amp, sigma = 1, Ngrid/100
f1 = f1 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f2 = f2 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# create plot
plt.ion()
fig, (ax1, ax2) = plt.subplots(2,1)

# density plot
x1, = ax1.plot(x, f1, 'ro')
ax1.set_xlim([0, dx*Ngrid+1])
ax1.set_ylim([0.8, 2.1])
ax1.set_xlabel('x')
ax1.set_ylabel('Density')

M = (f2/f1)/cs2 # Mach number

# Mach number plot
x2, = ax2.plot(x, M, 'ro')
ax2.set_xlabel('x')
ax2.set_ylabel('Mach Number')


fig.canvas.draw()

# animation
for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))
    
    # advect density and momentum
    f1 = adv(f1, u, dt, dx)
    f2 = adv(f2, u, dt, dx)
    
    # define sound speed (based on energy conservation derived in class)
    cs2 = (gamma - 1) * ((f3/f1) - 0.5*u[1:]**2)
    
    # define pressure
    P = f1 * cs2 / gamma
    
    # add the source term for momentum, which contains a pressure gradient 
    f2[1:-1] = f2[1:-1] - (dt / dx) * (P[2:] - P[:-2]) / gamma
    
    # re-calculate advection velocity 
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))
    
    # advect energy 
    f3 = adv(f3, u, dt, dx)     
    
    # add energy source term, which contains pressure flux gradient
    f3[1:-1] = f3[1:-1] - (dt / dx) * (P[2:] * u[2:-1] - P[:-2] * u[1:-2]) / gamma
    
    # correct for source term at the boundary (reflective)
    f2[0] = f2[0] - 0.5 * (dt / dx) * cs2[0] * (f1[1] - f1[0])
    f2[-1] = f2[-1] - 0.5 * (dt / dx) * cs2[0] * (f1[-1] - f1[-2])
    
    # re-calculate Mach number
    M = (f2/f1)/cs2

    # update the plot
    x1.set_ydata(f1)
    x2.set_ydata(M)
    fig.canvas.draw()
    plt.pause(0.001)