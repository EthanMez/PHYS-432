"""
PHYS 432 Assignment 3, Question 3
Author: Ethan Meszaros (adopted from code by Prof. Eve Lee)

Lava flow simulation. 
"""
import numpy as np
import matplotlib.pyplot as plt

#grid and step parameters 
Ngrid = 100
Nsteps = 5000
dt = 1
dy = 1

#setting up constants used in the matrix
g = 9.8 #accounting for gravity
al = np.pi/4

D = 0.1
beta = D*dt/dy**2
C = g*np.sin(al) #the extra constant due to the body force

#creating grid on which values will be updated 
y = np.arange(0, Ngrid*1., dy) / Ngrid #as per class code, 1. is to make sure floats are used, not integers

f = np.zeros(len(y))
f[int(Ngrid/2)] = 1

#setting up plot
plt.ion()
fig, ax = plt.subplots(1,1)

#initial state before animation for reference
ax.plot(f, y, 'k-')

#plot that will be updated
plt, = ax.plot(f, y, 'ro')

#axis limits
ax.set_xlim([0,1])
ax.set_ylim([0,2])

fig.canvas.draw()


# COULD NOT GET FOLLOWING TO WORK

'''
for ct in range(Nsteps):

    #matrix that solves for updated function 
    A1 = (np.eye(Ngrid) * (1.0 + 2.0 * beta1) + np.eye(Ngrid, k=1) * -beta1 + np.eye(Ngrid, k=-1) * -beta1 
          + np.eye(Ngrid, k=-1) * dt*C) #adding extra constant

    #using A1 to update function
    f = np.linalg.solve(A1, f)

    # update the plot
    plt.set_ydata(f)

    fig.canvas.draw()
    pl.pause(0.001)
'''
