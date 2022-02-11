'''
Author: Ethan Meszaros
Student ID: 260919728

PHYS 432 Assignment 2, Q3

Simulation of "leap frogging" between two vortices. 
'''

import numpy as np
import matplotlib.pyplot as pl

N_steps = 100
dt = 10

#setup initial positions 
y_v = np.array([8, 4, 10, 2])
x_v = np.array([2, 2, 8, 8])
k_v = np.array([1, -1, 1, -1])

#Set up the plot
pl.ion() #turn interactive mode on, which updates the plots after each command
fig, ax = pl.subplots(1, 1)
p, = ax.plot(x_v, y_v, 'k+', markersize = 10)

#create a grid that contains initial velocity streamlines
ngrid = 20
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
vel_x = np.zeros(np.shape(X))
vel_y = np.zeros(np.shape(Y))

#set limits of grid
ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

def velocity(x, y, y_i, x_i, k):
    '''(flt, flt, flt, flt, flt) -> flt, flt
    
    A function that determines the x and y velocities at a point (x, y)
    due to a vortex located at (x_i, y_i) that contains constant k. 
    '''
    
    y = y - y_i
    x = x - x_i
    
    #vortex modelled as line vortex in Cartesian
    u_x = -k*y/(x**2 + y**2) 
    u_y = k*x/(x**2 + y**2)
    return u_x, u_y

#add up the velocities due to each vortex at every point in space
for i in range(len(x_v)):
    v_x, v_y = velocity(X, Y, y_v[i], x_v[i], k_v[i])
    
    vel_x += v_x
    vel_y += v_y

stream = ax.streamplot(X, Y, vel_x, vel_y, arrowstyle = '->', density = [1,1])
     

#Evolution

count = 0 

#create empty arrays to be filled with updated velocities at each vortex    
U_x = np.array(np.zeros(4))
U_y = np.array(np.zeros(4))

while count < 4:
    for i in range(len(x_v)):
        #at the position of a vortex, we add the velocity contributions of 
        #the other vortices
        for j in range(len(x_v)):
            if i != j: #to prevent dividing by zero
                u_x, u_y = velocity(x_v[i], y_v[i], y_v[j], x_v[j], k_v[j])
            else:
                u_x, u_y = (0, 0)
            
            #adding the velocity contributions from other vortices 
            U_x[i] += u_x
            U_y[i] += u_y
        #updating the positions 
        x_v[i] += U_x[i]*dt
        y_v[i] += U_y[i]*dt
    
    ax.cla()
    
    p.set_xdata(x_v)
    p.set_ydata(y_v)
    
    ax.streamplot(X, Y, vel_x, vel_y, density = [1,1])
    
    fig.canvas.draw()
    pl.pause(0.001)
    count += 1
    
#Unfortunately, I couldn't finish the simulation
    
#y_v += U_y*dt 
#x_v += U_x*dt

#ax.plot(x_v, y_v, 'k+', markersize = 10)
    
    