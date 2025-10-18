"""
main.py
Author: Thomas Sutter
Description:
    Main script for lattice gas simulation scripts
"""

import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import Voronoi, voronoi_plot_2d
import scipy.special as scispec
import time
import tifffile
# Generates a "diffraction pattern" based on the crystal structure
# Expects an array of shape [# of particles, 2], where the axis=1 determines x/y in a cartesian coordinate system


def plot_reciprocal_lattice(cart_positions, q_range = 8, q_vals = 250, plot = True, fname = 'I_data'):
    
    x0, y0 = np.mean(cart_positions, axis = 0)
    
    # ~ condition = np.sqrt((cart_positions[:, 0] - x0)**2 +  (cart_positions[:, 1] - y0)**2) < 25
    # ~ cart_positions = cart_positions[condition]
    
    # ~ points_to_draw = 2*np.pi/np.sqrt(3)*np.array([[np.sqrt(3), -1],[-np.sqrt(3), +1],[0, 2],[0, -2],[np.sqrt(3), 1], [-np.sqrt(3), -1]])
    
    q1 = np.linspace(-q_range, q_range, q_vals)
    Q = np.meshgrid(q1,q1)
    
    f_data = np.zeros(np.shape(Q[0]),dtype='complex128')
    
    for n in range(len(cart_positions)):
        if n % int(0.1*len(cart_positions) + 1) == 0:
            print(f'{int(100*n/len(cart_positions))} %')
        f_data += np.exp(1.j*(cart_positions[n,0]*Q[0]+cart_positions[n,1]*Q[1]))
        
    I_data = np.abs(f_data)**2
    phi_data = np.angle(f_data)
    
    tifffile.imwrite(f'diffraction//{fname}.tiff', np.array(I_data, dtype = np.float32))
    

    fig, axs = plt.subplots(1, 2, figsize = (8, 4))
    axs[0].imshow(I_data,extent=[min(q1),max(q1),min(q1),max(q1)], vmax = 1e6, vmin = 0)
    
    axs[0].text(0,0,r'$\Gamma$',fontsize=20)
    axs[0].scatter([0],[0],c='red', alpha = 0.5)
        # ~ for i in range(len(points_to_draw)):
            # ~ ax.scatter(points_to_draw[i,0],points_to_draw[i,1],c='red', alpha = 0.5)
    
    # ~ axs[1].set_aspect('equal')
    axs[1].scatter(cart_positions[:, 0], cart_positions[:, 1])
    axs[1].set_xlim([min(cart_positions[:, 0]), max(cart_positions[:, 0])])
    axs[1].set_ylim([min(cart_positions[:, 1]), max(cart_positions[:, 1])])
        
    if plot == True:
        plt.show()
    else:
        plt.savefig(f'real_space//{fname}.png') #, dpi = 250)





def go():
    x_list = []
    y_list = []
    theta_list = []
    theta = 0

    
    N = 100
    M = 100
    
    for i in range(N):
        
        for j in range(M):
            x = i + j * np.sin(theta)
            y = -j * np.cos(theta)
        
            x_list.append(x)
            y_list.append(y)
        
        # move plane angle
        theta += np.random.normal(loc = 0, scale = 0.0100)

        
    cart_positions = np.array([x_list, y_list]).transpose()
    plot_reciprocal_lattice(cart_positions, plot = True)
    # ~ plt.scatter(x_list, y_list)
    # ~ plt.show()


go()
