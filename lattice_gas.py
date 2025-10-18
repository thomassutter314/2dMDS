"""
Description:
    Main script for lattice gas simulation scripts
"""

import numpy as np
from numba import njit
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay
import scipy.special as scispec
from scipy.interpolate import RBFInterpolator
import networkx as nx
import time
import tifffile
import os
import re


# ~ color_dict = {4 : '#008AAF', 5 : '#AFFF8F', 6 : 'black', 7 : '#DE8FFF', 8 : '#AF2500'}
color_dict = {5 : '#AFFF8F', 6 : 'black', 7 : '#DE8FFF'}


def extract_number(s):
    """Extracts the final number from a string, ignoring file extensions."""
    match = re.search(r'_(\d+)\.csv$', s)  # Extracts number before ".csv"
    return int(match.group(1)) if match else 0  # Convert to integer

@njit
def compute_kick(pos_1, pos_2, J_1 = 0.5, J_2 = 0.2, xi_1 = 2, xi_2 = 1, epsilon = 1e-3):
    D = pos_2 - pos_1
    R = np.sqrt(D[0]**2 + D[1]**2) + epsilon
    
    theta_ij = np.arctan2(D[1], D[0])

    return  -J_1 * D * np.exp(-R/xi_1) * (1/R**3 + 1/(xi_1 * R**2)) + J_2 * np.sin(6 * theta_ij) * np.array([-D[1], D[0]])/R * np.exp(-R/xi_2)

# J_2 = 0.005 used for paper
    
@njit
def time_evolve(par_positions, par_velocities, temperature, L, J_1, iterations_per_plot, dt = 0.05):
    indices = np.random.permutation(par_positions.shape[0])
    par_positions = par_positions[indices, :]
    par_velocities = par_velocities[indices, :]
    
    
    for i in range(iterations_per_plot):

        # ~ par_positions += 0.01*np.random.normal(loc = 0, scale = np.sqrt(temperature), size = par_positions.shape)
        if temperature > 0:
            par_velocities += 0.01 * np.random.normal(loc = 0, scale = np.sqrt(temperature), size = par_positions.shape)
        
        
        for n in range(par_positions.shape[0]):
            for m in range(par_positions.shape[0] - 1 - n):
                M = m + 1 + n
                
                kick = compute_kick(par_positions[n, :], par_positions[M, :], J_1 = J_1)
                
                par_velocities[n, :] += kick
                par_velocities[M, :] += -kick
            
            
            # Keep the particles inbounds 
            # ~ if par_positions[n, 0] > L:
                # ~ par_positions[n, 0] = L
                # ~ par_velocities[n, 0] *= -1
            # ~ if par_positions[n, 0] < 0:
                # ~ par_positions[n, 0] = 0
                # ~ par_velocities[n, 0] *= -1
            # ~ if par_positions[n, 1] > L:
                # ~ par_positions[n, 1] = L
                # ~ par_velocities[n, 1] *= -1
            # ~ if par_positions[n, 1] < 0:
                # ~ par_positions[n, 1] = 0
                # ~ par_velocities[n, 1] *= -1

            if np.sqrt((par_positions[n, 0] - 0.5 * L)**2 + (par_positions[n, 1] - 0.5 * L)**2) > 0.5 * L:
                theta = np.arctan2(par_positions[n, 1] - 0.5 * L, par_positions[n, 0] - 0.5 * L)
                par_velocities[n, 0] += -0.05 * np.cos(theta)
                par_velocities[n, 1] += -0.05 * np.sin(theta)
                

            # ~ if np.sqrt((par_positions[n, 0] - L/2)**2 + (par_positions[n, 1] - L/2)**2) > L/2:
                # ~ theta = np.arctan2(par_positions[n, 1] - L/2, par_positions[n, 0] - L/2)
                # ~ par_positions[n, 0] = L/2 * (1 + np.cos(theta))
                # ~ par_positions[n, 1] = L/2 * (1 + np.sin(theta))
        
        
        par_velocities += -0.1 * dt * par_velocities
        par_positions += dt * par_velocities
        
    return par_positions, par_velocities
   
@njit
def find_closest(X, x):
    X = np.asarray(X)  # Ensure X is a NumPy array
    idx = np.abs(X - x).argmin()  # Find the index of the closest value
    return idx  # Return the closest index

@njit
def compute_g6(par_positions, psi_6, n_bins = 100, L = 70):
    g6_r_sqr = np.linspace((L/n_bins)**2, L**2, n_bins)
    g6_r = np.sqrt(g6_r_sqr)
    
    g6 = np.zeros(len(g6_r), dtype = np.complex128)
    weights = np.zeros(len(g6_r))
    
    for i in range(par_positions.shape[0]):
        for j in range(par_positions.shape[0]):
            D = np.sqrt((par_positions[i, 0] - par_positions[j, 0])**2 + (par_positions[i, 1] - par_positions[j, 1])**2)
            
            idx = find_closest(g6_r, D)
            
            g6[idx] += (psi_6[i].conjugate() * psi_6[j])
            weights[idx] += 1
    
    g6[weights > 0] = g6[weights > 0]/weights[weights > 0]
    
    return g6[:-1], g6_r[:-1]
    
def locate_defects(par_pos):
    # Step 1: Compute Delaunay Triangulation
    tri = Delaunay(par_pos)
    neighbors = {i: set() for i in range(len(par_pos))}
    
    # Step 2: Count Neighbors
    for simplex in tri.simplices:
        for i in range(3):
            neighbors[simplex[i]].add(simplex[(i+1)%3])
            neighbors[simplex[i]].add(simplex[(i+2)%3])
    
    # Identify defects
    five_sites = set(i for i in neighbors if len(neighbors[i]) == 5)
    seven_sites = set(i for i in neighbors if len(neighbors[i]) == 7)

    # Step 3: Construct Graph of Defects
    G = nx.Graph()
    for i in five_sites:
        for j in seven_sites:
            if j in neighbors[i]:  # They are adjacent
                G.add_edge(i, j)
    
    # Step 4: Solve Min-Cost Perfect Matching (Blossom Algorithm)
    matching = nx.max_weight_matching(G, maxcardinality=True)
    # Convert set to list of tuples
    matching = list(matching)
    matching_flat = [item for tup in matching for item in tup]
    
    # Determine the lone sites
    defects = five_sites | seven_sites # | is union of two sets
    lone = []
    for i in defects:
        if not i in matching_flat:
            lone.append(i)
    
    # Measure psi_6
    psi_6 = np.zeros(par_pos.shape[0], dtype = np.complex128)
    for i in range(len(psi_6)):
        for j in neighbors[i]:
            theta_i = np.arctan2(par_pos[j, 0] - par_pos[i, 0], par_pos[j, 1] - par_pos[i, 1])
            psi_6[i] += 1/len(neighbors[i]) * np.exp(6j * theta_i)
        
    
    return matching, lone, psi_6
    
def count_defects_in_dir(mainDir):
    dir_list = os.listdir(mainDir)
    dir_list = [f for f in dir_list if '.csv' in f]
    
    dir_list = sorted(dir_list, key=extract_number)
    dir_list = dir_list[-30:]
    
    dislocation_count = []
    disclination_count = []
    
    for i in range(len(dir_list)):
        file_loc = f'{mainDir}//{dir_list[i]}'
        par_positions = np.loadtxt(file_loc, delimiter = ',')
        
        # ~ cond = np.sqrt((par_positions[:, 0] - 101/2)**2 + (par_positions[:, 1] - 101/2)**2) < 40
        # ~ par_positions = par_positions[cond]
        
        bonds, lone, psi_6 = locate_defects(par_positions)
        
        dislocation_count.append(len(bonds))
        disclination_count.append(len(lone))
        
    return np.mean(dislocation_count), np.mean(disclination_count), np.abs(np.mean(psi_6))

def plot_bond(p1, p2, ax):    
    x1, y1 = p1
    x3, y3 = p2
    
    x2 = 0.5*(x1 + x3)
    y2 = 0.5*(y1 + y3)
    
    ax.plot([x1, x2], [y1, y2], c = color_dict[5], linewidth = 4, zorder = 1, alpha = 1, solid_capstyle='butt')
    ax.plot([x2, x3], [y2, y3], c = color_dict[7], linewidth = 4, zorder = 1, alpha = 1, solid_capstyle='butt')          
        
    ax.scatter([x1], [y1], facecolor = color_dict[5], edgecolor = 'k', s = 50)
    ax.scatter([x3], [y3], facecolor = color_dict[7], edgecolor = 'k', s = 50)
    

# Generates a "diffraction pattern" based on the crystal structure
# Expects an array of shape [# of particles, 2], where the axis=1 determines x/y in a cartesian coordinate system
def plot_reciprocal_lattice(cart_positions, q_range = 3, q_vals = 250, plot = True, fname = 'I_data'):
    
    x0, y0 = np.mean(cart_positions, axis = 0)
    
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
    

    fig, axs = plt.subplots(1, 3, figsize = (12, 6))
    axs[0].imshow(I_data,extent=[min(q1),max(q1),min(q1),max(q1)], vmax = 1e4, vmin = 0)
    
    
    vor = Voronoi(cart_positions)
    tri = Delaunay(cart_positions)
    
    neighbors = {i: set() for i in range(len(cart_positions))}
    # Count Neighbors
    for simplex in tri.simplices:
        for i in range(3):
            neighbors[simplex[i]].add(simplex[(i+1)%3])
            neighbors[simplex[i]].add(simplex[(i+2)%3])
            
    vertex_count = np.zeros(cart_positions.shape[0])
    for i in range(len(vertex_count)):
        vertex_count[i] = len(neighbors[i])

    
    bonds, lone, psi_6 = locate_defects(cart_positions) #vertex_count, tri.simplices)
    
    g6, g6_r = compute_g6(cart_positions, psi_6)
    g6 = np.array(g6)
    g6 = g6.real
    
    for (b1, b2) in bonds:
        plot_bond(cart_positions[b1], cart_positions[b2], axs[1])
    
    # Assign colors to points based on their vertex count (default to gray if not in color_dict)
    colors = [color_dict.get(count, 'gray') for count in vertex_count]
    
    # ~ vor_plot = voronoi_plot_2d(vor, ax = axs[1], show_vertices = False, line_width = 2, point_size = 3)
    
    axs[1].triplot(cart_positions[:, 0], cart_positions[:, 1], tri.simplices, color = 'k', zorder = 0)
    
    axs[1].scatter(cart_positions[:, 0], cart_positions[:, 1], facecolors = colors, edgecolors='k', s = 20, zorder = 0)
    
    
    lone_colors = [colors[i] for i in lone]
    axs[1].scatter(cart_positions[lone, 0], cart_positions[lone, 1], facecolors = lone_colors, edgecolors='k', s = 50)
    axs[1].scatter(cart_positions[lone, 0], cart_positions[lone, 1], facecolors = 'none', edgecolors='k', s = 120)
    
    
    #####
    
    rr = np.linspace(min(g6_r), max(g6_r), 300)
    yy = 0.32*rr**(-0.25)
    axs[2].plot(rr, yy)
    
    # ~ print(g6_r)
    
    axs[2].plot(g6_r, np.abs(g6), '-ro')
    axs[2].set_ylim([0.01,0.4])
    # ~ axs[2].set_xscale('log')
    # ~ axs[2].set_yscale('log')
    
    
    # ~ axs[2].imshow(psi_6_interp.real)
    
    # ~ cmap = plt.cm.hsv
    # ~ colors = cmap(np.angle(psi_6)/(np.pi))
    # ~ axs[2].scatter(cart_positions[:, 0], cart_positions[:, 1], facecolors = colors, alpha = np.abs(psi_6), edgecolors='k', s = 20, zorder = 0)
    
    ####
    
    
    
    
    # Add shading to each region based on the number of sides
    # ~ for region_index, region in enumerate(vor.regions):
        # ~ if not -1 in region and len(region) > 0:  # Ignore regions with unbounded edges
            # ~ polygon = [vor.vertices[i] for i in region]
            
            # ~ # The number of sides is simply the number of vertices in the region
            # ~ num_sides = len(region)
            
            # ~ # Use the number of sides for the shading (normalize by the maximum number of sides)
            # ~ color = color_dict.get(num_sides, 'black')
            
            # ~ axs[1].fill(*zip(*polygon), color=color, edgecolor='black')
        
    axs[0].text(0,0,r'$\Gamma$',fontsize=20)
    axs[0].scatter([0],[0],c='red', alpha = 0.5)

    axs[1].set_xlim([min(cart_positions[:, 0]), max(cart_positions[:, 0])])
    axs[1].set_ylim([min(cart_positions[:, 1]), max(cart_positions[:, 1])])
    
    for ax in [axs[0], axs[1]]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_aspect('equal')
    
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.1, wspace=0.1, hspace=0.1)   
        
    if plot == True:
        plt.show()
    else:
        plt.savefig(f'real_space//{fname}.png') #, dpi = 250)



if __name__ == '__main__':
    # ~ par_positions = np.loadtxt(r"C:\Users\kogar\OneDrive\Documents\project_folders\UED_TaS2_liquid\figures\S1\temp_sequence_substrate\T=0.80\par_state_2000.csv", delimiter = ',')
    # ~ plot_reciprocal_lattice(par_positions, q_range = 5)
    
    
    # ~ 
    N1, N2, Psi_6 = count_defects_in_dir(r"C:\Users\kogar\OneDrive\Documents\project_folders\UED_TaS2_liquid\states")

    print(N1, N2, Psi_6)
