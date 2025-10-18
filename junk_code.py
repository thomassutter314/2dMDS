
def get_dislocation_indices_1(vertex_count, simplices):
    
    bonds = [] # list of tuples giving the indices of bonded pars
    bonded = [] # list that keeps track of which pars are already bonded
    lone = [] # list of disclination indices
    five_sites = []
    seven_sites = []
    
    for js in simplices:
        for s in range(3):
            # check whether to make bond
            if (not js[s] in bonded) and (not js[(s+1)%3] in bonded):
                if vertex_count[js[s]] == 5 and vertex_count[js[(s+1)%3]] == 7:
                    bonds.append((js[s], js[(s+1)%3]))
                    bonded.append(js[s])
                    bonded.append(js[(s+1)%3])
                    
                if vertex_count[js[s]] == 7 and vertex_count[js[(s+1)%3]] == 5:
                    bonds.append((js[(s+1)%3], js[s]))
                    bonded.append(js[s])
                    bonded.append(js[(s+1)%3])
                    
    for i, Nv in enumerate(vertex_count):
        if Nv == 5:
            five_sites.append(i)
        if Nv == 7:
            seven_sites.append(i)
            
        if (not i in bonded) and (Nv == 5 or Nv == 7):
            lone.append(i)
    

    
    # ~ print(nx.bipartite.sets(G))
        
    # Compute maximum matching
    # ~ matching = nx.bipartite.maximum_matching(G)
                    
    return bonds, lone
    
    
def get_dislocation_indices(vertex_count, simplices):
    
    bonds = []
    lone = []
    five_sites = []
    seven_sites = []
    
    for js in simplices:
        for s in range(3):
            if vertex_count[js[s]] == 5 and vertex_count[js[(s+1)%3]] == 7:
                bonds.append((js[s], js[(s+1)%3]))

            if vertex_count[js[s]] == 7 and vertex_count[js[(s+1)%3]] == 5:
                bonds.append((js[(s+1)%3], js[s]))
                
    for i, Nv in enumerate(vertex_count):
        if Nv == 5:
            five_sites.append(i)
        if Nv == 7:
            seven_sites.append(i)
            
    # Create a bipartite graph
    G = nx.Graph()
    G.add_nodes_from(five_sites, bipartite=0)  # Set of 5-sites
    G.add_nodes_from(seven_sites, bipartite=1) # Set of 7-sites  
            
    
    # Instead of relying on networkx to infer bipartiteness:
    five_check = set(five_sites)  # Nodes from one partition
    seven_check = set(seven_sites)  # Nodes from the other partition
    
    matching = nx.bipartite.maximum_matching(G, top_nodes=five_check)
    
    print(matching)
    
    return bonds, lone




    
    # Count how many times each vertex appears in the triangulation
    # ~ vertex_count = np.bincount(tri.simplices.flatten(), minlength = cart_positions.shape[0])



# Generates a "diffraction pattern" based on the crystal structure
# Expects an array of shape [# of particles, 2], where the axis=1 determines x/y in a cartesian coordinate system
def plot_reciprocal_lattice(cart_positions, q_range = 3, q_vals = 250, plot = True, fname = 'I_data'):
    
    x0, y0 = np.mean(cart_positions, axis = 0)
    
    # ~ condition = np.sqrt((cart_positions[:, 0] - x0)**2 +  (cart_positions[:, 1] - y0)**2) < 25
    # ~ cart_positions = cart_positions[condition]
    
    # ~ cart_positions = cart_positions[cart_positions[:, 0] > 20]
    # ~ cart_positions = cart_positions[cart_positions[:, 0] < 80]
    
    # ~ cart_positions = cart_positions[cart_positions[:, 1] > 20]
    # ~ cart_positions = cart_positions[cart_positions[:, 1] < 80]
    
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
    
    # interpolate the psi_6 values
    rbf_interpolator = RBFInterpolator(cart_positions, psi_6, kernel='thin_plate_spline')
    
    x_min, x_max = np.min(cart_positions[:, 0]), np.max(cart_positions[:, 0])
    y_min, y_max = np.min(cart_positions[:, 1]), np.max(cart_positions[:, 1])
    num_points = 100
    
    # Evaluate at new (x, y) grid points
    X_new, Y_new = np.meshgrid(np.linspace(x_min, x_max, num_points),
                               np.linspace(y_min, y_max, num_points))
    points_new = np.column_stack([X_new.ravel(), Y_new.ravel()])
    psi_6_interp = rbf_interpolator(points_new).reshape(X_new.shape)
    
    
    g6 = compute_c_2D(psi_6_interp)
    for i in range(2):
        g6 = np.roll(g6, g6.shape[i]//2, axis = i)
        
    # ~ g6 = g6[len(g6)//2:, 0].real
    # ~ g6 /= g6[0]
    # ~ g6_r = np.linspace(0, 0.5*(y_max-y_min), len(g6))
    
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
    axs[2].imshow(g6.real)
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





# Compute spatial correlation function for 2D images
@njit
def compute_c_2D(im):
    im_conj = np.conj(im)
    c = np.zeros(im.shape, dtype = np.complex128)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            for i1 in range(c.shape[0]):
                for j1 in range(c.shape[1]):
                    c[i,j] += im_conj[i1,j1]*im[(i+i1)%c.shape[0],(j+j1)%c.shape[1]]
    return c
 
