import matplotlib
matplotlib.use('TkAgg')
from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTk, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk)
from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import time
import tkinter as tk
from tkinter import filedialog as fd
import threading

import lattice_gas

# Function to update the Voronoi plot
def update_voronoi(vor, lines, new_points):
    # Recompute the Voronoi diagram with the new points
    vor = Voronoi(new_points)
    
    # ~ # Remove any previous edge lines
    # ~ for line in lines:
        # ~ line.remove()  # Remove old Voronoi edges
    
    # ~ # Clear the list of line objects
    # ~ lines.clear()
    
    # Update Voronoi edges
    for i, ridge in enumerate(vor.ridge_vertices):
        if -1 not in ridge:  # Ignore ridges going to infinity
            line = np.array([vor.vertices[ridge[0]], vor.vertices[ridge[1]]])
            lines[i].set_xdata(line[:, 0])  # Update x data for the line
            lines[i].set_ydata(line[:, 1])  # Update y data for the line

class App():
    def __init__(self):
        self.sim_live = False # GUI starts without the simulation initially running
        self.t1 = time.time() # record the start time
        
        # Some formatting stuff
        self.greenColor = '#50EBFF'
        self.yellowColor = '#FF00CA'
        self.blueColor = '#FCFF00'
        self.redColor = '#03F479'
        
        # Setup a tkinter gui window
        self.root = tk.Tk()
        self.root.geometry("1000x600")
        self.root.configure(bg='white')
        self.root.wm_title("RF GUI")
        # ~ self.root.iconify()
        
        framePady = 18
        
        self.controlsFrame = tk.Frame(master=self.root, width=20, height=20, padx=15, pady = framePady)
        self.controlsFrame.pack(side=tk.LEFT,anchor='nw')
        
        self.displayFrame = tk.Frame(master=self.root, width=20, height=20, padx=15, pady = framePady)
        self.displayFrame.pack(side=tk.RIGHT,anchor='nw')
        
        self.t0 = time.time()
        
        # Make buttons and entries
        self.startButton = tk.Button(
            master=self.controlsFrame,
            text="Start\nSimulation",
            width=10,
            height=2,
            fg='black', bg='white',
            command=self.runSim)
        def resetButtonFunc():
            self.sim_live = False
            self.reset = True
        self.resetButton = tk.Button(
            master=self.controlsFrame,
            text="Reset",
            width=10,
            height=2,
            fg='black', bg='white',
            command=resetButtonFunc)
        self.resetButton["state"] = "disabled"
        self.saveStateButton = tk.Button(
            master=self.controlsFrame,
            text="Save State",
            width=10,
            height=2,
            fg='black', bg='white',
            command=self.saveState)
        self.saveStateButton["state"] = "disabled"
        self.lattice_size_entry = tk.Entry(master=self.controlsFrame)
        self.particle_number_entry = tk.Entry(master=self.controlsFrame)
        self.initial_temp_entry = tk.Entry(master=self.controlsFrame)
        self.final_temp_entry = tk.Entry(master=self.controlsFrame)
        self.cooling_rate_entry = tk.Entry(master=self.controlsFrame)
        self.iterations_per_plot_entry = tk.Entry(master=self.controlsFrame)
        
        # Load the initial values for these entries
        settingsDict = {}
        with open('sim_settings.txt','r') as ssf:
            for line in ssf:
                settingsDict[line[:line.find(":")]] = line[line.find(":")+1:].replace('\n','')
            
        entryList = [self.lattice_size_entry,self.particle_number_entry,self.initial_temp_entry,self.final_temp_entry,self.cooling_rate_entry,self.iterations_per_plot_entry]
        entryCodes = ['lattice_size','particle_number','initial_temp','final_temp','cooling_rate','iterations_per_plot']
        entryType = [int, int, float, float, float, int]
        for key in settingsDict:
            for p in range(len(entryCodes)):
                if entryCodes[p] == key:
                    entryList[p].insert(0,entryType[p](settingsDict[key]))
                    break
        
        pady = 5
        padySep = 10
        sepHeight = 1
        sepWidth = 180
        # Pack buttons and entries
        self.startButton.pack(pady=pady)
        self.resetButton.pack(pady=pady)
        self.saveStateButton.pack(pady=pady)
        tk.Label(master=self.controlsFrame,text="Lattice Size (odd number)").pack(pady=0)
        self.lattice_size_entry.pack(pady=0)
        tk.Label(master=self.controlsFrame,text="Particle Number").pack(pady=0)
        self.particle_number_entry.pack(pady=0)
        tk.Label(master=self.controlsFrame,text="Initial Temperature").pack(pady=0)
        self.initial_temp_entry.pack(pady=0)
        tk.Label(master=self.controlsFrame,text="Final Temperature").pack(pady=0)
        self.final_temp_entry.pack(pady=0)
        tk.Label(master=self.controlsFrame,text="Cooling Rate (1e-6)").pack(pady=0)
        self.cooling_rate_entry.pack(pady=0)
        tk.Frame(master=self.controlsFrame, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        tk.Label(master=self.controlsFrame,text="Iteration Per Plot").pack(pady=0)
        self.iterations_per_plot_entry.pack(pady=0)
        
        
        # Exit sequence in which we rerun the simulation
        self.reset = False
        
        def on_closing():
            if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
                self.reset = False
                
                if self.sim_live:
                    self.sim_live = False
                else:
                    self.root.quit()
        
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
        self.root.destroy()
        
    def runSim(self):
        # Make sure the lattice size is odd
        lattice_size_temp = int(self.lattice_size_entry.get())
        if lattice_size_temp%2 ==0:
            self.lattice_size_entry.delete(0,tk.END)
            self.lattice_size_entry.insert(0,lattice_size_temp + 1)
            
        # Turn off the initial settings buttons
        self.startButton["state"] = "disabled"
        self.resetButton["state"] = "active"
        self.saveStateButton["state"] = "active"
        self.lattice_size_entry["state"] = "disabled"
        self.particle_number_entry["state"] = "disabled"
        self.initial_temp_entry["state"] = "disabled"
        self.cooling_rate_entry["state"] = "disabled"
        self.final_temp_entry["state"] = "disabled"
        
        # Save the initial settings
        settingsDict = {'lattice_size':int(self.lattice_size_entry.get()),
                    'particle_number':int(self.particle_number_entry.get()),
                    'initial_temp':float(self.initial_temp_entry.get()),
                    'final_temp':float(self.final_temp_entry.get()),
                    'cooling_rate':float(self.cooling_rate_entry.get()),
                    'iterations_per_plot':int(self.iterations_per_plot_entry.get()),}               
        with open('sim_settings.txt','w') as ssf:
            for key, value in settingsDict.items():
                ssf.write('%s:%s\n' % (key, value))
        
        # Get user entries
        self.lattice_length = int(self.lattice_size_entry.get())
        self.particle_number = int(self.particle_number_entry.get())
        self.temperature = float(self.initial_temp_entry.get())
        self.target_temperature = float(self.final_temp_entry.get())
        self.cooling_rate = float(self.cooling_rate_entry.get())
        self.coupling = 0.01
        self.target_coupling = 0.01 # 0.3 for paper
        self.coupling_rate = 5e-5
        
        # generate the particle positions
        # ~ self.par_positions = self.lattice_length * np.random.random(size = [self.particle_number, 2]) # 1st index is particle number, 2nd is x or y
        self.par_positions = np.loadtxt('initial_state.csv', delimiter = ',')
        self.par_velocities =  np.zeros(self.par_positions.shape)
        
        self.par_positions_avg = np.zeros(self.par_positions.shape)
        
        self.avg_weight = 0                                
                                             
        self.dm = DisplayManager(lattice_length=self.lattice_length)
        self.dm.plot_intercalants(self.par_positions)
        
        # Attach the matplotlib figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.dm.figure, master=self.root)
        #self.toolbar = NavigationToolbar2Tk( self.canvas, self.root )
        #self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.sim_live = True
        simBatchThr = threading.Thread(target=self.simulationBatch)
        simBatchThr.start() # Process the new data asynchronously
        
    # runs the iteration_per_save number of steps before plotting the results
    def simulationBatch(self):
        iterations_per_plot = int(self.iterations_per_plot_entry.get())
        state_save_index = 0
        # ~ time.sleep(1)       
         
        while self.sim_live:
            
            if state_save_index >= 600 and state_save_index % 100 == 0:
                self.saveState(ID = state_save_index)
            
            t0 = time.time()
            self.par_positions, self.par_velocities = lattice_gas.time_evolve(self.par_positions, self.par_velocities, self.temperature, self.lattice_length, self.coupling, iterations_per_plot)
            tf = time.time()
            
            if self.temperature > self.target_temperature:
                self.temperature += -1e-6 * self.cooling_rate
            else:
                self.temperature = self.target_temperature
            
            if state_save_index > 600:
                if self.coupling < self.target_coupling:
                    self.coupling += self.coupling_rate
                else:
                    self.coupling = self.target_coupling
                    

            state_save_index += 1
            self.dm.update_plot(self.par_positions)
            self.canvas.draw()
            
            
            self.dm.ax.set_title('Temp = %.3f' % self.temperature + ', J = %.3f' % self.coupling + ', CT = %.1f' % (1e6*(tf-t0)/iterations_per_plot) + ' us/iter')

        self.root.quit()
    
    def saveState(self, ID = 0):
        # ~ self.avg_weight += 1
        
        # ~ self.par_positions_avg = self.par_positions/self.avg_weight + (self.avg_weight - 1)/self.avg_weight * self.par_positions_avg
        
        print(f'SAVING STATE, ID = {ID}')
         
        np.savetxt(f'states//par_state_{ID}.csv',self.par_positions,delimiter=",")
        
        lattice_gas.plot_reciprocal_lattice(self.par_positions, plot = False, fname = f'I_data_{ID}')


class DisplayManager():
    def __init__(self, lattice_length):
        self.figure, self.ax = plt.subplots()
        self.ax.set_aspect(1)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.figure.tight_layout()
        
        self.lattice_length = lattice_length
        
        self.ax.set_xlim([0, self.lattice_length])
        self.ax.set_ylim([0, self.lattice_length])
        
    def plot_intercalants(self, par_positions):
        self.particlePlot, = self.ax.plot(par_positions[:, 0], par_positions[:, 1], 'o', color = 'black')
        
        # ~ # Compute initial Voronoi diagram
        # ~ self.vor = Voronoi(par_positions)
        
        # ~ # Plot Voronoi edges
        # ~ self.lines = []
        # ~ for ridge in self.vor.ridge_vertices:
            # ~ if -1 not in ridge:  # Ignore ridges going to infinity
                # ~ line = np.array([self.vor.vertices[ridge[0]], self.vor.vertices[ridge[1]]])
                # ~ line_plot, = self.ax.plot(line[:, 0], line[:, 1], color='blue')  # Plot edges
                # ~ self.lines.append(line_plot)
        
    def update_plot(self, par_positions):
        self.particlePlot.set_xdata(par_positions[:, 0])
        self.particlePlot.set_ydata(par_positions[:, 1])
        
        # ~ update_voronoi(self.vor, self.lines, par_positions)
        
        
if __name__ == '__main__':
    app = App()
    while app.reset:
        app = App()



    

    
