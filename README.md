# 2dMDS
Two dimensional molecular dynamics simulation for reproducing a solid to hexatic to liquid transition, as predicted in KTHNY theory.

For a technical description of the code please review https://doi.org/10.48550/arXiv.2505.04867

To get started run gui.py
Parameters such as lattice size, particle number, initial temperature, final temperature, and cooling rate can be adjusted from the GUI. The GUI will automatically generate a random set of particle positions.

To initiate the simulation with a previous particle position set, change the load_par_positions argument in App to True. This will load the initial_state.csv file in the main directory (this .csv file can be replaced with any set of particle positions).

By default, the J1 coupling is set by coupling = 0.3. A J1 quench can be implemented by changing target_coupling and coupling_rate. The former sets the J1 coupling to which the system will quench while the latter sets the J1 step per simulation step.
