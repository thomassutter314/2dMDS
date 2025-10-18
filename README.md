# 2dMDS
Two dimensional molecular dynamics simulation for reproducing a solid to hexatic to liquid transition, as predicted in KTHNY theory.

For a technical description of the code please review https://doi.org/10.48550/arXiv.2505.04867

To get started run gui.py
Parameters such as lattice size, particle number, initial temperature, final temperature, and cooling rate can be adjusted from the GUI. The GUI will automatically generate a random distribution of particle positions.
To implement a J1 quench, or start the simulation with a preconfigured set of particles positions, adjust the code accordingly on lines 191 and 195 respectively.
As the simulation evolves, particle states are periodically saved in the directory “states”. One can also save a particle state with the “save state” button on the GUI.
