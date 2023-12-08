# Cleveland-Clinic---Ecological-Effects-in-Multi-Drug-Resistance-Python-Code

This Repo is meant to simulate and visualize the emergence of multi-drug resistance accounting for the interaction between drugs, collateral effects, and ecological effects. 
This Repo is broken up into two main portions: a simulating one and a plotting one. In the Simulation folder, the Single Simulation.py file only runs through one simulation and saves the data using in memory storage. 
The other two files run through multiple simulations and writes the result of each one to a file, where you can save them all to the same folder. In the Plot folder, the Plots Trajectories Across Each Eco Effect and 
Their Averages.py file fetches the saved data for plotting. The Animation.py file creates an animation for a single simulation. In general using in memory storage, attach the plotting script as a function and use the 
assigned variable names for plotting. When using the two multiple plotting scripts after using pkl.load(), the specific data you want can be accessed like data[0][i][j]. Where the first index will always be 0, the second
corresponds to which count of the specific simulation you want, and the third is the returned data after the .simulate() function. You can slice further depending on the type of data stored in [j] as well as it's size.
