import numpy as np
from sympy import symbols, lambdify
import copy
import os
import pickle as pkl

#Same concept at the Multiple Simulations and Plot file, just copy this into the Single Simulation file
#This time it makes the starting drug concentrations along the 0.1 isogrowth lines for each respective drug landscape type
def main():
    INT_AB_values = [-0.75, 0, 1.5]
    INT_BA_values = [-0.75, 0, 1.5]
    H_A = 2
    H_B = 2
    EC50_A = 50
    EC50_B = 50
    EC50_INT_BA= 25
    EC50_INT_AB= 25

    # Set the extent and resolution you want your drug-drug lanscape
    a = np.linspace(0, 300, 1000) 
    b = np.linspace(0, 300, 1000)
    x, y = np.meshgrid(a, b)
    
    antag_x, antag_y = [55,170, 273], [55,170, 273]  #intial Drug 1 concentraction (between 0 to 200)
    add_x, add_y = [27, 74, 130], [27, 74, 130]  #intial Drug 2 concentraction (between 0 to 200)
    syn_x, syn_y = [17, 40, 79], [17, 40, 79]
    eco_effect_values = [0.5, 1, 1.5]
    timesteps = 200  #Number of time intervals to simulate over
    popsize = 100  #Number of mutant populations
    time = 400      #total time to simulate over
    mutantfreq = 0.01  #inital frequency of mutants

    cov = [[[0.01, 0.02], [0.02, 0.4]], [[0.5,0.02],[0.02,0.01]], [[0.3,0.21],[0.21,0.3]], [[0.3,-0.25],[-0.25,0.3]]]



    folder = 'Path to output folder'
    
    for INT_AB in INT_AB_values:
        for INT_BA in INT_BA_values:
            if INT_AB == INT_BA:
               if INT_AB == -0.75:
                   for x_conc in syn_x:
                       x_conc = int(x_conc)
                       for y_conc in syn_x:
                           y_conc = int(y_conc)
                           if (x_conc, y_conc) in [(17,79), (40,40), (79,17)]:
                                for eco_effect in eco_effect_values:
                                        for cov_index, cov_matrix in enumerate(cov):
                                            for count in range(100):
                                                example_DL = DrugLandscape(INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB)
                                                data = example_DL.simulate([1, 1], eco_effect, cov_matrix, timesteps, popsize, time, mutantfreq, x_conc, y_conc)

                                                filename = f"simulation_INT_AB_{INT_AB}_x_conc_{x_conc}_y_conc_{y_conc}_eco_effect{eco_effect}_cov{cov_index}_count{count}.pkl"
                                                file_path = os.path.join(folder, filename)

                                                parameters = {
                                                    'INT_AB': INT_AB,
                                                    'INT_BA': INT_BA,
                                                    'H_A': H_A,
                                                    'H_B': H_B,
                                                    'EC50_A': EC50_A,
                                                    'EC50_B': EC50_B,
                                                    'EC50_INT_BA': EC50_INT_BA,
                                                    'EC50_INT_AB': EC50_INT_AB,
                                                }
                                                with open(file_path, 'wb') as handle:
                                                    pkl.dump((parameters,[data]), handle)
  
            


    for INT_AB in INT_AB_values:
            for INT_BA in INT_BA_values:
                if INT_AB == INT_BA:
                    if INT_AB == 0:
                        for x_conc in add_x:
                            x_conc = int(x_conc)
                            for y_conc in add_x:
                                y_conc = int(y_conc)
                                if (x_conc, y_conc) in [(27,130), (74,74), (130,27)]:
                                    for eco_effect in eco_effect_values:
                                            for cov_index, cov_matrix in enumerate(cov):
                                                for cov_index, cov_matrix in enumerate(cov):
                                                    for count in range(100):
                                                        example_DL = DrugLandscape(INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB)
                                                        data = example_DL.simulate([1, 1], eco_effect, cov_matrix, timesteps, popsize, time, mutantfreq, x_conc, y_conc)

                                                        filename = f"simulation_INT_AB_{INT_AB}_x_conc_{x_conc}_y_conc_{y_conc}_eco_effect{eco_effect}_cov{cov_index}_count{count}.pkl"
                                                        file_path = os.path.join(folder, filename)

                                                        parameters = {
                                                            'INT_AB': INT_AB,
                                                            'INT_BA': INT_BA,
                                                            'H_A': H_A,
                                                            'H_B': H_B,
                                                            'EC50_A': EC50_A,
                                                            'EC50_B': EC50_B,
                                                            'EC50_INT_BA': EC50_INT_BA,
                                                            'EC50_INT_AB': EC50_INT_AB,
                                                        }
                                                        with open(file_path, 'wb') as handle:
                                                            pkl.dump((parameters,[data]), handle)
                                  
    for INT_AB in INT_AB_values:
        for INT_BA in INT_BA_values:
            if INT_AB == INT_BA:
                if INT_AB == 1.5:
                    for x_conc in antag_x:
                        x_conc = int(x_conc)
                        for y_conc in antag_x:
                            y_conc = int(y_conc)
                            if (x_conc, y_conc) in [(55,273), (170,170), (273,55)]:
                                for eco_effect in eco_effect_values:
                                        for cov_index, cov_matrix in enumerate(cov):
                                            for count in range(100):
                                                example_DL = DrugLandscape(INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB)
                                                data = example_DL.simulate([1, 1], eco_effect, cov_matrix, timesteps, popsize, time, mutantfreq, x_conc, y_conc)

                                                filename = f"simulation_INT_AB_{INT_AB}_x_conc_{x_conc}_y_conc_{y_conc}_eco_effect{eco_effect}_cov{cov_index}_count{count}.pkl"
                                                file_path = os.path.join(folder, filename)

                                                parameters = {
                                                    'INT_AB': INT_AB,
                                                    'INT_BA': INT_BA,
                                                    'H_A': H_A,
                                                    'H_B': H_B,
                                                    'EC50_A': EC50_A,
                                                    'EC50_B': EC50_B,
                                                    'EC50_INT_BA': EC50_INT_BA,
                                                    'EC50_INT_AB': EC50_INT_AB,
                                                }
                                                with open(file_path, 'wb') as handle:
                                                    pkl.dump((parameters,[data]), handle)
