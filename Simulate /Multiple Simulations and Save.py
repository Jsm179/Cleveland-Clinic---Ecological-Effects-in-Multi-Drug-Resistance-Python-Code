def main():
    #Stores the values that create synergistic, additive, and antagonistic landscapes respectively
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
  
    x_conc_values = [50, 150]  #intial Drug 1 concentractions (between 0 to max(a))
    y_conc_values = [50, 150] #intial Drug 2 concentractions (between 0 to max(b))
    eco_effect_values = [0.5, 1, 1.5] #Ecological interaction parameter between mutants and ancestor
    timesteps = 200  #Number of time intervals to simulate over
    popsize = 100  #Number of mutant populations
    time = 400      #total time to simulate over
    mutantfreq = 0.01  #inital frequency of the all mutant subpopulations
    #Stores the covariance matrices for all four mutant distribution cases in an array (Variance in Beta, Variance in Alpha, Positively Correlated, Negatively Correlated respectively)
    cov = [[[0.01, 0.02], [0.02, 0.4]], [[0.5,0.02],[0.02,0.01]], [[0.3,0.21],[0.21,0.3]], [[0.3,-0.25],[-0.25,0.3]]] 



    folder = '' #Add PATH to folder you want to save the data

  # Loops through all the different drug-drug landscapes, true external drug concentrations, ecological interactions, mutant distributions
  # Does this for each set of starting parameters 100 times
  # 3 landscapes* 3 external drug concentrations* 3 ecological interactions* 4 mutant distributions* 100 times = 10800 datasets
    for INT_AB in INT_AB_values:
        for INT_BA in INT_BA_values:
            if INT_AB == INT_BA: #This is to ensure the drug-drug landscape is symmetric 
               for x_conc in x_conc_values:
                    x_conc = int(x_conc)
                    for y_conc in y_conc_values:
                        y_conc = int(y_conc)
                        if (x_conc, y_conc) in [(50, 150), ( 150, 150), (150, 50)]:
                            for eco_effect in eco_effect_values:
                                for cov_index, cov_matrix in enumerate(cov):
                                    for count in range(100): 
                                        example_DL = DrugLandscape(INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB)
                                        data = example_DL.simulate([1, 1], eco_effect, cov_matrix, timesteps, popsize, time, mutantfreq, x_conc, y_conc)

                                        #Saves each simulation result based off its initial parameters and count as a.pkl file in the specified folder above
                                        filename = f"simulation_INT_AB_{INT_AB}_x_conc_{x_conc}_y_conc_{y_conc}_eco_effect{eco_effect}_cov{cov_index}_count{count}.pkl"  
                                        file_path = os.path.join(folder, filename)

                                        #Saves the drug-drug landscape parameters for plotting purposes later as well as the simulated data
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
