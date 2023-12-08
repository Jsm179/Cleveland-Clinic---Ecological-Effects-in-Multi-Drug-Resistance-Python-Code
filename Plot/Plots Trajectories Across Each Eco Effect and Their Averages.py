import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, lambdify
import os
import pickle as pkl


data_folder = 'Folder where you saved all your simulation data'
output_folder = 'Folder where you want to save output plots'  
    
class Landscape:
    
    #Recreates the drug landscape for each simulation
    def __init__(self, INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB):
        self.INT_AB = INT_AB
        self.INT_BA = INT_BA
        self.H_A = H_A
        self.H_B = H_B
        self.EC50_A = EC50_A
        self.EC50_B = EC50_B
        self.EC50_INT_BA = EC50_INT_BA
        self.EC50_INT_AB = EC50_INT_AB

        x = symbols('x') # C_A
        y = symbols('y') # C_B
        n = symbols('n') # H_A
        m = symbols('m') # H_B
        G = symbols('A') # INT_AB
        H = symbols('B') # INT_BA
        C = symbols('C') # EC50_A
        K = symbols('K') # EC50_B
        D = symbols('D') # EC50_INT_BA
        F = symbols('F') # EC50_INT_AB

        E_A = x**n / ( (C * (1 + (G*y)/(F+y)))**n + x**n)
        E_B = y**m / ( (K * (1 + (H*x)/(D+x)))**m + y**m)
        E_comb = 1-(E_A + E_B - E_A * E_B)
        E_comb = E_comb.subs([
                    (n, self.H_A),
                    (m, self.H_B),
                    (G, self.INT_AB),
                    (H, self.INT_BA),
                    (C, self.EC50_A),
                    (K, self.EC50_B),
                    (D, self.EC50_INT_BA),
                    (F, self.EC50_INT_AB),
                    ])

        self.f = lambdify((x, y), E_comb, 'numpy')


    def plot_combined(self, datas, plot_filename, label):
        fig, ax1 = plt.subplots(figsize=(30, 30))
        
        a = np.linspace(0, 300, 1000)
        b = np.linspace(0, 300, 1000)
        x, y = np.meshgrid(a, b)
        self.k = self.f(x, y)
         
        tot_a1 = []
        tot_b1 = []
        tot_a2 = []
        tot_b2 = []
        tot_a3 = []
        tot_b3 = []

        #Assigns a different color for each type of ecological effect's trajectories
        for data in datas:
            if data[0][0][0] < data[0][2][0]:                   #Negative Game
                color = 'red'
                tot_a1.append(data[0][4].sum(axis=0))
                tot_b1.append(data[0][5].sum(axis=0))
                ax1.plot(tot_a1[-1], tot_b1[-1], color = color)
            elif data[0][0][0] > data[0][2][0]:                 #Postive Game
                color = 'green'
                tot_a2.append(data[0][4].sum(axis=0))
                tot_b2.append(data[0][5].sum(axis=0))
                ax1.plot(tot_a2[-1], tot_b2[-1], color = color)
            else:                                               #No Game                  
                color = 'grey'
                tot_a3.append(data[0][4].sum(axis=0))
                tot_b3.append(data[0][5].sum(axis=0))
                ax1.plot(tot_a3[-1], tot_b3[-1], color = color)

        #Calculates the average across all trajectories for each ecological effect
        avg_a1 = np.mean(tot_a1, axis=0)
        avg_b1 = np.mean(tot_b1, axis=0)
        avg_a2 = np.mean(tot_a2, axis=0)
        avg_b2 = np.mean(tot_b2, axis=0)
        avg_a3 = np.mean(tot_a3, axis=0)
        avg_b3 = np.mean(tot_b3, axis=0)

        #Makes a contour plot with Isogrowth lines at 0.1, 0.25, 0.5, 0.75, 0.9
        #Plots the starting point (same across all simulations plotted on same figure)
        #Marks 
        blue_colors = [f'#{int(255 * i):02X}00FF' for i in np.linspace(0, 1, 6)]
        contour = ax1.contour(a / data[0][13], b / data[0][14], self.k, levels=[0.1, 0.25, 0.5, 0.75, 0.9], colors=blue_colors)
        ax1.clabel(contour, inline=True, fontsize=15, colors='black', fmt='%s')
        ax1.scatter(data[0][0][-1], data[0][1][-1], color='black', s=15)
        ax1.scatter(avg_a1[-1], avg_b1[-1], color = 'black', marker = 'x', s = 400, linewidth = 7.0, zorder = 100000)
        ax1.scatter(avg_a2[-1], avg_b2[-1], color = 'black', marker = 'x', s = 400, linewidth = 7.0, zorder = 100001)
        ax1.scatter(avg_a3[-1], avg_b3[-1], color = 'black', marker = 'x', s = 400, linewidth = 7.0, zorder = 100002)

        # Plots the average trajectories for each eco_effect group
        ax1.plot(avg_a1, avg_b1, color='black', linestyle=':', linewidth = 7.0)
        ax1.plot(avg_a2, avg_b2, color='black', linestyle='--', linewidth = 7.0,)
        ax1.plot(avg_a3, avg_b3, color='black', linestyle='-', linewidth = 7.0)
        


        ax1.set_xlabel('Alpha', fontsize = 30)
        ax1.set_ylabel('Beta', fontsize=30)
        ax1.tick_params(axis='both', labelsize=20)
        #Adjust axes for better visualization
        if self.INT_AB == 0 or self.INT_AB == -0.75:
            ax1.set_xlim(0, 300/ (data[0][13]*2))
            ax1.set_ylim(0, 300/ (data[0][14]*2))
        #ax1.legend()
        plt.savefig(plot_filename)  # Save the plot to the specified file
        plt.close(fig)



def main():
    # Create a dictionary to group files by their prefixes
    file_groups = {}

    #Opens each file and groups all simulation with same drug-drug interaction, starting drug concentration, and mutant distribution together
    for filename in os.listdir(data_folder):
        data_file_path = os.path.join(data_folder, filename)
        if os.path.isfile(data_file_path) and filename.endswith('.pkl'):
            with open(data_file_path, 'rb') as handle:
                parameters, datas = pkl.load(handle)
            parts = filename.split('_')  
            Int_AB = parts[3]
            x_conc = parts[6]
            y_conc = parts[9]
            cov = parts[12][-1]
            prefix = f'Int_AB_{Int_AB}_x_conc_{x_conc}_y_conc_{y_conc}_cov_{cov}'
            if prefix in file_groups:
                file_groups[prefix].append((parameters, datas, filename))
            else:
                file_groups[prefix] = [(parameters, datas, filename)]

    # Iterate through file groups and plot the combined data
    for prefix, file_data in file_groups.items():
        label = prefix # Use the prefix as the label
        plot_filename = os.path.join(output_folder, f'combined_plot_{prefix}.png')
        
        landscape = Landscape(*file_data[0][0].values())  # Use parameters from the first file in the group
        landscape.plot_combined([data[1] for data in file_data], plot_filename, label)

if __name__ == "__main__":
    main()





