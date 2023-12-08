import numpy as np
from sympy import symbols, lambdify
import copy

def main():
    # Initializes the starting landscape parameters. Edit out the two you are not using.
    """
    # Antagonistic
    INT_AB = 1.5
    INT_BA = 1.5
    H_A = 2
    H_B = 2
    EC50_A = 50
    EC50_B = 50
    EC50_INT_BA= 25
    EC50_INT_AB= 25
    """

    """
    # Synergistic
    INT_AB = -0.75
    INT_BA = -0.75
    H_A = 2
    H_B = 2
    EC50_A = 50
    EC50_B = 50
    EC50_INT_BA= 25
    EC50_INT_AB= 25
    """
    
    
    # Additive
    INT_AB = 0
    INT_BA = 0
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
    
    x_conc = 150  #intial Drug 1 concentraction (between 0 to max(a))
    y_conc = 150 #intial Drug 2 concentraction (between 0 to max(b))
    eco_effect1 = 1  #Ecological interaction parameter between mutants and ancestor ( =1 no interaction, <1 positive interaction, >1 negative interaction)
    timesteps = 400  #Number of time intervals to simulate over
    popsize = 100  #Number of mutant populations
    time = 200      #total time to simulate over
    mutantfreq = 0.01  #inital frequency of mutants

    # Variance in beta
    cov = [[0.01, 0.02], [0.02, 0.4]]

    # Variance in alpha
    cov2 = [[0.4,0.02],[0.02,0.01]]

    # Variance is correlated
    cov3 = [[0.3,0.21],[0.21,0.3]]

    # Variance is anti-correlated
    cov4 = [[0.3,-0.25],[-0.25,0.3]]

    #Create an instance of the class to simulate results based off above starting parameters
    example_DL = DrugLandscape(INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB)
    data = example_DL.simulate([1,1], eco_effect3, cov, timesteps, popsize, time, mutantfreq, x_conc, y_conc)
    datas = []
    datas.extend([data])

class DrugLandscape:
    def __init__(self, INT_AB, INT_BA, H_A, H_B, EC50_A, EC50_B, EC50_INT_BA, EC50_INT_AB):
        
        #Calculates the drug-drug landscape using two dose-response Hill, under the Bliss independence model 
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
        
    #Gets the growth value at a specific point
    def get_growth(self, x_coor, y_coor):
        return self.f(x_coor, y_coor)
        
    
    def simulate(self, mean, eco_effect, cov, timesteps, popsize, time, mutantfreq, x_conc, y_conc):
        # Generate mutant subpopulation from a multivariate normal distribution
        # with mean = [1,1] and covariance matrix 'cov'
        #Checks to see if any are negative and sets them to 0
        x_D, y_D = np.random.multivariate_normal(mean, cov, popsize-1).T
        for i in range(popsize-1):
            if x_D[i] < 0:
                x_D[i] = 0
            if y_D[i] < 0:
                y_D[i] = 0
        
        #Rescales the mutants position on drug-drug landscape based of the type of ecological interaction for our starting time at ~100% ancestor frequency
        x_C = eco_effect * x_D
        y_C = eco_effect * y_D

        # Adds the ancestor (point (1,1)) into our arrays at the 0-th indexed position
        x_D = np.insert(x_D, 0, 1)
        y_D = np.insert(y_D, 0, 1)
        x_C = np.insert(x_C, 0, 1)
        y_C = np.insert(y_C, 0, 1)

        x_D = np.array(x_D)
        x_C =  np.array(x_C)
        y_D = np.array(y_D)
        y_C = np.array(y_C)
        
        #Creates an array of the growth values for our mutants at ~0% and ~100% ancestor frequency
        g_D = self.f(x_D*x_conc, y_D*y_conc)
        g_C = self.f(x_C*x_conc, y_C*y_conc)

        #Creates arrays to store simulated values
        #Rows correspond to mutant subpopulations
        #Columns correspond to the number of time increments
        t = np.linspace(0,time,timesteps)                    #Not used in simulation, stored later on for plotting
        dt = (t[1]- t[0])/2                                  #Can adjust differentiation parameter as necessary: used to calculate time-derivative using np.gradient below
        alpha = np.zeros((popsize, timesteps))
        beta = np.zeros((popsize, timesteps))
        pop = np.zeros((popsize, timesteps))
        freq = np.zeros((popsize, timesteps))
        growth1 = np.zeros((popsize, timesteps))
        growth2 = np.zeros((popsize, timesteps))
        der_x = np.zeros((popsize, timesteps))
        der_y = np.zeros((popsize, timesteps))
        d_x = np.zeros((popsize, timesteps))
        d_y = np.zeros((popsize, timesteps))
        eff_x = np.zeros((popsize, timesteps))
        eff_y = np.zeros((popsize, timesteps))
        fracsens = np.zeros_like(t)
        
        #Initializes the starting time values
        ancestor_freq = 1- mutantfreq
        n0 = np.full(popsize-1, mutantfreq/popsize)
        n0= np.append(n0, ancestor_freq)
        pop[:,0] = copy.deepcopy(n0)
        fracsens[0] = pop[popsize-1,0] / np.sum(pop[:,0])
        eff_x[:,0] = x_D[:]*(1 - fracsens[0]) + x_C[:]*fracsens[0]
        eff_y[:,0] = y_D[:]*(1 - fracsens[0]) + y_C[:]*fracsens[0]
        

        #Calculates the next number of cells per populations based off previous ancestor frequency and each populations position on the drug-drug landscape
        for j in range(timesteps-1):
            for i in range(popsize):
                pop[i,j+1] = pop[i,j]*(np.exp(self.f(eff_x[i,j]*x_conc, eff_y[i,j]*y_conc)))
            
            fracsens[j+1] = pop[popsize-1,j+1]/(np.sum(pop[:,j+1]))
            eff_x[:,j+1] = x_D[:]*(1 - fracsens[j+1]) + x_C[:]*fracsens[j+1]
            eff_y[:,j+1] = y_D[:]*(1 - fracsens[j+1]) + y_C[:]*fracsens[j+1]

        #Calculates the frequency of each population at each timestep 
        #Uses that the calculate the corresponding alpha, beta, and growth weighted by the populations frequency at each timestep
        #Also calculates the unweighted growth values for each population over time
        for i in range(popsize):
            for j in range(timesteps):
                freq[i,j] = pop[i,j]/np.sum(pop[:,j])
                alpha[i,j] = freq[i,j]*eff_x[i,j]
                beta[i,j] = freq[i,j]*eff_y[i,j]
                growth1[i,j] = freq[i,j]*self.f(eff_x[i,j]*x_conc, eff_y[i,j]*y_conc)
                growth2[i,j] = self.f(eff_x[i,j]*x_conc, eff_y[i,j]*y_conc)

        
        #Calculates the time-derivative of each mutant's trajectory on the landscape
        #Ancestor set to 0 (position on landscape does not change)
        for i in range(popsize-1):
            for j in range(timesteps):
                der_x[i+1,j] = x_D[i+1]*(eco_effect-1)*fracsens[j]*(growth2[0,j]- growth1.sum(axis=0)[j])
                der_y[i+1,j] = y_D[i+1]*(eco_effect-1)*fracsens[j]*(growth2[0,j]- growth1.sum(axis=0)[j])
        
        #Calculates the time-derivative of each mutant's trajectory weighted by its frequency
        for i in range(popsize):
            for j in range(timesteps):
                d_x[i,j] = der_x[i,j]*freq[i,j]+eff_x[i,j]*freq[i,j]*(growth2[i,j]-growth1.sum(axis=0)[j])
                d_y[i,j] = der_y[i,j]*freq[i,j]+eff_y[i,j]*freq[i,j]*(growth2[i,j]-growth1.sum(axis=0)[j])
                
        #Gets the time-derivative of the mean alpha and beta 
        d_x = d_x.sum(axis=0)
        d_y = d_y.sum(axis=0)
        
        
        #Calulates the integral of the mean alpha and beta to recover the mean trajectory (integration constant is always 1)
        int_der_x = np.ones(timesteps)  
        int_der_y = np.ones(timesteps)
        
        for j in range(1, timesteps):
            int_der_x[j] = int_der_x[j-1] + d_x[j]
            int_der_y[j] = int_der_y[j-1] + d_y[j]
        

        #Calculates the gradients the mean alpha and beta based on simulated values
        deri_x = np.gradient(alpha.sum(axis=0), dt)
        deri_y = np.gradient(beta.sum(axis=0), dt)

        data = []
        
        data.extend([x_D, y_D, x_C, y_C, alpha, beta, growth1, popsize, fracsens, timesteps, time, t, alpha, beta, freq, x_conc, y_conc, int_der_x, int_der_y, d_x, deri_x, d_y, deri_y])

        return data
       

if __name__ == '__main__':
    main()   
