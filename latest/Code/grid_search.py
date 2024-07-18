#Parallelization
import joblib as par
from ligninkmc.kmc_functions import (run_kmc, generate_mol, find_fragments, gen_tcl)
from ligninkmc.create_lignin import (calc_rates, create_initial_monomers, create_initial_events,
                                     create_initial_state, analyze_adj_matrix,  adj_analysis_to_stdout)
from ligninkmc.kmc_common import (DEF_E_BARRIER_KCAL_MOL, ADJ_MATRIX, MONO_LIST, MONOMER, OX, GROW, C, Monomer, Event)
from ligninkmc.kmc_common import (BO4, B5, BB, B1, C5O4, AO4, C5C5, G, S, C,
                                   ADJ_MATRIX, BONDS, CHAIN_LEN, RCF_YIELDS)
from compare_simulation_to_exp import calculate_cost
from function_library import simulation, calculate_bond_distribution, generate_analysis, calculate_bond_distribution, adjust_energy_barriers
import numpy as np
import sys
import json
import os
import multiprocessing


def generate_grid_parameters(num_points):
    """
    Generate parameters in a grid based on conditions.

    Parameters:
        kin_mins (list): List of minimum values for each dimension.
        kin_maxs (list): List of maximum values for each dimension.
        num_points (int): Number of points along each dimension.
        variation_bools_kin (list): List of boolean values indicating whether to vary each dimension.

    Returns:
        numpy.ndarray: Array containing the generated parameters.
    """

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt", dtype = int)#Contains boolean values for which parameters are varied
    
    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("grid_search_max.txt")

    #minimum values of each parameter
    min_vals_kin = np.loadtxt("min_kin_vals.txt")
    N_kin_vals = max_vals_kin.size
    kin_mins= np.zeros(N_kin_vals)
    kin_maxs= np.zeros(N_kin_vals)

    for i in range(N_kin_vals):
        if variation_bools_kin[i] == 1:
            kin_mins[i] = (max_vals_kin[i] - min_vals_kin[i])*0.15
            kin_maxs[i] = max_vals_kin[i]*0.85
        else:
            kin_mins[i] = min_vals_kin[i]
            kin_maxs[i] = max_vals_kin[i]
#
    grid_points = []
    for i in range(N_kin_vals):
        if variation_bools_kin[i] != 0:
            grid_points.append(np.linspace(kin_mins[i], kin_maxs[i], num=num_points, endpoint=True))

            
        else:
            grid_points.append(kin_maxs[i])
    # Create grid using meshgrid
    if num_points >1:
        grid = np.meshgrid(*grid_points, indexing='ij')
        parameters = np.vstack([grid[i].reshape(-1) for i in range(N_kin_vals)]).T

    else:
        parameters = [grid_points[i][0] for i in range(len(grid_points))]

    # Reshape grid to (N, M) shape
    print("POINTS ARE: ", parameters)
    print("______________________")
    print("SHAPE IS: ", np.shape(parameters))
    return parameters


#Take in target bond distribution and SG ratio
def simulation_run_analysis_gridsearch(SG_ratio, simulation_time, monomer_addition_rate, minimum_num_monomers, maximum_num_monomers, reaction_rates, number_of_runs, filewriting, filepath):

    with multiprocessing.Pool(processes=int(sys.argv[2])) as pool:
        result = pool.starmap(simulation, [(SG_ratio, minimum_num_monomers, maximum_num_monomers, monomer_addition_rate, reaction_rates, simulation_time) for _ in range(number_of_runs)])
    data_dictionary = generate_analysis(result, number_of_runs, num_cores=int(sys.argv[2]), file_writing=False, savefile_path=filepath)
        

    bond_distribution  =  calculate_bond_distribution(data_dictionary)

    # We should also find a way to check the variance in the distribution, is this something that is important?
    
    return bond_distribution


def grid_search(biomass, N_points):
    """

    Completes a grid search with N_points to find the best starting conditions for the fitting algorithm

    Returns:
    Best_conditions - contains the best monomer addition rate, minimum monomer number and maximum monomer number

    """
    print("Number of points is: ", N_points)
    # Set up initial condition(s) for grid search
    initial_conditions = generate_grid_parameters(N_points)
    print(initial_conditions)
    costs = np.zeros(len(initial_conditions))
    # Find best starting conditions for grid search
    if len(initial_conditions) == 1:
        dir = "Gridsearch/Run_1"
        if not os.path.exists(dir):
            os.makedirs(dir)
            print(f"Directory '{dir}' created successfully.")

    for i in range(len(initial_conditions)):#
        dir = "Gridsearch/Run_" + str(i+1)
        print("CONDITIONS TO TEST ARE: ")
        print(initial_conditions[i])
        if not os.path.exists(dir):
            os.makedirs(dir)
        # Establish initial conditions for grid search point i
        monomer_rate, min_monomers, max_monomers, sg_ratio, scaling_GG, scaling_SS, scaling_SG = initial_conditions[i]
        monomer_rate = 10**monomer_rate
        min_monomers = int(min_monomers)
        max_monomers = int(max_monomers)
        # set output directories for simulation
        dir_path = "Gridsearch/Run_" + str(i+1) + "/Output/library.json"

        temp = 298.15  # K

        # Calculate reaction rates with the new energy barriers
        energy_barriers = adjust_energy_barriers(DEF_E_BARRIER_KCAL_MOL, scaling_GG, scaling_SS, scaling_SG)
        reaction_rates = calc_rates(temp, ea_kcal_mol_dict=energy_barriers)

        simulated_distribution = simulation_run_analysis_gridsearch(sg_ratio, 10000000, monomer_rate, int(min_monomers), int(max_monomers), reaction_rates, 20, filewriting = False, filepath=dir_path)
        costs[i] = calculate_cost(biomass, simulated_distribution)

    # Find the indices where the minimum value occurs
    indices = np.where(costs == np.min(costs))
    print("MINIMUM COST", np.min(costs))
    # Get the first occurrence of the minimum value (you can choose if there are multiple)
    min_index = (indices[0])
    print("MIN INDEX", min_index[0])
    #print("Min cost: ", np.min(costs))
    # Now you can use the indices to index the initial conditions
    best_conditions = initial_conditions[min_index[0]]
    cost = np.min(costs)
    print("BEST CONDITIONS: ", best_conditions)
    print("COST: ", cost)
    return(best_conditions, cost)

N_points = int(sys.argv[1])
N_cores = int(sys.argv[2])
path = "latest/Params/biomass_data.json"

# Open biomass JSON file for reading
with open(path, 'r') as file:
    biomass_data = json.load(file)
best_conditions, cost = grid_search(biomass_data, N_points)
np.savetxt("smallest_var.txt", np.array([cost]))

#print("costs are ", cost)
#print("BEST CONDITIONS", best_conditions)
np.savetxt("latest/Params/kinetic_parameters.txt", best_conditions, delimiter = "\t", fmt = "%1.8f")
np.savetxt("best_kin_specs.txt", best_conditions, delimiter = "\t", fmt = "%1.8f")



