from function_library import simulation, calculate_bond_distribution, generate_analysis, calculate_bond_distribution, adjust_energy_barriers
#Parallelization
from ligninkmc.kmc_functions import (run_kmc, generate_mol, find_fragments, gen_tcl)
from ligninkmc.create_lignin import (calc_rates, create_initial_monomers, create_initial_events,
                                     create_initial_state, analyze_adj_matrix,  adj_analysis_to_stdout)
from ligninkmc.kmc_common import (DEF_E_BARRIER_KCAL_MOL, ADJ_MATRIX, MONO_LIST, MONOMER, OX, GROW, C, Monomer, Event)
from ligninkmc.kmc_common import (BO4, B5, BB, B1, C5O4, AO4, C5C5, G, S, C,
                                   ADJ_MATRIX, BONDS, CHAIN_LEN, RCF_YIELDS)
import numpy as np
import sys
import json
import multiprocessing
# Lignin KMC set up code 



#Take in target bond distribution and SG ratio
def simulation_run_analysis(SG_ratio, simulation_time, monomer_addition_rate, minimum_num_monomers, maximum_num_monomers, reaction_rates, number_of_runs, filewriting, filepath):

    with multiprocessing.Pool(processes=int(sys.argv[2])) as pool:
        result = pool.starmap(simulation, [(SG_ratio, minimum_num_monomers, maximum_num_monomers, monomer_addition_rate, reaction_rates, simulation_time) for _ in range(number_of_runs)])
    data_dictionary = generate_analysis(result, number_of_runs, num_cores=int(sys.argv[2]), file_writing=False, savefile_path=filepath)
        

    bond_distribution  =  calculate_bond_distribution(data_dictionary)

    # We should also find a way to check the variance in the distribution, is this something that is important?
    
    return bond_distribution


            



"""#If total monomer number is too low, add more runs to improve statistics
if maximum_num_monomers <=20 and number_of_runs <= 50:
    number_of_runs *= 8
elif maximum_num_monomers <=50 and number_of_runs <= 50:
    number_of_runs *= 4
elif maximum_num_monomers <= 100 and number_of_runs <= 50:
    number_of_runs *= 2"""


params = np.loadtxt("Params/simulation_parameters.txt")
kin_params = np.loadtxt("Params/kinetic_parameters.txt")

biomass_path = "Params/biomass_data.json"
# Open biomass JSON file for reading
with open(biomass_path, 'r') as file:
    biomass_data = json.load(file)


monomer_addition_rate = kin_params[0]
#Scale parameter from search parameter to true value (Logarithmic scaling)
monomer_addition_rate = 10**monomer_addition_rate
minimum_num_monomers = int(kin_params[1])
maximum_num_monomers = int(kin_params[2])
SG_ratio = kin_params[3]

simulation_time = params[0]
filewriting = int(params[1])
number_of_runs = int(sys.argv[1])
scale_factor_GG = kin_params[4]
scale_factor_SS = kin_params[5]
scale_factor_SG = kin_params[6]


temp = 298.15  # K

# Calculate reaction rates with the new energy barriers
energy_barriers = adjust_energy_barriers(DEF_E_BARRIER_KCAL_MOL, scale_factor_GG, scale_factor_SS, scale_factor_SG)
rxn_rates = calc_rates(temp, ea_kcal_mol_dict=energy_barriers)

print("params", monomer_addition_rate, minimum_num_monomers, maximum_num_monomers, SG_ratio, scale_factor_GG, scale_factor_SS, scale_factor_SG)
filepath = "Output/library.json"

#print("sg ratio", SG_ratio)
#print("mon add rate", monomer_addition_rate)
simulation_run_analysis(SG_ratio, simulation_time, monomer_addition_rate, minimum_num_monomers, maximum_num_monomers, rxn_rates, number_of_runs, filewriting, filepath)

