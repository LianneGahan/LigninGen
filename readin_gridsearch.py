import matplotlib.pyplot as plt
import json 
import numpy as np
from latest.Code.figure_generation_functions import calculate_avg_bond_distribution, import_library
from latest.Code.compare_simulation_to_exp import calculate_cost

# Convert int64 objects to regular integers
def convert_int64(obj):
    if isinstance(obj, np.int64):
        return int(obj)
    raise TypeError(f'Object of type {obj.__class__.__name__} is not JSON serializable')


def gridsearch_scores(gridsearch_file, biomass_file):
    # Read in min and max kinetic parameters
    #Maximum values of each parameter
    max_vals_kin = np.loadtxt("max_kin_vals.txt")

    #Minimum values of each parameter
    min_vals_kin = np.loadtxt("min_kin_vals.txt")

    #Denotes which parameters should be varied
    variation_bools_kin = np.loadtxt("kin_params_to_randomize.txt",dtype = int)#Contains boolean values for which parameters are varied
    N_kin_vals = variation_bools_kin.size
    #Read in starting parameters if one of the variation bools is 0
    result = all(element == 1 for element in variation_bools_kin)
    if result==0:
        starting_parameters = np.loadtxt("latest/Params/kinetic_parameters.txt")

    # Read in file with simulated conditions
    with open(gridsearch_file, 'r') as json_file:
        gridsearch = json.load(json_file)

    # Read the data from the file
    with open(biomass_file, 'r') as file:
        biomass = json.load(file)

    ###################################################################

    labels=["monomer_addition_rate", "min_monomers", "max_monomers", "sg_ratio", "gg_scaling", "ss_scaling", "sg_scaling"]
    # Create array the size of the number of comparisons to make
    fitness = np.ones((len(gridsearch)))
    # Calculate cost of simulated bond distribution vs. target bond distribution
        # Do this for all simulated distributions
    for count, (key, dict) in enumerate(gridsearch.items()):
        # All check in range bools must switch to 1 in order for the kinetic param set to be used
        check_in_range_bool = np.zeros(N_kin_vals)
        for i in range(N_kin_vals):#
            if variation_bools_kin[i] == 1:
                # Find the parameters which fall inside the appropriate parameter range
                if (dict["kinetic_params"][labels[i]] >= min_vals_kin[i]) and (dict["kinetic_params"][labels[i]] <= max_vals_kin[i]):
                    check_in_range_bool[i] = 1
            # Need to also check that the parameter is even sampled
            elif variation_bools_kin[i] == 0:
                # If the parameter is not sampled, then we need to check if the chosen starting parameter is found in the gridsearch
                if dict["kinetic_params"][labels[i]]==starting_parameters[i]:
                    check_in_range_bool[i] = 1


        # If all "check in range" bools are 1
        result = all(element == 1 for element in check_in_range_bool)
        # Calculate cost between target and simulated bond distributions if all values fall in the range, else inflate fitness
        fitness[count] = calculate_cost(biomass, dict["bonds"]) if result else 99999
        #print("All parameters found in space") if result
    # Choose the lowest euclidian distance
    min = np.where(fitness == np.min(fitness))[0]
    print("Minimum is ", np.min(fitness))

    if np.min(fitness) > 100:
        print("Warning: All starting parameters are outside of the range of minimum and maximum kinetic parameters.")
    # Set best kin specs as the corresponding 
    keystring = "gridsearch_" + str(min[0])
    mon_add_rate = np.log10(gridsearch[keystring]["kinetic_params"]["monomer_addition_rate"])
    min_monomers = gridsearch[keystring]["kinetic_params"]["min_monomers"]
    max_monomers = gridsearch[keystring]["kinetic_params"]["max_monomers"]
    sg_ratio = gridsearch[keystring]["kinetic_params"]["sg_ratio"]
    gg_scaling = gridsearch[keystring]["kinetic_params"]["gg_scaling"]
    ss_scaling = gridsearch[keystring]["kinetic_params"]["ss_scaling"]
    sg_scaling = gridsearch[keystring]["kinetic_params"]["sg_scaling"]
    print("---------------------------------------------")
    #print("Starting parameters for fitting are set as: ")
    print("Monomer addition rate (10^x)", mon_add_rate)
    print("Minimum Monomers: ", min_monomers)
    print("Maximum Monomers: ", max_monomers)
    print("SG Ratio: ", sg_ratio)
    print("GG Scaling Parameter: ", gg_scaling)
    print("SS Scaling Parameter: ", ss_scaling)
    print("SG Scaling Parameter: ", sg_scaling)
    print("---------------------------------------------")


    new_kin_data_out = np.array([mon_add_rate, min_monomers, max_monomers, sg_ratio, gg_scaling, ss_scaling, sg_scaling])
    np.savetxt("best_kin_specs.txt", new_kin_data_out.reshape(1, new_kin_data_out.shape[0]), delimiter = "\t", fmt = "%1.8f")
    np.savetxt("latest/Params/kinetic_parameters.txt",new_kin_data_out.reshape(1, new_kin_data_out.shape[0]), delimiter = "\t", fmt = "%1.8f")
    np.savetxt("smallest_var.txt", np.array([fitness[min[0]]]))

    with open("kin_specs_history.txt", 'a') as f:
        for item in new_kin_data_out:
            f.write("%1.8f\t" % item)
        f.write("\n")

    return(gridsearch, fitness)

#gridsearch_scores("Gridsearch/gridsearch.json", "latest/Params/biomass_data.json")

