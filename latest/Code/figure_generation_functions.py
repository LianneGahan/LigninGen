import numpy as np
import json
from scipy.optimize import curve_fit


def import_library(json_file_name, string):
    """
    Input: File location of .json file containing the outputs of a Lignin KMC simulation
    Output: 
        _ids: The list of id keys associated to each dictionary (one lignin molecule)
        _dict: The library containing multiple dictionaries, each of which corresponds to the properties of one lignin molecule.
        Each dict contains:

        SMILES string of the molecule
        Bonds :"bo4" "bb" "b5" "b1" "a04" "55" "4o5" contributions
        func_grps
        MW
        DP

    """
# Create a defaultdict with a list as the default value
    _dict = []
    _ids = []
    # Load the JSON data
    with open(json_file_name, 'r') as json_file:
        json_data = json.load(json_file)

        for i in range(1, 1000000):
            id_string = string + str(i)
            if id_string in json_data.keys():
                _dict.append(json_data[id_string])
                _ids.append(id_string)
            else: break
    return(_ids, _dict)

def calculate_avg_bond_distribution(library):
    """
    Calculate the normalised distribution of specific bond types within a collection of simulated structures.

    Inputs:
        simulation_result (list): A collection of simulation data containing adjacency matrices and bond information.
        num_sims (int): The number of simulations or structures in the provided data.

    Returns:
        dict: A dictionary containing the normalized distribution of bond types within the simulated structures.
            The keys are bond type names, and the values are the normalized counts of each bond type."""
    #get all bond dictionaries
    Bonds = [inner_dict["Bonds"] for inner_dict in library if "Bonds" in inner_dict]
    #print(Bonds)
    BO4 = []
    BB = []
    B5 = []
    _55 = []
    _4O5 = []
    B1 = []
    AO4 =[]
    
    # get all bond types
    for dictionary in Bonds:
        if 'bo4' in dictionary:
            BO4.append(dictionary['bo4'])
        if 'bb' in dictionary:
            BB.append(dictionary['bb'])
        if 'b5' in dictionary:
            B5.append(dictionary['b5'])
        if '55' in dictionary:
            _55.append(dictionary['55'])
        if '5o4' in dictionary:
            _4O5.append(dictionary['5o4'])
        if 'b1' in dictionary:
            B1.append(dictionary['b1'])
        if 'ao4' in dictionary:
            AO4.append(dictionary['ao4'])
    #Normalise the distribution:
    norm = np.sum(BO4) + np.sum(BB) + np.sum(B5) + np.sum(_55) + np.sum(AO4) + np.sum(_4O5) + np.sum(B1)
    #print("norm ", norm)
    dict = {
            "bo4": np.sum(BO4)/norm,
            "bb": np.sum(BB)/norm,
            "b5": np.sum(B5)/norm,
            "ao4": np.sum(AO4)/norm,
            "5o4": np.sum(_4O5)/norm,
            "b1": np.sum(B1)/norm, 
            "55": np.sum(_55)/norm
            }
    return(dict)

def calculate_bond_distribution_fitness(simulated_distribution, target_distribution):
    """
    Inputs:
    simulated_distribution - a dictionary which contains a value for each of the dominant bonds from a generated library of lignin. 
                            This does not have to be normalised but it is better if it is. 

    target_distribution - a dictionary which contains a value for each of the dominant bonds. This does not have to be 
                            normalised but it is better if it is.
    """
    _target_distribution = np.array([target_distribution["bo4"], 
                                    target_distribution["bb"],
                                    target_distribution["b5"],
                                    target_distribution["b1"],
                                    target_distribution["4o5"]+
                                    target_distribution["ao4"]+
                                    target_distribution["55"]])
    _target_distribution = _target_distribution/np.sum(_target_distribution)
    _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                    simulated_distribution["5o4"]+
                                    simulated_distribution["ao4"]+
                                    simulated_distribution["55"]])
    _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)
    
    # Calculate the Euclidean distance
    euclidean_distance = np.linalg.norm(_target_distribution - _simulated_distribution)
    return(euclidean_distance)

def evaluate_minimum_maximum_DP(data):
    """
    Input: data -  A list containing multiple dictionaries. Each dictionary contains information about 1 lignin molecule. 

    Output: The minimum DP observed in a dataset and the maximum DP observed in a dataset
    """

    # Initialize lists to store the values
    values = []

    # Iterate through the list of dictionaries
    for record_dict in data:
        # Check if the 'value' key exists in the dictionary
        if 'DP' in record_dict:
            # Append the 'value' to the values list
            values.append(record_dict['DP'])
    # Find the minimum and maximum values using the min() and max() functions
    if values:
        min_value = min(values)
        max_value = max(values)
        #print("Minimum value:", min_value)
        #print("Maximum value:", max_value)
        return(min_value, max_value)
    else:
        print("No values found for the specified key.")

def exp_func(x, A, B):
    return B * np.exp(-A * x)

def exponential_fit(x, y):
    # Provide an initial guess for parameter A
    initial_guess_A = 0.4
    initial_guess_B = 4
    # Use curve fit with the initial guess for the fitting parameter
    popt, _ = curve_fit(exp_func, x, y, p0=(initial_guess_A, initial_guess_B))
    A = popt[0]
    B = popt[1]
    return A, B, popt

def calculate_mean_DP(library):
    """
    Input: The lignin library of one simulated biomass. This will contain structural information about each individual lignin molecule

    Output: Mean value of size distribution
    """
    # Initialize lists to store the values
    values = []
    
    # Iterate through the list of dictionaries
    for _dict in library:
        # Append the 'value' to the values list
        values.append(_dict["DP"])
        #print(values[-1])


    mean_DP = np.mean(values)
    return mean_DP
    

def readin_fitting_results(file_path):

    """
    Input: file_path - The file location of an output file containing each iteration of the fitting parameters in a fitting distribution

    euclidian distance, monomer addition rate, maximum monomer number and minimum monomer number

    Each row contains 1 instance of this. This function reads the final line of the file, indicating this is the end point of the fitting.
    """

    # Open the file in read mode
    with open(file_path, 'r') as file:
        # Read all lines into a list
        lines = file.readlines()

    # Check if the file is not empty
    if lines:
        # Get the last line
        last_line = lines[-1].strip()  # Remove any leading/trailing whitespace

        # Split the last line into four numbers using a delimiter (e.g., space)
        numbers = last_line.split()

        # Check if there are exactly four numbers
        if len(numbers) == 4:
            # Convert the strings to integers or floats as needed
            try:
                euclidian_distance = float(numbers[0])
                monomer_addition_rate = float(numbers[1])
                minimum_monomers = float(numbers[2])
                maximum_monomers = float(numbers[3])

                # Now you have your four numbers
                data = [euclidian_distance, monomer_addition_rate, minimum_monomers, maximum_monomers]
                #print("Data:", data) 
            except ValueError:
                print("Error: Could not convert all elements to numbers.")
        else:
            print("Error: The last line does not contain exactly four numbers.")
        return(data)
    else:
        print("Error: The file is empty.")
    

