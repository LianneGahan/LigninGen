from __future__ import print_function
from common_wrangler.common import (MAIN_SEC, GOOD_RET, INPUT_ERROR, KB, PLANCK_CONST_JS, KCAL_MOL_TO_J_PART,
                                    INVALID_DATA, OUT_DIR, InvalidDataError, warning, process_cfg, make_dir,
                                    create_out_fname, str_to_file, round_sig_figs)
from ligninkmc.kmc_functions import (run_kmc, generate_mol, find_fragments, fragment_size, break_bond_type, gen_tcl)
from ligninkmc.create_lignin import (calc_rates, create_initial_monomers, create_initial_events,
                                     create_initial_state, analyze_adj_matrix,  adj_analysis_to_stdout)
from ligninkmc.kmc_common import (DEF_E_BARRIER_KCAL_MOL, ADJ_MATRIX, MONO_LIST, MONOMER, OX, GROW, C, Monomer, Event)
from ligninkmc.kmc_common import (BO4, B5, BB, B1, B1_ALT, C5O4, AO4, C5C5, G, S, C,
                                   ADJ_MATRIX, BONDS, CHAIN_LEN, RCF_YIELDS)

from scipy.sparse import dok_matrix
from scipy.optimize import curve_fit
import cProfile
import pstats

from rdkit.Chem.Draw import MolToFile
from rdkit.Chem.rdMolInterchange import MolToJSON
from rdkit.Chem import AllChem, Descriptors
from rdkit import Chem

import matplotlib.pyplot as plt
import json
import numpy as np
import copy
#Parallelization
import multiprocessing

# For performance
import time


temp = 298.15  # K
rxn_rates = calc_rates(temp, ea_kcal_mol_dict=DEF_E_BARRIER_KCAL_MOL)

TCL_NAME = "psfgen.tcl"
PSF_NAME = 'lignin'
TOPPAR_DIR = "toppar/"

#fun = par.delayed(run_kmc)


def calculate_bond_distribution(data_dictionary):
    """
    Calculate the normalised distribution of specific bond types within a collection of simulated structures.

    Inputs:
        simulation_result (dictionary): A collection of simulation data containing adjacency matrices and bond information.
        num_sims (int): The number of simulations or structures in the provided data.

    Returns:
        dict: A dictionary containing the normalized distribution of bond types within the simulated structures.
            The keys are bond type names, and the values are the normalized counts of each bond type."""
    #get all bond dictionaries
    Bonds = [inner_dict["Bonds"] for inner_dict in data_dictionary.values() if "Bonds" in inner_dict]
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
        #print(dictionary)
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


def reassign_unique_ids(list_of_list_of_dicts):
    unique_id = 1  # Start with ID 1
    new_dictionary = {}
    for s in list_of_list_of_dicts:

        for num, dict in enumerate(s):
            #print(s[num])
            new_dictionary["ligninkmc_" + str(unique_id)] = s[num]
            unique_id += 1  # Increment unique ID
    return new_dictionary

def count_aromatic_rings(molecule):
    """
    Takes in an RDKIT molecule, and calculates the number of 6 sided rings there are in a structure. This is a solid calculation for how many 
    monomers there are in a given molecule of lignin
    """
        # Generate a molecular graph with aromatic perception
    Chem.Kekulize(molecule)
    sssr = Chem.GetSymmSSSR(molecule)
    aromatic_rings = [ring for ring in sssr if len(ring) == 6 and all(molecule.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring)]
    return len(aromatic_rings)

def split_molecule_with_aromaticity(molecule):
    # Convert the SMILES string into an RDKit molecule object with aromaticity detection
    #molecule = Chem.MolFromSmiles(smiles_string, sanitize=False)

    # Ensure aromaticity is perceived correctly
    Chem.SanitizeMol(molecule, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)

    # Split the molecule into separate molecules (fragments)
    molecules = Chem.GetMolFrags(molecule, asMols=True)
    #print("Number of molecules", len(molecules))
    return molecules


def SMILES_pipeline(adjacency_matrix, monomer_list):
    # Default out is SMILES, which requires getting an rdKit molecule object; also required for everything
    #    except the TCL format
    dict={}
    block = generate_mol(adjacency_matrix, monomer_list)
    molecule_total = Chem.MolFromMolBlock(block)
    molecules = split_molecule_with_aromaticity(molecule_total)
    params = Chem.SmilesWriteParams()
    params.allHsExplicit = True
    params.canonical = True
    params.doIsomericSmiles = True
    """molecules = Chem.MolFromMolBlock(block)
    try:
        smi_str = Chem.MolToSmiles(molecules) + '\n'
    except:
        raise InvalidDataError("Error in producing SMILES string.")
        smi_str = Non
        # if SMI is to be saved, don't output to stdout
    #str_to_file(smi_str, "molecule_smiles.txt", print_info=False)
    #Save SMILES string for individual molecule
    if smi_str is not None:
        smiles = smi_str.split(".")"""

    smiles = [Chem.MolToSmiles(molecule, params) for molecule in molecules]
    number_rings = [count_aromatic_rings(mol) for mol in molecules]

    #smiles = [Chem.MolToSmiles(molecule) for molecule in molecules]
    #Write dictionaryoutput to return to workflow
    for i in range(len(smiles)):
        dict[i] = {
                "smilestring": smiles[i],
                "rings": number_rings[i]}
                #"MW": mw}
    return(dict)

def generate_analysis_parallel(adjacency_matrix, monomer_list, file_writing, savefile_path):
    lignin_dictionaries = []

    # Complete work flow to calculate bond distribution and DP for each lignin molecule
    partial_dict_A = separate_molecule_properties(adjacency_matrix, monomer_list)
    #print("Bond Analysis complete!")

    partial_dict_B = SMILES_pipeline(adjacency_matrix, monomer_list)
    
    smiles_dictionary_id = 0
    bond_dictionary_id = 0

    b1count = 0
    smiles=""
    for i in range(len(partial_dict_A)):
    # Cycling through bond info dictionary, but need to check smiles dict if B-1 bonds are present
        if partial_dict_A[bond_dictionary_id]['Bonds']["b1"] != 0:  
            b1count +=  partial_dict_A[bond_dictionary_id]['Bonds']["b1"] 
            nrings = np.sum([partial_dict_B[smiles_dictionary_id + k]["rings"] for k in range(partial_dict_A[bond_dictionary_id]['Bonds']["b1"]+1)])
            smiles = ""
            for k in range((partial_dict_A[bond_dictionary_id]['Bonds']["b1"]+1)):
                #print(partial_dict_B[smiles_dictionary_id + k]["smilestring"])
                smiles+= (partial_dict_B[smiles_dictionary_id + k]["smilestring"] + ".")
            smiles=smiles[:-1]
            partial_dict_A[bond_dictionary_id]["smilestring"] = smiles


            #smiles = [partial_dict_B[smiles_dictionary_id + k]["smilestring"] + "." for k in range(partial_dict_A[bond_dictionary_id]['Bonds']["b1"])]
            #smiles = ".".join([partial_dict_B[smiles_dictionary_id + k]["smilestring"] for k in range(partial_dict_A[bond_dictionary_id]['Bonds']["b1"])])
            #print(smiles)
            smiles_dictionary_id += partial_dict_A[bond_dictionary_id]['Bonds']["b1"]
            #print("DP: ", partial_dict_A[bond_dictionary_id]['DP'])
            #print("multi estimated DP: ", nrings)
        else:
            smiles = partial_dict_B[smiles_dictionary_id]["smilestring"]
            partial_dict_A[bond_dictionary_id]["smilestring"] = smiles

            #print("DP: ", partial_dict_A[bond_dictionary_id]['DP'])
            #print("estimated DP: ", partial_dict_B[smiles_dictionary_id]["rings"])
        bond_dictionary_id+=1
        smiles_dictionary_id+=1
    #print("Analysis complete!")

    #combined dictionaries
    return(partial_dict_A)

def save_structures(id_list, molecules, dict):
    """
    Saves individual molecule images of type ".png" and ".svg".
    Saves a SMILES string of the molecule
    """

    return()

def separate_molecule_properties(adjacency_matrix, monomer_list):

    bonding_dict = {(4, 8): BO4, (8, 4): BO4, (8, 1): B1, (1, 8): B1, (8, 8): BB, (5, 5): C5C5,
                    (8, 5): B5, (5, 8): B5, (7, 4): AO4, (4, 7): AO4, (5, 4): C5O4, (4, 5): C5O4}
    dict = {}


    #Provides tuples full of monomer indexes, each tuple contains monomers for one lignin molecule. Also provides number of branches in each molecule
    connected, branches = find_fragments(adjacency_matrix)
    # for each molecule
    #Looping over each molecule 
    
    #print(len(connected))
    for number, molecule in enumerate(connected):
        #print("Molecule Analysis complete!")
        molecule = list(molecule)


        #initialise the dictionary  
        bond_count_dict = {BO4: 0,  BB: 0, B5: 0, B1: 0, C5O4: 0, AO4: 0, C5C5: 0}

        #make an empty dok matrix to copy values over to from the global adj matrix
        molecule_adj_matrix = dok_matrix((len(molecule), len(molecule)), dtype=int)
        #make an empty monomer list to copy values over from the global monolist
        #to calculate SG ratio properties
        s_lignol_sum, g_lignol_sum, c_lignol_sum = 0, 0, 0

        for i in range(0, len(molecule)):
            # Confirm the identity of the monomer    
            if monomer_list[molecule[i]].type == 'syringyl':
                s_lignol_sum += 1
            elif monomer_list[molecule[i]].type == 'guaiacol':
                g_lignol_sum += 1
            else: c_lignol_sum += 1

            for j in range(i+1, len(molecule)):
                #print(i, j)
                mol1 = molecule[i]
                mol2 = molecule[j]
                #copy value of adjacency matrix at index mol1 and mol2
                molecule_adj_matrix[i, j] = adjacency_matrix[(mol1, mol2)]
                #Check number at location of ij in adj matrix
                if i != j:
                    #Access bond type by looking at bound carbon location in adjacency matrix
                    bond = adjacency_matrix[(mol1, mol2)], adjacency_matrix[(mol2, mol1)]

                    #so long as the values are not zero, match them according to the bonding dictionary
                    if all(bond):
                        #add one count to the correct dictionary
                        bond_count_dict[bonding_dict[bond]] += 1


        #calculate sg ratio of molecule:
        if g_lignol_sum == 0:
            sg_ratio = s_lignol_sum
        elif s_lignol_sum == 0:
            sg_ratio = 0
        else: sg_ratio = s_lignol_sum/g_lignol_sum

        #Identify number of branches that a particular lignin molecule has
        branching_number = branches[number]

        dict[number] = {
            "Bonds": bond_count_dict,
            "sg_ratio": sg_ratio,
            "DP": len(molecule),
            "Branches": branching_number,
            "Monolignols": {"S": s_lignol_sum,
                            "G": g_lignol_sum,
                            "C": c_lignol_sum},
        }
    return dict

def adjust_energy_barriers(energy_barrier_dict, scale_factor_GG, scale_factor_SS, scale_factor_SG):
    adjusted_dict = copy.deepcopy(energy_barrier_dict)
    # Iterate over each key in the dictionary
    for key, value in adjusted_dict.items():
        # Check if the value is a dictionary
        if isinstance(value, dict):
            # Iterate over each key in the nested dictionary
            for nested_key, nested_value in value.items():
                if nested_key == (G, G):  
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_GG
                if nested_key == (S, G): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SG
                if nested_key == (G, S): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SG
                if nested_key == (S, S): 
                    for subnested_key, subnested_value in nested_value.items():
                        # Assuming subnested_value is an integer or float
                        value[nested_key][subnested_key] *= scale_factor_SS
    return(adjusted_dict)

def generate_analysis(results, num_sims, num_cores, file_writing, savefile_path):
    """
    An analysis of a series of lignin KMC simulations which uses 2 separate pipelines to calculate 
    bond distributions, smiles structures, branching coefficients and chemical functions.

    Input: 
    results - the raw form of simulation output from run_kmc
    num_sums - the number of repeats taken of the simulation scheme 
    """

    #Simulation needs to be split into different adj matrices for different structures 
    #Also need monomer lists from each simulation
    #adjacency matrix for each simulation
    cur_adjs = [results[j][ADJ_MATRIX] for j in range(num_sims)] 
    #list of monomers in each simulation
    monolist = [results[j][MONO_LIST] for j in range(num_sims)]
    #assign an id to each simulation 
    #simulations = ["ligninkmc_" + str(i) for i in range(len(results))]

    # Multiprocessing code is directly placed here
    with multiprocessing.Pool(processes=num_cores) as pool:
        # Lignin libraries is a list of list of dictionaries
        lignin_libraries = pool.starmap(generate_analysis_parallel, [(cur_adjs[i], monolist[i], file_writing, savefile_path) for i, result in enumerate(results)])

    # Let's flatten the dictionary and reassign the IDs
    final_lignin_dictionary = reassign_unique_ids(lignin_libraries)

    #if file_writing:
        #save_structures(simulation_IDs, simulation_molecules, simulation_dicts)

    print("Dictionaries built!") 
    print("Converting and writing to json file...", savefile_path)
    json_data = json.dumps(final_lignin_dictionary)
    with open(savefile_path, "w") as outfile:
        outfile.write(json_data)

    print("Your work here is done!")
    return(final_lignin_dictionary)


def simulation(sg_ratio, ini_num_monos, max_num_monos, mono_add_rate, simulation_reaction_rates, t_max):
    """
    The core functions which run the simulations via lignin KMC

    Input:
    SG_ratio (integer) to describe the ratio of S and G monolignols in a simulation. 
    """
    pct_s = sg_ratio / (1 + sg_ratio)
    # Make choices about what kinds of monomers there are and create them
    # Use a random number and the given sg_ratio to determine the monolignol types to be initially modeled
    np.random.seed()
    monomer_draw = np.random.rand(ini_num_monos)
    initial_monomers = create_initial_monomers(pct_s, monomer_draw)

    # Initialize the monomers, events, and state
    initial_events = create_initial_events(initial_monomers, simulation_reaction_rates)
    initial_state = create_initial_state(initial_events, initial_monomers)
    initial_events.append(Event(GROW, [], rate=mono_add_rate))
    result = run_kmc(simulation_reaction_rates, initial_state, initial_events, n_max = max_num_monos, t_max = t_max, 
                                            sg_ratio = sg_ratio)
    return(result)



