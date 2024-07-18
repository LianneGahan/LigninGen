import numpy as np

def calculate_cost(biomass, simulated_distribution):
    """ Calculates the Euclidian distance between simulated and desired distribution 
    
    Input: 

    Biomass - a JSON file containing a dictionary with information about the biomass. This will
                contain a bond distribution, sg ratio, average molecular weight and the level of
                detail in the biomass data. i.e. Complete, Klose or Terrell type

    simulated distribution - A python dictionary containing propotions of bo4, bb, b5, b1, ao4, 4o5 and 55
                bond types


    Output:

    Cost - The euclidian distance between the simulated and target bond distribution
    """


    if biomass["bond_detail_level"] == "klose":
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"]])
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)
        
    elif biomass["bond_detail_level"] == "Complete":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"], 
                                        biomass["ao4"],
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                   simulated_distribution["5o4"],
                                    simulated_distribution["ao4"], 
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)

    elif biomass["bond_detail_level"] == "terrell":
        # Calculate the difference between simulated and desired distribution
        _target_distribution = np.array([biomass["bo4"], 
                                        biomass["bb"],
                                        biomass["b5"],
                                        biomass["b1"],
                                        biomass["4o5"] + 
                                        biomass["ao4"] +
                                        biomass["55"]])
        
        _simulated_distribution = np.array([simulated_distribution["bo4"], 
                                    simulated_distribution["bb"],
                                    simulated_distribution["b5"],
                                    simulated_distribution["b1"],
                                    simulated_distribution["5o4"] +
                                    simulated_distribution["ao4"] +
                                    simulated_distribution["55"]])
        
        # Normalise so we keep the ratios of bond types
        _target_distribution = _target_distribution/np.sum(_target_distribution)
        _simulated_distribution = _simulated_distribution/np.sum(_simulated_distribution)



    cost = np.linalg.norm(_target_distribution - _simulated_distribution)
    #print("Euclidian distance is: ", cost)
    return(cost)

