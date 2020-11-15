import numpy as np
import numpy.random as npr

def sample_transmission_matrices(ND, NW, td, tw, pair_mixing=0):
    #drivers are paired 1-2, 3-4, etc. Have interaction (1-pair_mixing)*td
    #certain on site staff have regular contact with a subset of drivers (i.e. sign in and out)
        # transmission rate tS ~ small, brief interaction
    #all other on site staff have some background rates proportional to:
        #number of onsite staff infected * tW / (NW-1)
        #number of drivers infected  * tWD / ND
    #drivers have bg rates proportional to:
        #number of drivers infected * (pair_mixing + normal_mixing)*td/(ND-1) 
        #number of on site staff infected * tWD / NW 
    #could also consider fixed warehouse teams / pairs
    
    