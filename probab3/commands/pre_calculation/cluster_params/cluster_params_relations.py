"""This module contains functions that relate the cluster parameters to each other."""
from probab3.commands.common.constants import *
import numpy as np

SLOPE = 1.86276 
M_START = 10**(3.42374) 

def select_rh_for_Mc(Mc):
    """Return rh in parsecs"""
    G_in_km_MSun = G *((MSun)/((1000)**3))
    rh_in_km = (6*G_in_km_MSun*(((np.pi**2)/8) -1))*(M_START**(2/(SLOPE)))/(np.pi*(Mc**((2/(SLOPE)) - 1)))
    return rh_in_km*(1000)/parsec

def select_Mc_for_sigma(sigma):
    return M_START*(sigma)**(SLOPE)