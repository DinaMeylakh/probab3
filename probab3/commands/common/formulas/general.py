import numpy as np
from probab3.commands.common.constants import *
from probab3.commands.common.general_code import LOG_FILE_NAME
import logging

logging.basicConfig(filename=LOG_FILE_NAME, filemode='a', level=logging.INFO, format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

def M(ma, mb, ms):
    return ma + mb + ms

def m(ma, mb, ms):
    return ((ma + mb)*ms)/(ma + mb + ms)

def mB(ma, mb):
    return ma + mb

def CurlyM(ma, mb):
    return (ma*mb)/(ma + mb)

def R(ma, mb, EB, alpha):
    return (alpha*G*ma*mb)/(-2*EB)

def EB(ma, mb, aB):
    return -G*ma*mb/(2*aB)

def qB_from_aBeB(aB, eB):
    return aB*(1-eB)

def eB_from_aBqB(aB, qB):
    return 1 - (qB/aB)

def qB_from_EBeB(ma, mb, EB, eB):
    return (-G*ma*mb*(1-eB))/(2*EB)

def qB_from_EBLB(ma, mb, EB, LB): 
    return ((G*ma*mb)/(2*(-EB)))*(1-(1-((LB**2)*2*(-EB)*mB(ma, mb))/((G**2)*((ma*mb)**3)))**(1/2))

def aB(ma, mb, EB):
    return G*ma*mb/(-2*EB)

def period(ma, mb, EB):
    return (((np.pi**2)*(G**2)*(ma*mb)**3)/(-2*(EB**3)*(ma + mb)))**(1/2)

def LB_from_aBeB(ma, mb, aB, eB):
    return CurlyM(ma,mb)*((G*mB(ma, mb)*aB*(1-eB**2))**(1/2))

def LB_from_EBeB(ma, mb, EB, eB):
    return (((G*CurlyM(ma, mb))**2*ma*mb*mB(ma, mb)*(1-eB**2))/(-2*EB))**(1/2)

def LB_fromEBqB(ma, mb, EB, qB):
    return CurlyM(ma, mb)*((-2*(G*mB(ma, mb)*qB + (EB*(qB**2))/CurlyM(ma,mb)))**(1/2))

def Ls_from_LBL0(LB, CB, L0):
    return (LB**2 - 2*LB*L0*CB + L0**2)**(1/2)

def eB_from_EBLB(ma, mb, EB, LB):
    return (1 + (2*(LB**2)*EB)/((CurlyM(ma, mb)**2)*(G**2)*ma*mb*mB(ma, mb)))**(0.5)

def heaviside(x):
    return 1 * (x > 0)

def Ls_bounds_from_L0max(ma, mb, ms, E0, LB, CB, alpha):
    L0_max_value = L0_max(ma, mb, ms, E0, alpha)
    plus_root = -LB*CB + np.sqrt((L0_max_value**2) - (LB**2)*(1-(CB**2)))
    minus_root = -LB*CB - np.sqrt((L0_max_value**2) - (LB**2)*(1-(CB**2)))
    return minus_root, plus_root

def L0_max(ma, mb, ms, E0, alpha):
    return (((G**2)*ma*mb)/np.abs(E0))**(1/2)*(m(ma, mb, ms)*((M(ma, mb, ms)*alpha)**(1/2)) + CurlyM(ma, mb)*((mB(ma, mb)/2)**(1/2)))

def mass_dep_L0_max(ma, mb, ms, alpha):
    return ((ma*mb))**(1/2)*(m(ma, mb, ms)*((M(ma, mb, ms)*alpha)**(1/2)) + CurlyM(ma, mb)*((mB(ma, mb)/2)**(1/2))) 

def WD_radius(m_WD):
    return 2.9*(10**(6))*((MSun/m_WD)**(1/3))

def star_radius(m_star):
    if m_star < MSun:
        return RSun*((m_star/MSun)**(0.8))
    else:
        return RSun*((m_star/MSun)**(0.55))

def r_TDE(small_mass, small_radius, big_mass):

    return small_radius*((big_mass/small_mass)**(1/3))

def r_schwarzschild(m_BH):
    return (2*G*m_BH)/(c**2)

def CBs_max(ma, mb, ms, E0, LB, EB, alpha):
    L0_max_value = L0_max(ma, mb, ms, E0, alpha)
    CBs_max_L0_LB = - (1 - ((L0_max_value/LB)**2))**(1/2)
    return CBs_max_L0_LB

def calc_Cs_from_CBd_CBs_plus(CBd, CBs):
    return CBd*CBs + np.sqrt((1-CBd**2)*(1-CBs**2))

def calc_Cs_from_CBd_CBs_minus(CBd, CBs):
    return CBd*CBs - np.sqrt((1-CBd**2)*(1-CBs**2))

def calc_CBs_from_CBd_Cs_plus(CBd, Cs):
    return CBd*Cs + np.sqrt((1-CBd**2)*(1-Cs**2))

def calc_CBs_from_CBd_Cs_minus(CBd, Cs):
    return CBd*Cs - np.sqrt((1-CBd**2)*(1-Cs**2))

def calc_C0(L0, LB, Ls, CB, Cs):
    return (LB*CB + Ls*Cs)/L0
    
def choose_C(C_plus, C_minus, C_max=1.0):
    if C_max > 1.0:
        C_max = 1.0
    
    if ((C_plus < C_max) and (-1 < C_plus)) and ((C_minus < C_max) and (-1 < C_minus)):
        return np.random.choice(a=[C_plus, C_minus], p=[0.5, 0.5])
    elif (C_plus < C_max) and (-1 < C_plus):
        return  C_plus
    elif (C_minus < C_max) and (-1 < C_minus):
        return C_minus
    
    raise Exception(f"C options are out of bounds. C_plus={C_plus}, C_minus={C_minus}, C_max={C_max}")
 

def calc_CBd(C0, CB):
    CBd_plus = C0*CB + np.sqrt((1-C0**2)*(1-CB**2))
    CBd_minus = C0*CB - np.sqrt((1-C0**2)*(1-CB**2))
    return choose_C(CBd_plus, CBd_minus)
    
def calc_CBs(CBd, Cs, CBs_max_value=1.0):
    CBs_plus = calc_CBs_from_CBd_Cs_plus(CBd, Cs)
    CBs_minus = calc_CBs_from_CBd_Cs_minus(CBd, Cs)
    return choose_C(CBs_plus, CBs_minus, CBs_max_value)

def calc_Cs_bounds_old(CBd, CBs_max_value=1.0):
    minus_root = calc_Cs_from_CBd_CBs_minus(CBd, CBs_max_value)
    plus_root = calc_Cs_from_CBd_CBs_plus(CBd, CBs_max_value)
    roots = [minus_root, plus_root]
    sorted_roots = sorted(roots)

    CBs_method_plus = calc_CBs_from_CBd_Cs_plus
    CBs_method_minus = calc_CBs_from_CBd_Cs_minus
    critical_Cs_value = CBs_max_value/CBd
    positive_bounds = [[-1.0, min(sorted_roots[0], 1.0)], [max(sorted_roots[1], -1.0), 1.0]]
    negative_bounds = [[max(sorted_roots[0], -1.0), min(sorted_roots[1], 1.0)]] 
    negative_bounds_crit_Cs = [[max(sorted_roots[0], -1.0, critical_Cs_value), min(sorted_roots[1], 1.0)]]
    positive_bounds_crit_Cs_left = [[-1.0, min(sorted_roots[0], 1.0, critical_Cs_value)], 
                                    [max(sorted_roots[1], -1.0), 1.0]]
    positive_bounds_crit_Cs_right = [[-1.0, min(sorted_roots[0], 1.0)], 
                                    [max(sorted_roots[1], -1.0, critical_Cs_value), 1.0]]
    if CBd >= 0.0:
        if critical_Cs_value > 1.0:
            return positive_bounds, [CBs_method_plus, CBs_method_plus]
        if critical_Cs_value < -1.0:
            return negative_bounds, [CBs_method_minus]
        if critical_Cs_value <= sorted_roots[1]:
            if critical_Cs_value >= sorted_roots[0]:
                return positive_bounds + negative_bounds_crit_Cs, [CBs_method_plus, CBs_method_plus, CBs_method_minus]
            return positive_bounds_crit_Cs_left + negative_bounds, [CBs_method_plus, CBs_method_plus, CBs_method_minus]
        return positive_bounds_crit_Cs_right + negative_bounds, [CBs_method_plus, CBs_method_plus, CBs_method_minus]
    
def calc_Cs_bounds(CBd, CBs_max_value=1.0):
    minus_root = calc_Cs_from_CBd_CBs_minus(CBd, CBs_max_value)
    plus_root = calc_Cs_from_CBd_CBs_plus(CBd, CBs_max_value)
    roots = [minus_root, plus_root]
    sorted_roots = sorted(roots)

    CBs_method_plus = calc_CBs_from_CBd_Cs_plus
    CBs_method_minus = calc_CBs_from_CBd_Cs_minus
    critical_Cs_value = CBs_max_value/CBd
    
    CBd_positive_bounds = [[-1.0, min(sorted_roots[0], critical_Cs_value, 1.0)],
                           [max(sorted_roots[1], -1.0), min(1.0, critical_Cs_value)], 
                           [max(sorted_roots[0], -1.0, critical_Cs_value), min(sorted_roots[1], 1.0)]]
    
    CBd_negative_bounds = [[max(-1.0, critical_Cs_value), min(sorted_roots[0], 1.0)],
                           [max(sorted_roots[1], -1.0, critical_Cs_value), 1.0], 
                           [max(sorted_roots[0], -1.0), min(sorted_roots[1], critical_Cs_value, 1.0)]]
    
    CBs_methods = [CBs_method_plus, CBs_method_plus, CBs_method_minus]

    Cs_bounds = CBd_positive_bounds if CBd >= 0.0 else CBd_negative_bounds

    Cs_valid_bounds = []
    CBs_valid_methods = []
    for region, CBs_method in zip(Cs_bounds, CBs_methods) :
        if region[0] < region[1]:
            Cs_valid_bounds.append(region)
            CBs_valid_methods.append(CBs_method)
    
    return Cs_valid_bounds, CBs_valid_methods
 

def vs_min_of_LB(ma, mb, ms, LB, EB, alpha):
    if LB < L0_max(ma, mb, ms, EB, alpha):
        return 0 

    denominator = LB - m(ma, mb, ms)*((2*G*R(ma, mb, EB, alpha)*M(ma, mb, ms))**(1/2))
    vs_min_value = np.sqrt(((-2*EB)/m(ma, mb, ms)) - (2/m(ma, mb, ms))*(((G*mass_dep_L0_max(ma, mb, ms, alpha))/(denominator))**2), dtype=complex)
    return 0 if np.imag(vs_min_value) != 0 else np.real(vs_min_value) 

def ms_verge(ma, mb, EB, LB, alpha):
    char_mass = (LB**2)/(2*G*R(ma, mb, EB, alpha)*(mB(ma, mb)**2))
    return (1/2)*(char_mass + (char_mass**2 + 4*char_mass*mB(ma, mb))**(1/2))

def p_chaotic_ms(ma, mb, ms, EB, LB, alpha, vrms):
    denominator = LB - m(ma, mb, ms)*((2*G*R(ma, mb, EB, alpha)*M(ma, mb, ms))**(1/2))
    delta_vs = ((1/(-2*EB))*((G*mass_dep_L0_max(ma, mb, ms, alpha))/(denominator))**2) 
    char_vs = (((-2*EB)/(m(ma, mb, ms)))**(1/2))*(1 - (1/2)*delta_vs)

    prob_of_vs = delta_vs*((1/np.sqrt(2*np.pi*(vrms**2)))*np.exp((-1/2)*((char_vs**2)/(vrms**2))))

    prob_of_CBs = (1/2) # probability of CBs corresponding to LB > L0_max
    
    return prob_of_vs*prob_of_CBs