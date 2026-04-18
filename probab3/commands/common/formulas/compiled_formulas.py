from probab3.commands.common.constants import G, MSun
import numpy as np
from numba import jit, njit

@njit
def hyperbollic_dsigma_dEBdLBdCB(EBm, LBm, CB, mam, mbm, msm, alpha, E0m, L0m):
    EB = EBm*MSun
    LB = LBm*MSun
    ma = mam*MSun
    mb = mbm*MSun
    ms = msm*MSun
    E0 = E0m*MSun
    L0 = L0m*MSun

    Ls = np.sqrt(LB**2 - 2*LB*L0*CB + L0**2)
    current_R = (alpha*G*ma*mb)/(-2*EB)
    m = ((ma + mb)*ms)/(ma + mb + ms)
    mB = ma + mb
    M = mB + ms

    heaviside_term =(2*current_R*G*M*(m**2)
                    +2*m*(E0-EB)*(current_R**2)-Ls**2)
        
    if heaviside_term <= 0:
        return 0

    const_prefactor = (2*(np.pi**4)*(G**2)*(M**(5/2))
                       *mB
                       )/((ma*mb*ms)**(3/2))
    term1 = LB/(Ls*((E0-EB)**(3/2))*((-EB)**(3/2)))

    arccosh_inside = (1+((2*current_R*(E0-EB))/(G*mB*ms))
                     )/((1+((2*M*(E0-EB)*(Ls**2))/((G**2)*((mB*ms)**3))))**(1/2)) 

    arccosh_term = np.arccosh(arccosh_inside)

    radical_term_inside = (((2*M*(E0-EB))/((G**2)*(ms**3)*(mB**3)))
                    *(2*current_R*G*M*(m**2)
                      +2*m*(E0-EB)*(current_R**2)-Ls**2))
        
    radical_term = np.sqrt(radical_term_inside)
        

    return const_prefactor*term1*(radical_term - arccosh_term)

@njit
def elliptic_dsigma_dEBdLBdCB(EBm, LBm, CB, mam, mbm, msm, alpha, E0m, L0m):
    EB = EBm*MSun
    LB = LBm*MSun
    ma = mam*MSun
    mb = mbm*MSun
    ms = msm*MSun
    E0 = E0m*MSun
    L0 = L0m*MSun
    m = ((ma + mb)*ms)/(ma + mb + ms)
    mB = ma + mb
    M = mB + ms
    
    Ls = np.sqrt(LB**2 - 2*LB*L0*CB + L0**2)
    current_R = (alpha*G*ma*mb)/(-2*EB)
    heaviside_term = -(Ls**2 -2*current_R*G*M*(m**2)
                        +2*m*(EB-E0)*(current_R**2))
    if heaviside_term <= 0:
        return 0

    const_prefactor = (2*(np.pi**4)*(G**2)*(M**(5/2))*mB)/((ma*mb*ms)**(3/2))
    term1_energies = (EB-E0)**(3/2) 
    term1 = LB/(Ls*(term1_energies)*((-EB)**(3/2)))

    inner_arccos = (1-((2*current_R*(EB - E0))/(G*mB*ms)))/(np.sqrt(1-((2*M*(EB-E0)*(Ls**2))/((G**2)*((mB*ms)**3)))))

    arccos_term = np.arccos(inner_arccos)

    inner_radical_term = ((-2*M*(EB - E0))/((G**2)*(ms**3)*(mB**3))
                            *(Ls**2 -2*current_R*G*M*(m**2)
                            +2*m*(EB-E0)*(current_R**2)))
    radical_term = np.sqrt(inner_radical_term)
        
    result = const_prefactor*term1*(arccos_term - radical_term)
    return result