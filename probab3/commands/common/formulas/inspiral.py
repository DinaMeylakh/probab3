from probab3.commands.common.constants import *
from probab3.commands.common.formulas.general import *
import numpy as np

class Inspiral:

    def __init__(self, m1, m2, S1 = None, S2 = None):
        if m1 > m2:
            raise Exception("The masses should follow m1 <= m2")
        self.m1 = m1
        self.m2 = m2
        # The unitless spin parameter
        self.S1 = S1 if S1 else self.calc_spin_param(m1)
        # Take std as difference if the same m, else take the different mean values
        da = self.calc_spin_param_std(m1) if m1 == m2 else 0
        self.S2 = S2 if S2 else (self.calc_spin_param(m2) + da) 
        self.dm = (m1 - m2)/(m1 + m2)
        self.q = m1/m2
        self.ni = (m1*m2)/((m1 + m2)**2)

    
    @staticmethod
    def calc_spin_param(m, p1 = 0.86, p2=0.13, p3=29.5):
        return ((p1 - p2)/2)*np.tanh(p3 - m/MSun) + ((p1 + p2)/2)

    def calc_spin_param_std(self, m):
        alpha_max = self.calc_spin_param(m, p1=(0.86+0.06), p2=(0.13+0.13), p3=(29.5+8.5))
        alpha_min = self.calc_spin_param(m, p1=(0.86-0.06), p2=(0.13-0.13), p3=(29.5-8.5))
        return np.abs(alpha_max - alpha_min)/np.sqrt(12)
    
    def calc_S_tilde_parallel(self, theta1=0, theta2=0):
        return (self.S2*np.cos(theta2) + self.S1*np.cos(theta1)*(self.q**2))/((1 + self.q)**2)
    
    def calc_Delta_tilde_parallel(self, theta1=0, theta2=0):
        return (self.S2*np.cos(theta2) - self.q*self.S1*np.cos(theta1))/(1 + self.q)
    
    def calc_v_m(self):
        """Returns value in m/2"""
        A = -8712
        B = -6516
        C = 3907
        return (self.ni**2)*self.dm*(A + B*self.dm +C*(self.dm**4))*1000

    def calc_v_orth(self, theta1 = 0, theta2 = 0):
        H = 7499.115
        H_2a = -1.736510
        H_2b = -0.598144
        H_3a = -0.318117
        H_3b = -0.748613
        H_3c = -1.749784 
        H_3d = -0.011247
        H_3e = -0.920198
        H_4a = -0.434318
        H_4b = -1.716134
        H_4c = 0.619181
        H_4d = 1.633127
        H_4e = -2.253606
        H_4f = -0.028194

        Delta_par = self.calc_Delta_tilde_parallel(theta1, theta2)
        S_par = self.calc_S_tilde_parallel(theta1, theta2)

        return H*(self.ni**2)*(Delta_par +
                               + H_2a*S_par*self.dm +
                               + H_2b*Delta_par*S_par +
                               + H_3a*(Delta_par**2)*self.dm +
                               + H_3b*(S_par**2)*self.dm +
                               + H_3c*Delta_par*(S_par**2) +
                               + H_3d*(Delta_par**3) +
                               + H_3e*Delta_par*(self.dm**2) +
                               + H_4a*(S_par)*(Delta_par**2)*self.dm +
                               + H_4b*(S_par**3)*self.dm +
                               + H_4c*S_par*(self.dm**3) +
                               + H_4d*Delta_par*S_par*(self.dm**2) +
                               + H_4e*Delta_par*(S_par**3) +
                               + H_4f*(Delta_par**3)*S_par)*1000
    
    def calc_v_recoil_abs(self, theta1, theta2):
        return np.sqrt((self.calc_v_m()**2) + self.calc_v_orth(theta1, theta2)**2)


class Bounds():

    @staticmethod
    def Fofe(e):
        return (e**(12/19)/(1+e))*(1+(121/304)*(e**2))**(870/2299)
    
    @staticmethod
    def qB_EM(ma, mb, f=F_LIGO, ef=0.1):
        qf = (G*mB(ma, mb)/((np.pi*f)**2))**(1/3)
        return (qf/(2*Bounds.Fofe(ef)))*((425/304)**(870/2299))
        #return 2.7*(G*mB(ma, mb)/((np.pi*f)**2))**(1/3)

    @staticmethod
    def qB_merge(ma, mb, ms, E0, EB):
        const_term = (85*np.pi/3)**(2/7)*(G/2)
        energies_term = (EB*ma*mb/(((E0-EB)**3)))**(1/7)
        mass_term = ((mB(ma, mb)**4)*(ms**3)/(c**10))**(1/7)
        return const_term*energies_term*mass_term 

    @staticmethod
    def EB_GW_analytic(ma, mb, eB, T_ref):
        return ((-G*ma*mb)/2)*((5*(c**5)*
                (1 + 0.27*(eB**(10)) + 0.33*(eB**(20)) + 0.2*(eB**(1000)))*
                (1 - eB**2)**(7/2))/(4*T_ref*64*(G**3)*ma*mb*mB(ma, mb)))**(1/4)

    @staticmethod
    def EB_EM(ma, mb, eB, f=F_LIGO, ef=0.1):
        return (-G*ma*mb*(1-eB))/(2*Bounds.qB_EM(ma=ma, 
                                                 mb=mb,
                                                 f=f,
                                                 ef=ef))
