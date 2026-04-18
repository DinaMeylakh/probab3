from probab3.commands.common.constants import *
from probab3.commands.common.formulas.general import *
from probab3.commands.common.data_classes import MassObj
from probab3.commands.common.formulas.inspiral import Inspiral
from scipy import stats, integrate, special
import random
import pickle
import os

class ClusterCalc:

    def __init__(self, rho_c, vrms, escape_vel):
        self.rho_c = rho_c
        self.vrms = vrms
        self.escape_vel = escape_vel

    def aHB(self, ma, mb, ms, k=1.4831497279200792):
        # the 3 is added since MB distirbution <v^2> = 3*sigma^2
        return (1-(1/k))*(G*ma*mb)/(3*(m(ma, mb, ms))*(self.vrms**2))

    def sample_vs(self, lower_limit, upper_limit, size=1):

        scale = self.vrms

        if upper_limit == np.inf:
            norm = 1
        else: 
            mb_cdf = lambda x: special.erf(x/(np.sqrt(2)*scale)) - np.sqrt(2/np.pi)*(x/scale)*np.exp(-(x**2)/(2*(scale**2)))
            norm = mb_cdf(upper_limit) - mb_cdf(lower_limit) 
        
        class MaxwellBotzmannDist(stats.rv_continuous):
            def _pdf(self, x):
                return (1/norm)*(np.sqrt(2/np.pi))*(x**2/(scale**3))*np.exp(-(x**2)/(2*(scale**2)))

        maxwell_boltzman_dist = MaxwellBotzmannDist(a=lower_limit, b=upper_limit)
        samples = maxwell_boltzman_dist.rvs(size=size)
        return samples

    def ejection_energy(self, ma, mb, ms):
        return (mB(ma, mb)*(1+(mB(ma, mb)/ms))*(self.escape_vel**2))/(2)

    def aB_for_tidal_disruption(self, ma, mb):
        return (mB(ma, mb)/self.rho_c)**(1/3)
    
    def sample_a_in_main_sequence(self, ma=MSun, mb=MSun):
        upper_domain_bound = self.aB_for_tidal_disruption(ma, mb) 
        lower_domain_bound = (star_radius(ma) + star_radius(mb))
        norm = np.log(upper_domain_bound) - np.log(lower_domain_bound)
        class InverseDist(stats.rv_continuous):
            def _pdf(self, x):
                return (1/norm)*(1/x)
        
        inverse_dist = InverseDist(name="inverse", a=lower_domain_bound, b=upper_domain_bound)        
        return inverse_dist.rvs(size=1)[0]

    def sample_e_in_thermal(self, lower_bound=0, upper_bound=1):

        norm = (1/2)*((upper_bound**2) - (lower_bound**2))
        class ThermalDist(stats.rv_continuous):
            def _pdf(self, x):
                return (1/norm)*x
        
        thermal_dist = ThermalDist(name="inverse", a=lower_bound, b=upper_bound)        
        return thermal_dist.rvs(size=1)[0]

    def sample_initial_binary(self, ma, mb):
        aB_sample = self.sample_a_in_main_sequence(ma, mb)
        eB_sample = self.sample_e_in_thermal(upper_bound=(1-((star_radius(ma) + star_radius(mb))/(aB_sample)))) 
        
        return eB_sample, aB_sample

    def stayed_in_cluster_aligned(self, m_BH):
        equal_mass_inspiral = Inspiral(m_BH, m_BH)
        v_recoil = equal_mass_inspiral.calc_v_recoil_abs(0, 0)

        if v_recoil >= self.escape_vel:
            return False
        return True 


class StarCluster(ClusterCalc):

    def __init__(self, Mc, rh, rc, average_star_mass=0.38*MSun):
        self.Mc = Mc
        self.rh = rh
        self.rc = rc
        self.average_star_mass = average_star_mass

        super().__init__(rho_c=self.core_density(),
                         vrms=self.calc_vc(),
                         escape_vel=self.calc_escape_vel())

    def potential(self, r):
        if r == 0:
            return -((4*np.pi*self.rho_c*G*(self.rc**2)*(self.rh**2))/((self.rh**2) - (self.rc**2)))*((1/2)*np.log((self.rh**2)/(self.rc**2))) 
        return -((4*np.pi*self.rho_c*G*(self.rc**2)*(self.rh**2))/((self.rh**2) - (self.rc**2)))*((self.rh/r)*np.arctan(r/self.rh)-(self.rc/r)*np.arctan(r/self.rc)+(1/2)*np.log((r**2 + self.rh**2)/(r**2 + self.rc**2)))

    def number_density(self, r):
        return self.density(r)/self.average_star_mass
    
    def density(self, r):
        return (self.rho_c)/((1+((r**2)/(self.rc**2)))*(1+((r**2)/(self.rh**2))))

    def mass_enclosed(self, r):
        return (((4*np.pi*(self.rc**2)*(self.rh**2)*self.rho_c)/((self.rh**2) - (self.rc**2)))*(self.rh*np.arctan(r/self.rh) - self.rc*np.arctan(r/self.rc)))
        
    def calc_vc(self):
        return ((6*G*self.Mc*(((np.pi**2)/8)-1))/(np.pi*self.rh))**(1/2)

    def calc_escape_vel(self):
        return 2*(((G*self.Mc*np.log(self.rh/self.rc))/(np.pi*(self.rh-self.rc)))**(1/2))

    def calc_escape_vel_from_potential(self, r=0):
        return (-2*self.potential(r))**(1/2)

    def core_density(self):
        return self.Mc*(self.rh+self.rc)/(2*(np.pi**2)*(self.rc**2)*(self.rh**2))
    
    def core_number_density(self):
        return self.core_density() / self.average_star_mass

    def sigma_squared(self, r):
        func = lambda r_tag: (1/(r_tag**2))*(self.mass_enclosed(r_tag))*(self.density(r_tag))
        integral_value = integrate.quad(func, r, 10*parsec)
        return ((3*G)/(self.density(r)))*integral_value[0]

    def binary_single_scattering_time_of_r(self, ma, mb, R, r):
        m_star = self.average_star_mass
        return np.sqrt(self.sigma_squared(r))*m_star/(2*np.pi*R*G*(ma + mb)*self.density(r))
    
    def binary_single_scattering_time(self, ma, mb, R):
        return self.vrms/(2*np.pi*R*G*(ma + mb)*self.core_number_density())



class TwoComponentCluster(ClusterCalc):
    """Cluster class for a cluster containing BHs sub cluster within a stellar cluster."""

    def __init__(self, Mc_stars, rh_stars, rc_stars, 
                 fraction_BH=10**(-3), average_star_mass=0.38*MSun, 
                 average_BH_mass=20*MSun, BH_concentration=3.5):
        self.stars_cluster = StarCluster(Mc_stars, rh_stars, rc_stars)
        Mc_BH = fraction_BH*average_BH_mass*(Mc_stars/average_star_mass)
        rh_BH = (Mc_BH/Mc_stars)*rh_stars
        rc_BH = (1/BH_concentration)*rh_BH
        self.BH_cluster = StarCluster(Mc_BH, rh_BH, rc_BH, average_star_mass=average_BH_mass)
        
        super(TwoComponentCluster, self).__init__(rho_c=(self.stars_cluster.rho_c + self.BH_cluster.rho_c),
                                                  vrms=self.calc_vc(),
                                                  escape_vel=self.calc_escape_vel())
    
    def calc_escape_vel(self, r=0):
        return (-2*(self.stars_cluster.potential(r) + self.BH_cluster.potential(r)))**(1/2)
    
    def calc_vc(self):
        return (self.sigma_squared(0))**(1/2)

    def density(self, r):
        return self.BH_cluster.density(r) + self.stars_cluster.density(r)

    def number_density(self, r):
        return self.BH_cluster.number_density(r) + self.stars_cluster.number_density(r)
    
    def core_number_density(self):
        return self.BH_cluster.core_number_density() + self.stars_cluster.core_number_density()
    
    def mass_enclosed(self, r):
        return self.BH_cluster.mass_enclosed(r) + self.stars_cluster.mass_enclosed(r)

    def sigma_squared(self, r):
        func = lambda r_tag: (1/(r_tag**2))*(self.mass_enclosed(r_tag))*(self.density(r_tag))
        integral_value = integrate.quad(func, r, 10*parsec)
        return ((3*G)/(self.density(r)))*integral_value[0]

    def binary_single_scattering_time(self, ma, mb, R):
        return self.BH_cluster.binary_single_scattering_time(ma, mb, R)

class UnequalMassCluster(TwoComponentCluster):

    pdmf_updated = False

    def update_pdmf(self):
        if self.pdmf_updated:
            return
        self.pdmf = ClusterCorePDMF(light_obj_num_dens=self.stars_cluster.core_number_density(),
                                    light_obj_disper=self.stars_cluster.vrms,
                                    massive_obj_num_dens=self.BH_cluster.core_number_density(),
                                    massive_obj_disper=self.BH_cluster.vrms)

        self.pdmf_updated = True

    def sample_star_mass(self) -> MassObj:
        return MassObj(mass=self.pdmf.sample_star_mass(), identity=MassIdentity.STAR.value)

    def sample_BH_mass(self) -> MassObj:
        return MassObj(mass=self.pdmf.sample_BH_mass(), identity=MassIdentity.BH.value)

    def sample_mass(self) -> MassObj:
        return self.pdmf.sample_mass_obj()

    def binary_single_scattering_time(self, ma, mb, R):
        return self.vrms/(2*np.pi*R*G*(ma + mb)*self.core_number_density())

class PDMF:

    imf_norm = 1/(0.6592+7.0518+0.0147+0.0025)
    prob_of_object = {
        MassIdentity.WD.value: imf_norm*0.6592,
        MassIdentity.STAR.value: imf_norm*7.0518,
        MassIdentity.BH.value: imf_norm*0.0147,
        MassIdentity.NS.value: imf_norm*0.0025
    }
    INITIAL_FINAL_PATH = os.path.join('probab3', 'commands', 'pre_calculation', 'inital_final_BH_mass', 'initial_final_interp.pck')

    def __init__(self):
        
        with open(self.INITIAL_FINAL_PATH, 'rb') as interp_file:
            loaded_initial_final_func = pickle.load(interp_file)
        
        self.initial_final_BH_relation = loaded_initial_final_func 


    def sample_mass_obj(self, size=1):
        mass_objs = []
        mass_identities = self.sample_mass_identity(size)
        sample_funcs = {
            MassIdentity.WD.value: self.sample_WD_mass,
            MassIdentity.STAR.value: self.sample_star_mass,
            MassIdentity.BH.value: self.sample_BH_mass,
            MassIdentity.NS.value: self.sample_NS_mass 
        }
        for mass_identity in mass_identities:
            relevant_sample_func = sample_funcs[mass_identity]
            mass_objs.append(MassObj(mass=relevant_sample_func(), identity=mass_identity))
        
        return mass_objs[0] if len(mass_objs) == 1 else mass_objs

    def sample_mass_identity(self, size=1):
        return random.choices(list(self.prob_of_object.keys()), list(self.prob_of_object.values()), k=size)

    def sample_WD_mass(self):
        return 0.6*MSun
    
    def sample_NS_mass(self):
        return 1.4*MSun
    
    def sample_star_mass(self):
        norm = 0.1418
        class StarsPDMF(stats.rv_continuous):
            def _pdf(self, x):
                if x < 0.5:
                    return (2*norm)*(x**(-1.3))
                return (norm)*(x**(-2.3))
        
        star_pdmf = StarsPDMF(name="star_pdmf", a=0.08, b=0.92)        
        star_mass_in_MSun = star_pdmf.rvs(size=1)[0]
        return star_mass_in_MSun*MSun 

    def sample_BH_mass(self):
        norm = 440
        class BHProjenitorIMF(stats.rv_continuous):
            def _pdf(self, x):
                return (norm)*(x**(-2.7))
        
        BH_projenitor_imf = BHProjenitorIMF(name="BH_pro_imf", a=25, b=110)
        BH_projenitor_in_MSun = BH_projenitor_imf.rvs(size=1)[0]
        BH_mass_in_MSun = self.initial_final_BH_relation(BH_projenitor_in_MSun) 
        return BH_mass_in_MSun*MSun


class ClusterCorePDMF(PDMF):

    def __init__(self, light_obj_num_dens, light_obj_disper, massive_obj_num_dens, massive_obj_disper):

        super(ClusterCorePDMF, self).__init__()

        all_objs_norm = light_obj_disper*light_obj_num_dens + massive_obj_num_dens*massive_obj_disper
        lighter_prob = (light_obj_disper*light_obj_num_dens)/all_objs_norm 
        massive_prob = (massive_obj_num_dens*massive_obj_disper)/all_objs_norm 

        self.prob_of_object = {
            MassIdentity.WD.value: lighter_prob*((0.6592)/(0.6592 + 7.0518)),
            MassIdentity.STAR.value: lighter_prob*((7.0518)/(0.6592 + 7.0518)),
            MassIdentity.BH.value: massive_prob*((0.0147)/(0.0147 + 0.0025)),
            MassIdentity.NS.value: massive_prob*((0.0025)/(0.0147 + 0.0025))
        }
