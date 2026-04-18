from typing import Optional
import random
from scipy.stats import rv_continuous
from scipy.integrate import quad
from matplotlib import pyplot as plt
import numpy as np
import dataclasses
from enum import Enum
import click
import abc
import logging
import json
import os
import time
from colorama import init as colorama_init
from colorama import Fore
from colorama import Style
import traceback

from probab3.commands.common.formulas import star_cluster, triple_phase_space 
from probab3.commands.common.formulas.general import *
from probab3.commands.common.general_code import * 
from probab3.commands.common.formulas.inspiral import Bounds
from probab3.commands.common.data_classes import *

logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)
colorama_init()


class BinaryTeritiarySampler():
    """A class to sample a binary and a teritiary star from a star cluster."""

    def __init__(self, alpha: float, Mc: int, rh: int, rc: int):
        """Initialize the BinaryTeritiarySampler. Use SI units.
            
        Args:
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
        """
        self.ms = None
        self.alpha = alpha
        self.init_cluster(Mc, rh, rc) 

    def init_cluster(self, Mc, rh, rc):
        """Initialize the cluster."""
        self.cluster = star_cluster.TwoComponentCluster(Mc_stars=Mc, rh_stars=rh, rc_stars=rc)

    def sample_CB(self, bounds=[[-1.0, 1.0]]):
        """Sample the cosine of the binary inclination angle.

        Default distribution is uniform [-1, 1], matches spherical clusters.

        Args:
            upper_bound (float): The upper bound of the distribution.

        Return:
            (float) The cosine of the binary inclination angle.
        """
        region_lengths = []
        for region in bounds:
            if region[0] > region[1]:
                raise Exception(f"Invalid bounds when sampling C: {bounds}")
            region_lengths.append(region[1] - region[0])

        # choose region by size
        region = bounds[0]
        if len(bounds) > 1:
            region_num = np.random.choice(len(bounds), p=[region_length/sum(region_lengths) for region_length in region_lengths])
            region = bounds[region_num]

            # sample from region and return region num
            return np.random.uniform(low=region[0], high=region[1]), region_num
            
        # sample from region
        return np.random.uniform(low=region[0], high=region[1]), 0

    def sample_ms(self):
        """Sample the mass object of the teritiary.

        Default behaviour is to always return the same ms as was given in initialization.
        
        Return:
            (MassObj) The mass object of the teritiary.
        """
        return self.ms
    
    def sample_vs(self, binary_state, ms, enforce_limits=False):
        """Sample the velocity of the teritiary.
        
        Args:
            ma (int): The mass of the first star in the binary.
            mb (int): The mass of the second star in the binary.
            ms (int): The mass of the teritiary.
            EB (int): The energy of the binary.

        Return:
            (float) The velocity of the teritiary.
        """
        ma_mass = binary_state.ma.mass
        mb_mass = binary_state.mb.mass
        lower_limit = 0
        if enforce_limits:
            lower_limit = vs_min_of_LB(ma_mass, mb_mass, ms.mass, binary_state.LB, binary_state.EB, self.alpha) 
        
        return np.abs(self.cluster.sample_vs(lower_limit=lower_limit, upper_limit=np.inf))

    def sample_b(self, ma, mb, ms, EB, vs, E0, LB, CBs):
        """Sample the impact parameter of the teritiary.
        
        Args:
            ma (MassObj): The mass of the first star in the binary.
            mb (MassObj): The mass of the second star in the binary.
            ms (MassObj): The mass of the teritiary.
            EB (int): The energy of the binary.
            vs (float): The velocity of the teritiary.
            E0 (int): The total energy of the triple system.
            LB (int): The angular momentum of the binary.
            CBs (float): The cosine of the teritiary-binary center of mass inclination angle.
        """
        Ls_lower_bound, Ls_upper_bound = Ls_bounds_from_L0max(ma.mass, mb.mass, ms.mass, E0, LB, CBs, self.alpha)
        domain_upper_E = (2*R(ma.mass, mb.mass, EB, self.alpha)*G*M(ma.mass, mb.mass, ms.mass))**(1/2)/vs
        domain_upper_L = Ls_upper_bound/(m(ma.mass, mb.mass, ms.mass)*vs)
        domain_upper = min(domain_upper_E, domain_upper_L)
        domain_lower = max(0, Ls_lower_bound/(m(ma.mass, mb.mass, ms.mass)*vs)) 
        if domain_upper <= domain_lower:
            logger.warning(f"Could not sample b, domain_upper={domain_upper} <= domain_lower={domain_lower}"
                           f" with Ls_lower_bound={Ls_lower_bound}, Ls_upper_bound={Ls_upper_bound}, "
                           f"domain_upper_E={domain_upper_E}, domain_upper_L={domain_upper_L}")
            return -1, False
        domain_length = domain_upper - domain_lower

        class Linear(rv_continuous):
            def _pdf(this, x):
                return x/((1/2)*domain_length**2)

        linear_dist = Linear(name="linear", a=domain_lower, b=domain_upper)        
        return linear_dist.rvs(size=1)[0], True
    
    def __str__(self):
        return f"Sampling DBS with alpha={self.alpha}, cluster={self.cluster}"

    def calculate_dbs_state(self, binary_state, teritiary_state) -> TripleState:
        """Calculate the triple state from the binary and teritiary states.
        
        Args:
            binary_state (BinaryState): The binary state.
            teritiary_state (TeritiaryState): The teritiary state.
        
        Return:
            (TripleState) The triple state.
        """
        CBs = teritiary_state.CBs
        new_L0 = (binary_state.LB**2 + teritiary_state.Ls**2 + 2*binary_state.LB*teritiary_state.Ls*CBs)**(1/2) 
        C0_new = calc_C0(L0=new_L0, LB=binary_state.LB, Ls=teritiary_state.Ls, CB=binary_state.CBd, Cs=teritiary_state.Cs)
        return TripleState(
            ma=binary_state.ma,
            mb=binary_state.mb,
            ms=teritiary_state.ms,
            EB=binary_state.EB,
            LB=binary_state.LB,
            CB=CBs,
            E0=binary_state.EB + teritiary_state.Es,
            L0=new_L0,
            C0=C0_new)
    
    def sample_teritiary_CBs(self, binary_state):
        
        sampled_Cs, _ = self.sample_CB()  
        sampled_CBs = calc_CBs(binary_state.CBd, sampled_Cs)

        return sampled_CBs, sampled_Cs

 
    def sample_teritiary(self, binary_state, crude_time=0) -> TeritiaryState:
        """Sample a new teritiary state impacting a binary.
        
        Args:
            binary_state (BinaryState): The binary state.
            crude_time (float): The time passed since the beginning of the simulation.
        
        Return:
            (TeritiaryState): The new teritiary state.
            (float): The new crude time.
        """

        def sample_ms_state(ms, binary_state, enforce_limits=False):
            vs = None
            b = None
            CBs = None
            E0_value = None
            sample_success = False
            try:
                vs = self.sample_vs(binary_state=binary_state, ms=ms, enforce_limits=enforce_limits)[0]
                m_triple = m(binary_state.ma.mass, binary_state.mb.mass, ms.mass) 
                Es_value = (1/2)*m_triple*(vs**2)
                E0_value = binary_state.EB + Es_value 
                CBs, Cs = self.sample_teritiary_CBs(binary_state)
                if E0_value <= 0:
                    b, sample_success = self.sample_b(ma=binary_state.ma, mb=binary_state.mb, ms=ms, EB=binary_state.EB, vs=vs, 
                                                      E0=E0_value, LB=binary_state.LB, CBs=CBs)
                    if not sample_success:
                        logger.error(f'{Fore.RED}{crude_time + time_passed}::Sampled bad teritiary: vs={vs}, E0={E0_value}, binary_state={binary_state}, ms={ms}, CBs={CBs}, b={b}.{Style.RESET_ALL}')
                        raise Exception(f"Sampled bad teritiary: vs={vs}, E0={E0_value}, binary_state={binary_state}, ms={ms}, CBs={CBs}, b={b}")
                else:
                    b = 0
                    sample_success = True

            except Exception as err:
                logger.error(f'{Fore.RED}Error when sampling teritiary: vs={vs}, E0={E0_value}, binary_state={binary_state}, ms={ms}, CBs={CBs}, b={b}.{Style.RESET_ALL}')
                raise err
            return TeritiaryState(ms=ms, Es=Es_value, Ls=m_triple*b*vs, CBs=CBs, Cs=Cs)
        
        ms = self.sample_ms()
        time_passed = self.cluster.binary_single_scattering_time(binary_state.ma.mass, binary_state.mb.mass, binary_state.aB())
        sample = sample_ms_state(ms=ms, binary_state=binary_state, enforce_limits=False)

        return sample, crude_time + time_passed


class TripleSampler(abc.ABC):
    """Sample binary properties from a triple system during its excursions or escape."""

    def __init__(self, triple_state: TripleState, rejection: bool=True, alpha: float=2.5):
        """Initialize the TripleSampler.
        
        Args:
            triple_state (TripleState): The triple state.
            rejection (bool): Whether to use rejection sampling.
            alpha (float): The alpha parameter for the chaotic triple range.
        """
        self.triple_state = triple_state
        self.alpha = alpha
        self.triple = triple_phase_space.TripleSystem(m1=int(triple_state.ma.mass), m2=int(triple_state.mb.mass), m3=int(triple_state.ms.mass),
                                                     E0=int(triple_state.E0), L0=int(triple_state.L0), alpha=float(self.alpha), calc_phase_space=True)
        logger.info(f"Constructed triple {self.triple}")
        self.rejection = rejection
    
    def sample_fs_dist(self, size=1):
        """Sample the final state distribution.
        
        Args:
            size (int): The amount of samples to return.

        Return:
            (List) The samples.
        """
        pass

    def sample_ims_dist(self, size=1):
        """Sample the intermediate state distribution.
        
        Args:
            size (int): The amount of samples to return.

        Return:
            (List) The samples.
        """
        pass

class EqualMassTripleSampler(TripleSampler):
    """An Equal Mass implementation of the TripleSampler class."""

    def sample_fs_dist(self, size=1):
        """Sample the final state distribution.
        
        Args:
            size (int): The amount of samples to return.

        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary) format.
        """
        sample_points = self.triple.hyper123.sample(size=size, rejection=self.rejection)
        return self.construct_binary_ms(sample_points)
        
    def sample_ims_dist(self, size=1):
        """Sample the intermediate state distribution.
        
        Args:
            size (int): The amount of samples to return.

        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary mass object) format.
        """
        sample_points = self.triple.ellip123.sample(size=size, rejection=self.rejection)
        return self.construct_binary_ms(sample_points)
    
    def construct_binary_ms(self, sample_points):
        """Construct (binary, teritiary) tuples from sampled points.

        Sample the sign of sin(I_0)sin(I_Bs) and calculate CBd from CB.
        
        Args:
            sample_points [(int, int, float),]: The sampled points from the triple ConeDist in (EB, LB, CB) format.
        
        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary mass object) format.
        """
        sample_binaries_mss = [] 
        for sample_point in sample_points:
            sample_binaries_mss.append((BinaryState(ma=self.triple_state.ma, 
                                                    mb=self.triple_state.mb,
                                                    EB=sample_point[1], 
                                                    LB=sample_point[2], 
                                                    CB=float(sample_point[0]), 
                                                    CBd=float(calc_CBd(self.triple_state.C0, float(sample_point[0])))), 
                                        self.triple_state.ms))
        return sample_binaries_mss
        

class UnequalMassTripleSampler(TripleSampler):
    """An Unequal Mass implementation of the TripleSampler class."""

    def sample_fs_dist(self, size=1):
        """Sample the final state distribution.
        
        Args:
            size (int): The amount of samples to return.

        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary) format.
        """
        dists = {
            "123": self.triple.hyper123, 
            "231": self.triple.hyper231,
            "312": self.triple.hyper312
        }
        probabilities = self.triple.get_fs_probabilities()
        return self.get_dist_binary_state(dists, probabilities, size)
    
    def sample_ims_dist(self, size=1):
        """Sample the intermediate state distribution.

        Args:
            size (int): The amount of samples to return.

        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary) format.
        """
        dists = {
            "123": self.triple.ellip123, 
            "231": self.triple.ellip231,
            "312": self.triple.ellip312
        }
        probabilities = self.triple.get_ims_probabilities()
        return self.get_dist_binary_state(dists, probabilities, size)

    def get_dist_binary_state(self, dists, probabilities, size):
        """Sample which mass iteration distribution to sample from and sample from it.
        
        Args:
            dists (dict): The mass iteration distributions to sample from.
            probabilities (dict): The probabilities of sampling from each distribution.
            size (int): The amount of samples to return.
        
        Return:
            [(BinaryState, MassObj),] The samples in (binary, teritiary) format.
        """
        dist_name_to_mass_ids = {
            "123": (self.triple_state.ma, self.triple_state.mb, self.triple_state.ms), 
            "231": (self.triple_state.mb, self.triple_state.ms, self.triple_state.ma),
            "312": (self.triple_state.ms, self.triple_state.ma, self.triple_state.mb)
        }
        sample_binaries_mss = [] 

        for sample_num in range(size):

            dist, chosen_dist_name = self.sample_which_dist(dists, probabilities)
            sample_point = dist.sample(size=1, rejection=self.rejection)[0]
            logger.debug(f"dist probabilities for {sample_num}/{size}::{probabilities}")                                                                                         
            logger.debug(f"chosen_dist_name for {sample_num}/{size}::{chosen_dist_name}")

            mass_ids = dist_name_to_mass_ids[chosen_dist_name]
            state_sample =  BinaryState(ma=mass_ids[0], mb=mass_ids[1], 
                                        EB=sample_point[1], 
                                        LB=sample_point[2], 
                                        CB=sample_point[0],
                                        CBd=calc_CBd(self.triple_state.C0, sample_point[0]))
            sample_binaries_mss.append((state_sample, mass_ids[2]))
        return sample_binaries_mss 

    def sample_which_dist(self, dists, probabilities):
        """Sample which mass iteration distribution to sample from.
        
        Args:
            dists (dict): The mass iteration distributions to sample from.
            probabilities (dict): The probabilities of sampling from each distribution.

        Return:
            (ConeDist, str) The sampled distribution and its name.
        """
        chosen_dist_name = np.random.choice(a=list(dists.keys()), p=probabilities)
        return dists[chosen_dist_name], chosen_dist_name


class ClusterBinaryTeritiarySampler(BinaryTeritiarySampler):
    """A class to evolve a dynamical binary sequence inside a two-component star cluster."""
    SCRAMBLE_CHUNK_SIZE = 50

    def __init__(self, alpha: float, Mc: int, rh: int, rc: int, 
                 rejection: bool=True, sampler_class: TripleSampler=EqualMassTripleSampler):
        """Initialize the ClusterBinaryTeritiarySampler.
        
        Args:
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            rejection (bool): Whether to use rejection sampling.
            sampler_class (TripleSampler): The triple sampler class to use. Default is EqualMassTripleSampler.
        """
        super(ClusterBinaryTeritiarySampler, self).__init__(alpha, Mc, rh, rc)
        self.sampler_class = sampler_class
        self.binary_init_state = None # expected to be initialized in child class
        self.rejection = rejection

    def get_ejection_energy(self, fs_binary_state, fs_ms):
        """Get the ejection energy of the binary from the cluster.
        
        Args:
            fs_binary_state (BinaryState): The binary state.
            fs_ms (MassObj): The mass object of the teritiary.

        Return:
            (float) The ejection energy.
        """
        return self.cluster.ejection_energy(ma=fs_binary_state.ma.mass, mb=fs_binary_state.mb.mass, ms=fs_ms.mass)

    def get_evolution_state(self, binary_state: BinaryState, 
                            ms: MassObj, evolution_state: str, 
                            crude_time: int,  dbs_state: Optional[TripleState]=None, 
                            teritiary_state: Optional[TeritiaryState]=None, 
                            P_dis: Optional[float]=None,
                            scramble_num: Optional[int]=None,
                            num_of_scrambles: Optional[int]=None) -> MetastableState:
        """Get the metastable state from the evolution data.

        Args:
            binary_state (BinaryState): The binary state.
            ms (MassObj): The mass object of the teritiary.
            evolution_state (DBSEvolutionState): The evolution state.
            crude_time (int): The crude time of the evolution.
            dbs_state (TripleState): Optional. The triple state.
            teritiary_state (TeritiaryState): Optional. The teritiary state.
            P_dis (float): Optional. The disintegration probability.
            scramble_num (int): Optional. The scramble number.
            num_of_scrambles (int): Optional. The number of scrambles.
        
        Return:
            (MetastableState) The metastable state.
        """

        EB_EM = None
        EB_merger = None
        dbs_kick = None
        ejection_energy = None
        qB_merge = None
        qB_EM = None
        ims_period = None

        if evolution_state in [DBSEvolutionState.FS.value, DBSEvolutionState.EJECTED_FS.value, DBSEvolutionState.MANUAL_EXCHANGE.value]:

            ejection_energy = self.cluster.ejection_energy(ma=binary_state.ma.mass, mb=binary_state.mb.mass, ms=ms.mass)

            if evolution_state in [DBSEvolutionState.FS.value, DBSEvolutionState.MANUAL_EXCHANGE.value]:
                T_ref = self.cluster.binary_single_scattering_time(
                    ma=binary_state.ma.mass, mb=binary_state.mb.mass, R=binary_state.aB())
                dbs_kick = dbs_state.E0 - binary_state.EB
                
            if evolution_state == DBSEvolutionState.EJECTED_FS.value: 
                T_ref = HUBBLE_TIME

            EB_merger = Bounds.EB_GW_analytic(ma=binary_state.ma.mass, mb=binary_state.mb.mass, eB=binary_state.eB(), 
                                            T_ref=T_ref)
            EB_EM = Bounds.EB_EM(ma=binary_state.ma.mass, mb=binary_state.mb.mass, eB=binary_state.eB())
            qB_merge = qB_from_EBeB(ma=binary_state.ma.mass, mb=binary_state.mb.mass, EB=EB_merger, eB=binary_state.eB())
            qB_EM = qB_from_EBeB(ma=binary_state.ma.mass, mb=binary_state.mb.mass, EB=EB_EM, eB=binary_state.eB())
        elif evolution_state == DBSEvolutionState.IMS.value:
            qB_merge = Bounds.qB_merge(ma=binary_state.ma.mass, mb=binary_state.mb.mass, ms=ms.mass, 
                                       E0=dbs_state.E0, EB=binary_state.EB)
            qB_EM = Bounds.qB_EM(ma=binary_state.ma.mass, mb=binary_state.mb.mass)
            ims_period = period(mB(binary_state.ma.mass, binary_state.mb.mass), ms.mass, (dbs_state.E0 - binary_state.EB))

        a_s = None
        if evolution_state == DBSEvolutionState.IMS.value:
            a_s = aB(ms.mass, mB(binary_state.ma.mass, binary_state.mb.mass), dbs_state.E0 - binary_state.EB)

        return MetastableState(dbs_state=dbs_state,
                               binary_state=binary_state,
                               teritiary_state=teritiary_state,
                               crude_time=crude_time,
                               state=evolution_state, 
                               ms=ms,
                               r_TDE = binary_state.r_TDE() if evolution_state != DBSEvolutionState.NEW_MS.value else None,
                               r_collision = binary_state.r_collision() if evolution_state != DBSEvolutionState.NEW_MS.value else None, 
                               qB=binary_state.qB() if evolution_state != DBSEvolutionState.NEW_MS.value else None,
                               aB=binary_state.aB() if evolution_state != DBSEvolutionState.NEW_MS.value else None,
                               qB_merger=qB_merge,
                               qB_EM=qB_EM,
                               period=ims_period,
                               scramble_num=scramble_num,
                               dbs_kick=dbs_kick,
                               ejection_energy=ejection_energy,
                               EB_merger=EB_merger, 
                               EB_EM=EB_EM,
                               P_dis=P_dis,
                               num_of_scrambles=num_of_scrambles,
                               a_s=a_s,
                               a_s_cluster_tides=self.cluster.aB_for_tidal_disruption(ms.mass, mB(binary_state.ma.mass, binary_state.mb.mass)),
                               EB_cluster_tides=EB(binary_state.ma.mass, binary_state.mb.mass, self.cluster.aB_for_tidal_disruption(binary_state.ma.mass, binary_state.mb.mass)))
    
    def get_evolution_state_resolution(self, evolution_state: MetastableState, 
                                       merger_condition: bool, em_condition: bool, merger_value: MetastableFinalState, 
                                       em_value: MetastableFinalState, collision_value: MetastableFinalState, 
                                       tde_value: MetastableFinalState):
        """Get the resolution of the evolution state.

        Args:
            evolution_state (MetastableState): The evolution state.
            merger_condition (bool): Whether the merger condition is met.
            em_condition (bool): Whether the eccentric merger condition is met.
            merger_value (MetastableFinalState): The merger metastable final state name.
            em_value (MetastableFinalState): The eccentric merger metastable final state name.
            collision_value (MetastableFinalState): The collision metastable final state name.
            tde_value (MetastableFinalState): The TDE metastable final state name.
        
        Return:
            (MetastableFinalState) The resolution of the evolution state.
        """

        binary = evolution_state.binary_state
        binary_identities = frozenset({binary.ma.identity, binary.mb.identity})

        # Check for GW merger
        if binary_identities in {frozenset({MassIdentity.BH.value, MassIdentity.BH.value}),
                                 frozenset({MassIdentity.BH.value, MassIdentity.NS.value}),
                                 frozenset({MassIdentity.NS.value, MassIdentity.NS.value}),
                                 frozenset({MassIdentity.BH.value, MassIdentity.WD.value}),
                                 frozenset({MassIdentity.NS.value, MassIdentity.WD.value}),
                                 frozenset({MassIdentity.WD.value, MassIdentity.WD.value})}:
            if merger_condition:
                if em_condition:
                    return em_value 
                return merger_value 

        # Check for collision
        if binary_identities in {frozenset({MassIdentity.STAR.value, MassIdentity.STAR.value}),
                                 frozenset({MassIdentity.WD.value, MassIdentity.STAR.value})}:
            if evolution_state.qB < evolution_state.r_collision:
                return collision_value

        # Check for TDEs 
        if binary_identities in {frozenset({MassIdentity.BH.value, MassIdentity.STAR.value}),
                                 frozenset({MassIdentity.NS.value, MassIdentity.STAR.value}),
                                 frozenset({MassIdentity.WD.value, MassIdentity.STAR.value}),
                                 frozenset({MassIdentity.BH.value, MassIdentity.WD.value}),
                                 frozenset({MassIdentity.NS.value, MassIdentity.WD.value}),
                                 frozenset({MassIdentity.WD.value, MassIdentity.WD.value})}:
            if evolution_state.qB < evolution_state.r_TDE:
                return tde_value        
        
        return None
    
    def get_resolution(self, evolution: list, evolution_state: MetastableState):
        """Get the resolution of the evolution state and append it to the evolution list.
        
        Args:
            evolution ([MetastableState,]): The evolution list.
            evolution_state (MetastableState): The last evolution state.

        Return:
            (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]) 
                The resolution and the evolution list with {"final": resolution, "crude_time": crude time} appended at the end.
        """

        resolution = None 

        if evolution_state.state == DBSEvolutionState.NEW_MS.value:
            return None, evolution
        
        if evolution_state.state == DBSEvolutionState.IMS.value:
            resolution = self.get_evolution_state_resolution(evolution_state=evolution_state, 
                                                             merger_condition=evolution_state.qB < evolution_state.qB_merger,
                                                             em_condition=evolution_state.qB < evolution_state.qB_EM,
                                                             merger_value=MetastableFinalState.IMS_MERGER.value,
                                                             em_value=MetastableFinalState.IMS_EM_MERGER.value,
                                                             collision_value=MetastableFinalState.IMS_COLLISION.value,
                                                             tde_value=MetastableFinalState.IMS_TDE.value)
            
            # Check for binary tidal disruption
            if not resolution and evolution_state.binary_state.EB > evolution_state.EB_cluster_tides:
                resolution = MetastableFinalState.IMS_BINARY_TD.value
            
        
        if evolution_state.state == DBSEvolutionState.FS.value:
            resolution = self.get_evolution_state_resolution(evolution_state=evolution_state,
                                                             merger_condition=evolution_state.binary_state.EB < evolution_state.EB_merger,
                                                             em_condition=evolution_state.binary_state.EB < evolution_state.EB_EM,
                                                             merger_value=MetastableFinalState.FS_MERGER.value,
                                                             em_value=MetastableFinalState.FS_EM_MERGER.value,
                                                             collision_value=MetastableFinalState.FS_COLLISION.value,
                                                             tde_value=MetastableFinalState.FS_TDE.value)
            # Check for binary tidal disruption 
            if not resolution and evolution_state.binary_state.EB > evolution_state.EB_cluster_tides:
                    resolution = MetastableFinalState.FS_BINARY_TD.value
            
        if evolution_state.state == DBSEvolutionState.EJECTED_FS.value:
            resolution = self.get_evolution_state_resolution(evolution_state=evolution_state,
                                                             merger_condition=evolution_state.binary_state.EB < evolution_state.EB_merger,
                                                             em_condition=evolution_state.binary_state.EB < evolution_state.EB_EM,
                                                             merger_value=MetastableFinalState.EJECTED_FS_MERGER.value,
                                                             em_value=MetastableFinalState.EJECTED_FS_EM.value,
                                                             collision_value=MetastableFinalState.EJECTED_FS_COLLISION.value,
                                                             tde_value=MetastableFinalState.EJECTED_FS_TDE.value)

        if resolution:
            evolution.append({"final": resolution, "crude_time": evolution_state.crude_time})
        
        return resolution, evolution
    

    def keep_heaviest_binary(self, binary_state: BinaryState, ms: MassObj):
        """Keep the heaviest binary and teritiary star.
        
        Args:
            binary_state (BinaryState): The binary state.
            ms (MassObj): The mass object of the teritiary.

        Return:
            (BinaryState) The binary state and the evolution state.
        """

        if binary_state.ma.mass > ms.mass and binary_state.mb.mass > ms.mass:
            return binary_state
        
        triple_bodies = [binary_state.ma, binary_state.mb, ms]
        triple_bodies.sort(key=lambda x: x.mass, reverse=True)
        new_ma = triple_bodies[0]
        new_mb = triple_bodies[1]

        return BinaryState(ma=new_ma,
                           mb=new_mb,
                           EB=binary_state.EB,
                           LB=binary_state.LB,
                           CB=binary_state.CB,
                           CBd=binary_state.CBd)


    def evolve_binary(self) -> MetastableFinalState: 
        """Evolve the binary inside the cluster. This is the heart of the code.

        Return:
            (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]) 
                The resolution and the evolution list with 
                    {"final": resolution, "crude_time": crude time} appended at the end.
        """
        dbs_kick = 0
        cluster_ejection_energy = dbs_kick + 1
        binary_state = self.binary_init_state
        crude_evolution_time = 0
        evolution = []
        triple_td_sampled = False

        logger.info(f"Starting, binary_state = {binary_state}")
       
        # While the kick the binary recieved is smaller than the cluster ejection energy, sample a teritiary and evolve the triple.
        while dbs_kick < cluster_ejection_energy:
            triple_td_sampled = False
            teritary_state, new_crude_time = self.sample_teritiary(binary_state=binary_state, crude_time=float(crude_evolution_time))
            if not teritary_state or new_crude_time >= HUBBLE_TIME:
                evolution.append({"final": MetastableFinalState.IN_CLUSTER_BINARY.value,
                                  "crude_time": new_crude_time})
                return MetastableFinalState.IN_CLUSTER_BINARY.value, evolution
             
            dbs_state = self.calculate_dbs_state(binary_state, teritary_state) 
            crude_evolution_time = new_crude_time

            logger.info(f"{Fore.GREEN}{crude_evolution_time}::sampled teritiary::{teritary_state} with binary {binary_state}{Style.RESET_ALL}")

            if dbs_state.E0 >= 0:
                logger.info(f"{Fore.LIGHTCYAN_EX}{crude_evolution_time}::triple E0 > 0 and (Es/EB)^(1/2) = {np.sqrt(-teritary_state.Es/binary_state.EB)} {Style.RESET_ALL}")
                evolution_state = self.get_evolution_state(binary_state=binary_state, teritiary_state=teritary_state, 
                                                           ms=teritary_state.ms, dbs_state=dbs_state,
                                                           evolution_state=DBSEvolutionState.NEW_MS.value, 
                                                           crude_time=crude_evolution_time) 
                evolution.append(evolution_state)

                if np.sqrt(-teritary_state.Es/binary_state.EB) < 1.5:
                    # If the triple is not chaotic and not scattered, keep the heaviest binary
                    binary_state = self.keep_heaviest_binary(binary_state, teritary_state.ms)
                    evolution_state = self.get_evolution_state(binary_state=binary_state, teritiary_state=teritary_state,
                                                               ms=teritary_state.ms, dbs_state=dbs_state,
                                                               evolution_state=DBSEvolutionState.MANUAL_EXCHANGE.value,
                                                               crude_time=crude_evolution_time)
                    evolution.append(evolution_state)
                    continue
                
                # Teritiary ionized the binary
                evolution.append({"final": MetastableFinalState.IONIZED_TRIPLE.value,
                                  "crude_time": crude_evolution_time})
                return MetastableFinalState.IONIZED_TRIPLE.value, evolution


            logger.info(f"{Fore.GREEN}{crude_evolution_time}::initializing triple::{dbs_state}{Style.RESET_ALL}")                                                                                         

            sampler = self.sampler_class(dbs_state, rejection=self.rejection, alpha=self.alpha)
            triple = sampler.triple
            disintegration_prob = triple.disintegration_probability()
            num_of_scrambles = np.random.geometric(p=disintegration_prob)
            
            evolution_state = self.get_evolution_state(binary_state=binary_state, teritiary_state=teritary_state, 
                                                       ms=teritary_state.ms, dbs_state=dbs_state,
                                                       evolution_state=DBSEvolutionState.NEW_MS.value, 
                                                       crude_time=crude_evolution_time,
                                                       P_dis=disintegration_prob,
                                                       num_of_scrambles=num_of_scrambles) 
            evolution.append(evolution_state)

            logger.info(f"{Fore.LIGHTYELLOW_EX}{crude_evolution_time}::evolution state::{evolution_state}{Style.RESET_ALL}")
            logger.info(f"{Fore.LIGHTYELLOW_EX}{crude_evolution_time}::scrambles::num_of_scrambles is {num_of_scrambles} for p={triple.disintegration_probability()}{Style.RESET_ALL}")
            num_of_ims_scrambles = num_of_scrambles - 1
            scrambles_chunks_nums = [self.SCRAMBLE_CHUNK_SIZE for _ in range(num_of_ims_scrambles // self.SCRAMBLE_CHUNK_SIZE)] + [num_of_ims_scrambles % self.SCRAMBLE_CHUNK_SIZE]            

            total_scramble_num = 1
            for num_chunk_scrambles in scrambles_chunks_nums:
                if num_chunk_scrambles < 1 or triple_td_sampled:
                    continue
                
                scrambles = sampler.sample_ims_dist(size=num_chunk_scrambles) 

                for scramble_num in range(0, num_chunk_scrambles):
                    if crude_evolution_time >= HUBBLE_TIME:
                        evolution.append({"final": MetastableFinalState.IN_CLUSTER_BINARY.value,
                                          "crude_time": crude_evolution_time})
                        return MetastableFinalState.IN_CLUSTER_BINARY.value, evolution
                
                    scramble, scramble_ms = scrambles[scramble_num]
                    evolution_ims_state = self.get_evolution_state(binary_state=scramble, 
                                                                   ms=scramble_ms, dbs_state=dbs_state,
                                                                   evolution_state=DBSEvolutionState.IMS.value, 
                                                                   crude_time=crude_evolution_time,
                                                                   scramble_num=total_scramble_num + scramble_num)
                        
                    logger.info(f"{Fore.YELLOW}{crude_evolution_time}::evolution_ims_state::{evolution_ims_state}{Style.RESET_ALL}")

                    evolution.append(evolution_ims_state)

                    # Check for GW merger/collision/TDE
                    resolution, evolution = self.get_resolution(evolution=evolution, evolution_state=evolution_ims_state)
                    if resolution:
                        return resolution, evolution

                    elif evolution_ims_state.a_s > evolution_ims_state.a_s_cluster_tides:
                        # Triple is tidally disrupted
                        crude_evolution_time = crude_evolution_time + (1/2)*evolution_ims_state.period
                        triple_td_state = self.get_evolution_state(binary_state=scramble,
                                                                   ms=scramble_ms, dbs_state=dbs_state,
                                                                   evolution_state=DBSEvolutionState.TRIPLE_TD.value, 
                                                                   crude_time=crude_evolution_time)
                        evolution.append(triple_td_state)
                        logger.info(f"{Fore.LIGHTBLUE_EX}{crude_evolution_time}::triple_td_state::{triple_td_state}{Style.RESET_ALL}")
                        binary_state = scramble 
                        triple_td_sampled = True
                        break

                    crude_evolution_time = crude_evolution_time + evolution_ims_state.period
                
                total_scramble_num = total_scramble_num + num_chunk_scrambles 

            if triple_td_sampled:
                continue
                
            fs_binary_state, fs_ms = sampler.sample_fs_dist()[0]
            evolution_fs_state = self.get_evolution_state(binary_state=fs_binary_state, 
                                                          ms=fs_ms, dbs_state=dbs_state,
                                                          evolution_state=DBSEvolutionState.FS.value, 
                                                          crude_time=crude_evolution_time)            

            evolution.append(evolution_fs_state)
            logger.info(f"{Fore.BLUE}{crude_evolution_time}::fs_final_state::{evolution_fs_state}{Style.RESET_ALL}")
            
            # Check for GW merger/collision/TDE
            resolution, evolution = self.get_resolution(evolution=evolution, evolution_state=evolution_fs_state)
            if resolution:
                return resolution, evolution

            dbs_kick = dbs_state.E0 - fs_binary_state.EB 
            cluster_ejection_energy = self.get_ejection_energy(fs_binary_state, fs_ms)
            logger.info(f"{Fore.MAGENTA}{crude_evolution_time}::dbs_kick = {dbs_kick} ? {cluster_ejection_energy} = cluster ejection energy{Style.RESET_ALL}")
            binary_state = fs_binary_state

        fs_ej_state, fs_ms = sampler.sample_fs_dist()[0]
        evolution_fs_ej_state = self.get_evolution_state(binary_state=fs_ej_state, 
                                                           ms=fs_ms, dbs_state=dbs_state,
                                                           evolution_state=DBSEvolutionState.EJECTED_FS.value, 
                                                           crude_time=crude_evolution_time)            
        evolution.append(evolution_fs_ej_state)
        logger.info(f"{Fore.MAGENTA}{crude_evolution_time}::evolution_fs_ej_state::{evolution_fs_ej_state}")

        # Check for GW merger/collision/TDE
        resolution, evolution = self.get_resolution(evolution=evolution, evolution_state=evolution_fs_ej_state)
        if resolution:
            return resolution, evolution
        
        # End with ejected binary
        evolution.append({"final": MetastableFinalState.EJECTED_BINARY.value,
                          "crude_time": crude_evolution_time})
        return MetastableFinalState.EJECTED_BINARY.value, evolution


class StarsBinaryToBHsBinary(ClusterBinaryTeritiarySampler):
    """A class to sample a star binary and a teritiary BH from a star cluster and evolve it to a binary of BHs."""

    def __init__(self, Mc: int, rh: int, rc: int, star_mass: int = MSun, BH_mass: int = 20*MSun, alpha: float = 2.5, rejection: bool=True):
        """Initialize the StarsBinaryToBHsBinary.
        
        Args:
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            star_mass (int): The mass of the stars in the cluster.
            BH_mass (int): The mass of the BHs in the cluster.
            alpha (float): The alpha parameter for the chaotic triple range.
            rejection (bool): Whether to use rejection sampling.
        """
        self.BH_mass = BH_mass
        self.star_mass = star_mass
        super(StarsBinaryToBHsBinary, self).__init__(alpha, Mc, rh, rc)
        self.rejection = rejection
        self.binary_state = self.binary_init_state
    
    def get_star_mass(self):
        return MassObj(mass=self.star_mass, identity=MassIdentity.STAR.value)
    
    def get_BH_mass(self):
        return MassObj(mass=self.BH_mass, identity=MassIdentity.BH.value)
    
    def init_stars_binary(self):
        """Initialize the first stars binary."""
        CB_sampled, _ = self.sample_CB() 
        ma = self.get_star_mass()
        mb = self.get_star_mass()
        while star_radius(ma.mass) + star_radius(mb.mass) > self.cluster.aB_for_tidal_disruption(ma.mass, mb.mass):
            logger.info(f"Sampling initial star masses again. Previously: ma={ma} mb={mb}")
            ma = self.get_star_mass()
            mb = self.get_star_mass()
        self.e_in, self.a_in = self.cluster.sample_initial_binary(ma.mass, mb.mass)
        self.binary_init_state = BinaryState(ma=ma, 
                                             mb=mb, 
                                             EB=EB(ma=ma.mass, mb=mb.mass, aB=self.a_in), 
                                             LB=LB_from_aBeB(ma=ma.mass, mb=mb.mass, aB=self.a_in, eB=self.e_in), 
                                             CB=CB_sampled,
                                             CBd=CB_sampled)
        
    def sample_ms(self):
        """Sample the mass object of the teritiary.
        
        Return:
            (MassObj) The mass object of the teritiary.
        """
        return self.get_BH_mass()
    
    def sample_vs(self, binary_state, ms, enforce_limits=False):
        """Sample the velocity of the teritiary.
        
        Args:
            ma (int): The mass of the first star in the binary.
            mb (int): The mass of the second star in the binary.
            ms (int): The mass of the teritiary.
            EB (int): The energy of the binary.

        Return:
            (float) The velocity of the teritiary.
        """
        ma_mass = binary_state.ma.mass
        mb_mass = binary_state.mb.mass
        upper_limit = ((-2*binary_state.EB)/m(ma_mass, mb_mass, ms.mass))**(1/2)
        lower_limit = 0
        if enforce_limits:
            lower_limit = vs_min_of_LB(ma_mass, mb_mass, ms.mass, binary_state.LB, binary_state.EB, self.alpha) 
        
        return np.abs(self.cluster.sample_vs(lower_limit=lower_limit, upper_limit=upper_limit))
    
    def sample_BH_and_evolve_until_outcome(self, allowed_outcome_objects=None, crude_time=0):
        """Sample a BH and evolve the triple system until the outcome is in the allowed_outcome_objects.
        
        Args:
            allowed_outcome_objects ([{MassIdentity, MassIdentity},]): The allowed outcome objects to stop the evolution at.
        
        Return:
            (Bool, [MetastableState,]) Flag if binary remains, the evolution states.
        """
        crude_evolution_time = crude_time

        if allowed_outcome_objects is None:
            allowed_outcome_objects = [{MassIdentity.STAR.value, MassIdentity.BH.value}]

        evolution = []
        outcome_objects = {self.binary_state.ma.identity, self.binary_state.mb.identity}

        while outcome_objects not in allowed_outcome_objects:
            # Introduce BH until the BH will be in the leftover binary
            BH_sampled, new_crude_time = self.sample_teritiary(self.binary_state, crude_evolution_time)
            if not BH_sampled or new_crude_time >= HUBBLE_TIME:
                evolution.append({"final": MetastableFinalState.IN_CLUSTER_BINARY.value,
                                  "crude_time": new_crude_time})
                logger.info(f"{Fore.CYAN}{crude_evolution_time}::Pre-evolution ended with {MetastableFinalState.IN_CLUSTER_BINARY.value}{Style.RESET_ALL}")
                return False, evolution 
            triple_state = self.calculate_dbs_state(self.binary_state, BH_sampled)
            crude_evolution_time = new_crude_time

            evolution_state = self.get_evolution_state(binary_state=self.binary_state, teritiary_state=BH_sampled, 
                                                       ms=BH_sampled.ms, dbs_state=triple_state,
                                                       evolution_state=DBSEvolutionState.NEW_MS.value, 
                                                       crude_time=crude_evolution_time) 
            evolution.append(evolution_state)

            logger.info(f"{Fore.GREEN} NEW MS evolution_state::{evolution_state}{Style.RESET_ALL}")
            
            #sample which binary survives from probabilities
            sampler = UnequalMassTripleSampler(triple_state, self.alpha)
            binary_state_sample, ms_obj = sampler.sample_fs_dist()[0]
            evolution_fs_state = self.get_evolution_state(binary_state=binary_state_sample, 
                                                          ms=ms_obj, dbs_state=triple_state,
                                                          evolution_state=DBSEvolutionState.FS.value, 
                                                          crude_time=crude_evolution_time)            

            evolution.append(evolution_fs_state)
            logger.info(f"{Fore.BLUE}{crude_evolution_time}::FS evolution_state::{evolution_fs_state}{Style.RESET_ALL}")
            self.binary_state = binary_state_sample
            outcome_objects = {binary_state_sample.ma.identity, binary_state_sample.mb.identity}

            # Check for GW merger/collision/TDE
            resolution, evolution = self.get_resolution(evolution=evolution, evolution_state=evolution_fs_state)
            if resolution:
                logger.info(f"{Fore.CYAN}{crude_evolution_time}::Pre-evolution ended with {resolution}{Style.RESET_ALL}")
                return False, evolution 

        return True, evolution 
    
    def evolve_stars_to_BHs(self):
        """Evolve the stars binary to a binary of BHs.
        
        Return:
            ([MetastableState,], MetastableState) The evolution states and the BH-star outcome metastable state.
        """
        failed_pre_evolutions = []
        sample_success = False
        while not sample_success:
            logger.info(f"{Fore.CYAN}Sampling new initial binary{Style.RESET_ALL}") 
            self.init_stars_binary()
            self.binary_state = self.binary_init_state
            logger.info(f"{Fore.CYAN}Sampled initial binary::{self.binary_state}{Style.RESET_ALL}")
            binary_remains, evolution1 = self.sample_BH_and_evolve_until_outcome(allowed_outcome_objects=[{MassIdentity.STAR.value, 
                                                                                                       MassIdentity.BH.value}])
            logger.info(f"Recieved first pre-evolution::binary remained={binary_remains}")
            if binary_remains:
                sample_success, evolution2 = self.sample_BH_and_evolve_until_outcome(allowed_outcome_objects=[{MassIdentity.BH.value, 
                                                                                                       MassIdentity.BH.value}])
                logger.info(f"Recieved second pre-evolution::binary remained={sample_success}")
                if not sample_success:
                    failed_pre_evolutions.append(evolution1 + evolution2[1:])
            else:
                failed_pre_evolutions.append(evolution1)
            

        logger.info(f"{Fore.CYAN}Pre-evolution ened successfully{Style.RESET_ALL}")
        return evolution1 + evolution2[1:], evolution2[0], failed_pre_evolutions 


class OneTripleSampler:

    def __init__(self, ma, mb, ms, E0, L0, alpha: float, 
                 rejection: bool=True, sampler_class: TripleSampler=UnequalMassTripleSampler, crude_time: int=0):
        
        self.triple_state = TripleState(ma=ma,
                                        mb=mb,
                                        ms=ms,
                                        E0=E0, L0=L0, C0=0)
        
        self.rejection = rejection
        self.sampler = sampler_class(triple_state=self.triple_state, rejection=self.rejection, alpha=alpha)
        self.crude_time = crude_time

    def evolve_triple(self):

        evolution = []
        triple = self.sampler.triple
        disintegration_prob = triple.disintegration_probability()
        num_of_scrambles = np.random.geometric(p=disintegration_prob)
        crude_time = self.crude_time

        evolution_state = MetastableState(dbs_state=self.triple_state,
                               crude_time=crude_time,
                               state=DBSEvolutionState.NEW_MS.value, 
                               num_of_scrambles=num_of_scrambles,
                               P_dis=disintegration_prob)
            
        evolution.append(evolution_state)

        logger.info(f"{Fore.GREEN}{crude_time}::initializing triple::{self.triple_state}{Style.RESET_ALL}") 

        scrambles = self.sampler.sample_ims_dist(size=num_of_scrambles)

        for scramble_num, scramble_tuple in enumerate(scrambles):
            scramble, scramble_ms = scramble_tuple
            evolution_ims_state = MetastableState(binary_state=scramble,
                                                  ms=scramble_ms,
                                                  dbs_state=self.triple_state,
                                                  crude_time=crude_time,
                                                  state=DBSEvolutionState.IMS.value, 
                                                  scramble_num=scramble_num,
                                                  period=period(mB(scramble.ma.mass, scramble.mb.mass), 
                                                                scramble_ms.mass, 
                                                                (self.triple_state.E0 - scramble.EB)),
                                                  aB=scramble.aB(),
                                                  qB=scramble.qB())
                        
            logger.info(f"{Fore.YELLOW}{crude_time}::evolution_ims_state::{evolution_ims_state}{Style.RESET_ALL}")

            evolution.append(evolution_ims_state)
            crude_time = crude_time + evolution_ims_state.period


        fs_binary_state, fs_ms = self.sampler.sample_fs_dist()[0]
           
        evolution_fs_state = MetastableState(binary_state=fs_binary_state,
                                             ms=fs_ms,
                                             dbs_state=self.triple_state,
                                             crude_time=crude_time,
                                             state=DBSEvolutionState.FS.value,
                                             aB=fs_binary_state.aB(),
                                             qB=fs_binary_state.qB())

        evolution.append(evolution_fs_state)
        logger.info(f"{Fore.BLUE}{crude_time}::fs_final_state::{evolution_fs_state}{Style.RESET_ALL}")
        evolution.append({"final": "unknown", "crude_time": crude_time})

        return evolution
