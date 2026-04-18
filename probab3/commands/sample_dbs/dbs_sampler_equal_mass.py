from probab3.commands.sample_dbs.dbs_sampler import * 
from probab3.commands.sample_dbs.sample_and_save import *

logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)
colorama_init()


class SampleDynamicalBinarySequence(ClusterBinaryTeritiarySampler):
    """A class to sample an equal mass dynamical binary sequence system in a star cluster."""

    def __init__(self, ma: MassObj, mb: MassObj, ms: Optional[MassObj], alpha: float, Mc: int, rh: int, rc: int, initial_binary: BinaryState=None, 
                 rejection: bool=True):
        """Initialize the sampler.
        
        Args:
            ma (MassObj): Mass object of the binary.
            mb (MassObj): Mass object of the binary.
            ms (MassObj): Optional. Initial mass of the teritiary.
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            rejection (bool): Whether to use rejection sampling.            
        """ 
        super(SampleDynamicalBinarySequence, self).__init__(alpha, Mc, rh, rc, 
                                                     rejection, 
                                                     EqualMassTripleSampler)
        self.ms = ms
        if initial_binary:
            self.a_in = initial_binary.aB()
            self.binary_init_state = initial_binary
        else:
            # Above the hard binary limit, binaries tend to seperate, below they tend to become harder. 
            self.a_in = self.cluster.aHB(ma.mass, mb.mass, ms.mass) # k=1.47) # fix later
            CB_sampled, _= self.sample_CB() 
            self.binary_init_state = BinaryState(ma=ma, mb=mb, EB=EB(ma=ma.mass, mb=mb.mass, aB=self.a_in), 
                                                 LB=LB_from_aBeB(ma=ma.mass, mb=mb.mass, aB=self.a_in, eB=0), # default is eccentricity 0 
                                                 CB=CB_sampled, CBd=CB_sampled)

        self.cluster_ejection_energy = self.cluster.ejection_energy(ma=ma.mass, mb=mb.mass, ms=ms.mass) 

    def __str__(self):
        return f"Sampling DBS with alpha={self.alpha}, cluster={self.cluster}, binary_init_state={self.binary_init_state}"

    def get_ejection_energy(self, fs_binary_state, fs_ms):
        """Return the ejecton energy of the first calculated sustem. This is done to reduce calculations."""
        return self.cluster_ejection_energy


class SampleMSs2BHs2DBS:
    """A class to evolve an initially sampled stars binary to a BHs binary and then to a dynamical binary sequence system."""

    def __init__(self, alpha: float, Mc: int, rh: int, rc: int, 
                 rejection: bool=True, star_mass: int=MSun, BH_mass: int=20*MSun):
        """Initially sample the stars binary, then evolve it to a BHs binary and construct a dynamical binary sequence system.
        
        Args:
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            star_mass (int): The mass of the stars in the cluster.
            BH_mass (int): The mass of the BHs in the cluster.
        """ 
        self.MSs_to_BHs_sampler = StarsBinaryToBHsBinary(Mc, rh, rc, star_mass, BH_mass, alpha, rejection) 
        MSs2BHs_evolution, star_BH_state, pre_evolution_fails = self.MSs_to_BHs_sampler.evolve_stars_to_BHs()
        self.MSs_binary_init_state = MSs2BHs_evolution[0].binary_state
        self.BHs_binary_init_state = MSs2BHs_evolution[-1].binary_state
        self.pre_evolution = {"MSs2BHs_evolution": MSs2BHs_evolution,
                              "star_BH_state": star_BH_state,
                              "pre_evolution_fails": pre_evolution_fails}

        self.dbs_sampler = SampleDynamicalBinarySequence(MassObj(mass=BH_mass, identity=MassIdentity.BH.value), 
                                                 MassObj(mass=BH_mass, identity=MassIdentity.BH.value), 
                                                 MassObj(mass=BH_mass, identity=MassIdentity.BH.value), 
                                                 alpha, Mc, rh, rc, 
                                                 self.BHs_binary_init_state, rejection)

    def evolve(self):
        """Evolve the dynamical binary sequence system until its end.
        
        Return:
            (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]) 
                The resolution and the evolution list with 
                    {"final": resolution, "crude_time": crude time} appended at the end.
        """
        return self.dbs_sampler.evolve_binary()


@sample_and_save
def basic_equal_mass_sample(m1, m2, m3, alpha, Mc, rh, rc, output_file_path, size=1, rejection=True):
    """Sample a triple system N times and save the sample to a file. Use SI units.

    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file. 
          The decorator is defined in probab3/commands/sample_dbs/sample_and_save.py.

    Args:
        m1 (int): The mass of the first object.
        m2 (int): The mass of the second object.
        m3 (int): The mass of the third object.
        alpha (float): The alpha parameter for the chaotic triple range.
        Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
        rh (int): Half mass radius of the stars cluster.
        rc (int): Core radius of the stars cluster.
        output_file_path (str): The path to the output file.
        size (int): The number of evolutions to generate.
        rejection (bool): Whether to use rejection sampling.
    
    Return:
        {str: (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]),
         str: MetastableFinalState}
        A dictionary with the evolution and the end state of the format {"evolution": evolution, "end": end}
    """ 
    dbs_sample_obj = SampleDynamicalBinarySequence(ma=MassObj(mass=m1, identity=MassIdentity.BH.value),
                                           mb=MassObj(mass=m2, identity=MassIdentity.BH.value),
                                           ms=MassObj(mass=m3, identity=MassIdentity.BH.value),
                                           alpha=alpha, Mc=Mc, rh=rh, rc=rc, rejection=rejection)
    
    end, evolution = dbs_sample_obj.evolve_binary()

    return {
        "evolution": evolution,
        "end": end
    }

@sample_and_save
def stars_to_BHs_sample(star_mass, BH_mass, alpha, Mc, rh, rc, output_file_path, size=1, rejection=True):
    """Evolve a star binary to BH binary N times and save the sample to a file. Use SI units.
    
    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file.

    Args:
        star_mass (int): The average mass of a star in the cluster.
        BH_mass (int): The average mass of a BH in the cluster.
        alpha (float): The alpha parameter for the chaotic triple range.
        Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
        rh (int): Half mass radius of the stars cluster.
        rc (int): Core radius of the stars cluster.
        output_file_path (str): The path to the output file.
        size (int): The number of evolutions to generate.
        rejection (bool): Whether to use rejection sampling.

    Return:
        {str: (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]),
         str: MetastableFinalState}
        A dictionary with the evolution and the end state of the format {"evolution": evolution, "end": end}
    """

    ms_stars = StarsBinaryToBHsBinary(Mc, rh, rc, star_mass, BH_mass, alpha, rejection)
    evolution, star_BH_state, pre_ev_fails = ms_stars.evolve_stars_to_BHs()

    return {
        "evolution": evolution,
        "star_BH_state": star_BH_state,
        "pre_evolution_fails": pre_ev_fails
    }

@sample_and_save
def ms_to_BH_to_end_sample(alpha: float, Mc: int, rh: int, rc: int, rejection: bool=True, 
                              star_mass: int=MSun, BH_mass: int=20*MSun, size=1, 
                              output_file_path=""):
    """Evolve a stars binary to a BHs binary and then evolve a triple system in the cluster.
    Do this N times and save the sample to a file. Use SI units.
    
    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file.

    Args:
        alpha (float): The alpha parameter for the chaotic triple range.
        Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
        rh (int): Half mass radius of the stars cluster.
        rc (int): Core radius of the stars cluster.
        rejection (bool): Whether to use rejection sampling.
        star_mass (int): The average mass of a star in the cluster.
        BH_mass (int): The average mass of a BH in the cluster.
        size (int): The number of evolutions to generate.
        output_file_path (str): The path to the output file.

    Return:
        {str: (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]),
         str: MetastableFinalState}
        A dictionary with the evolution and the end state of the format {"evolution": evolution, "end": end}
    """
    dbs_sample_obj = SampleMSs2BHs2DBS(alpha, Mc, rh, rc, rejection, star_mass, BH_mass)
    pre_evolution = dbs_sample_obj.pre_evolution

    logger.info(f"{Fore.MAGENTA}Pre-evolution finished, starting mt evolution {Style.RESET_ALL}")

    end, evolution = dbs_sample_obj.evolve()
    return {
        "pre_evolution": pre_evolution,
        "evolution": evolution,
        "end": end
    }


@sample_and_save
def one_triple_equal_mass_sample(m1, alpha, E0, L0, output_file_path, size=1, rejection=True):
    """Sample a triple system N times and save the sample to a file. Use SI units.

    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file. 
          The decorator is defined in probab3/commands/sample_dbs/sample_and_save.py.

    Args:
        m1 (int): The mass of the any object.
        alpha (float): The alpha parameter for the chaotic triple range.
        E0 (int): The initial energy of the system.
        L0 (int): The initial angular momentum of the system.
        output_file_path (str): The path to the output file.
        size (int): The number of evolutions to generate.
        rejection (bool): Whether to use rejection sampling.
    
    Return:
        ({str: [MetastableState, ..., MetastableState]}):
        A dictionary with the evolution of the format {"evolution": evolution}
    """ 
    dbs_sample_obj = OneTripleSampler(ma=MassObj(mass=m1, identity=MassIdentity.POINT.value),
                                      mb=MassObj(mass=m1, identity=MassIdentity.POINT.value),
                                      ms=MassObj(mass=m1, identity=MassIdentity.POINT.value),
                                      E0=E0, L0=L0,
                                      alpha=alpha, rejection=rejection, sampler_class=EqualMassTripleSampler)
    
    evolution = dbs_sample_obj.evolve_triple()

    return {
        "evolution": evolution
    }