from probab3.commands.sample_dbs.dbs_sampler import * 
from probab3.commands.sample_dbs.sample_and_save import *

logging.basicConfig(level=logging.INFO, filename=LOG_FILE_NAME, filemode='a', format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)
colorama_init()

class SampleDBSUnequalMassTriple(ClusterBinaryTeritiarySampler):
    """A class to sample an un-equal mass dynamical binary sequence system in a star cluster."""

    def __init__(self, alpha: float, Mc: int, rh: int, rc: int,  
                 rejection: bool=True, initial_binary: BinaryState=None,):
        """Initialize SampleDBSUnequalMassTriple object. Use SI units.
        
        Args:
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            rejection (bool, optional): Whether to reject binaries that are too hard or too soft. Defaults to True.
            rejection (bool): Whether to use rejection sampling.  
        """
        super(SampleDBSUnequalMassTriple, self).__init__(alpha, Mc, rh, rc, 
                                                        rejection,
                                                        sampler_class=UnequalMassTripleSampler)
        if initial_binary:
            self.binary_init_state = initial_binary
        else:
            ma = self.cluster.sample_star_mass()
            mb = self.cluster.sample_star_mass() 
            self.init_stars_binary(ma, mb)
        
        self.binary_state = self.binary_init_state

    def init_cluster(self, Mc, rh, rc):
        """Initialize the cluster."""
        self.cluster = star_cluster.UnequalMassCluster(Mc_stars=Mc, rh_stars=rh, rc_stars=rc)
        self.cluster.update_pdmf() 

    def init_stars_binary(self, ma, mb):
        """Initialize the inital binary state with the given mass objects.
        
        Assume thermal eccentricity distribution and Opik's law for the initial semi-major axis distribution.

        Args:
            ma (MassObj): The mass of the first object.
            mb (MassObj): The mass of the second object.
        """
        self.e_in = self.cluster.stars_cluster.sample_e_in_thermal()
        self.a_in = self.cluster.stars_cluster.sample_a_in_main_sequence(self.e_in)
        CB_sampled, _ = self.sample_CB() 
        self.binary_init_state = BinaryState(ma=ma, 
                                             mb=mb, 
                                             EB=EB(ma=ma.mass, mb=mb.mass, aB=self.a_in), 
                                             LB=LB_from_aBeB(ma=ma.mass, mb=mb.mass, aB=self.a_in, eB=self.e_in), 
                                             CB=CB_sampled,
                                             CBd=CB_sampled)
        
    def sample_ms(self):
        """Sample the mass of the new tertiary. Use a present day mass function (PDMF) to do so."""
        return self.cluster.sample_mass()

    def __str__(self):
        return f"Sampling DBS with alpha={self.alpha}, cluster={self.cluster}, binary_init_state={self.binary_init_state}"

class StarsBinaryToBHsBinaryUnequalMass(StarsBinaryToBHsBinary):
    """A class to sample a star binary and a teritiary BH from a star cluster and evolve it to a binary of BHs."""

    def __init__(self, Mc: int, rh: int, rc: int, alpha: float = 2.5, rejection: bool=True):
        """Initialize the StarsBinaryToBHsBinaryUnequalMass.
        
        Args:
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            alpha (float): The alpha parameter for the chaotic triple range.
            rejection (bool): Whether to use rejection sampling.
        """
        super(StarsBinaryToBHsBinaryUnequalMass, self).__init__(Mc, rh, rc, None, None, alpha, rejection)
    
    def init_cluster(self, Mc, rh, rc):
        """Initialize the cluster."""
        self.cluster = star_cluster.UnequalMassCluster(Mc_stars=Mc, rh_stars=rh, rc_stars=rc)
        self.cluster.update_pdmf()
    
    def get_star_mass(self):
        """Sample star mass from cluster pdmf."""
        return self.cluster.sample_star_mass()
    
    def get_BH_mass(self):
        """Sample BH mass from cluster pdmf."""
        return self.cluster.sample_BH_mass()

class SampleMSs2BHs2DBSUnequalMass:
    """A class to evolve an initially sampled stars binary to a BHs binary and then to a dynamical binary sequence system. 
    Sample all masses from cluster's pdmf.
    """

    def __init__(self, alpha: float, Mc: int, rh: int, rc: int, 
                 rejection: bool=True):
        """Initially sample the stars binary, then evolve it to a BHs binary and construct a dynamical binary sequence system.
        
        Args:
            alpha (float): The alpha parameter for the chaotic triple range.
            Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
            rh (int): Half mass radius of the stars cluster.
            rc (int): Core radius of the stars cluster.
            star_mass (int): The mass of the stars in the cluster.
            BH_mass (int): The mass of the BHs in the cluster.
        """ 
        self.MSs_to_BHs_sampler = StarsBinaryToBHsBinaryUnequalMass(Mc, rh, rc, rejection) 
        MSs2BHs_evolution, star_BH_state, pre_evolution_fails = self.MSs_to_BHs_sampler.evolve_stars_to_BHs()
        self.MSs_binary_init_state = MSs2BHs_evolution[0].binary_state
        self.BHs_binary_init_state = MSs2BHs_evolution[-1].binary_state
        self.pre_evolution = {"MSs2BHs_evolution": MSs2BHs_evolution,
                              "star_BH_state": star_BH_state,
                              "pre_evolution_fails": pre_evolution_fails}

        self.dbs_sampler = SampleDBSUnequalMassTriple(alpha=alpha, Mc=Mc, rh=rh, rc=rc, 
                                                    initial_binary=self.BHs_binary_init_state, 
                                                    rejection=rejection)

    def evolve(self):
        """Evolve the dynamical binary sequence system until its end.
        
        Return:
            (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]) 
                The resolution and the evolution list with 
                    {"final": resolution, "crude_time": crude time} appended at the end.
        """
        return self.dbs_sampler.evolve_binary()


@sample_and_save
def unequal_mass_sample(alpha: float, Mc: int, rh: int, rc: int, rejection: bool=True, 
                        size: int=1, output_file_path: str=""):
    """Sample an unequal mass triple system N times and save the sample to a file.
    
    Args:
        alpha (float): The alpha parameter for the chaotic triple range.
        Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
        rh (int): Half mass radius of the stars cluster.
        rc (int): Core radius of the stars cluster.
        rejection (bool, optional): Whether to reject binaries that are too hard or too soft. Defaults to True.
        size (int, optional): The number of samples to generate. Defaults to 1.
        output_file_path (str, optional): The path to the output file. Defaults to "".

    Return:
        {str: (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]),
         str: MetastableFinalState}
        A dictionary with the evolution and the end state of the format {"evolution": evolution, "end": end}
    """
    dbs_sample_obj = SampleDBSUnequalMassTriple(alpha, Mc, rh, rc, rejection)
    end, evolution = dbs_sample_obj.evolve_binary()
    return {
        "evolution": evolution,
        "end": end
    }
    

@sample_and_save
def unequal_mass_bbh_sample(alpha: float, Mc: int, rh: int, rc: int, rejection: bool=True, 
                            size: int=1, output_file_path: str=""):
    """Evolve a stars binary to a BHs binary and then evolve a triple system in the cluster.
    Do this N times and save the sample to a file. Use SI units. Sample masses from cluster's pdmf.
    
    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file.

    Args:
        alpha (float): The alpha parameter for the chaotic triple range.
        Mc (int): Cluster mass of stars, BH subcluster properties will be calculated via TwoComponentCluster.
        rh (int): Half mass radius of the stars cluster.
        rc (int): Core radius of the stars cluster.
        rejection (bool): Whether to use rejection sampling.
        size (int): The number of evolutions to generate.
        output_file_path (str): The path to the output file.

    Return:
        {str: (MetastableFinalState, [MetastableState, ..., MetastableState, {str: MetastableFinalState, str: int}]),
         str: MetastableFinalState}
        A dictionary with the evolution and the end state of the format {"evolution": evolution, "end": end}
    """
    dbs_sample_obj = SampleMSs2BHs2DBSUnequalMass(alpha, Mc, rh, rc, rejection)
    pre_evolution = dbs_sample_obj.pre_evolution

    logger.info(f"{Fore.MAGENTA}Pre-evolution finished, starting mt evolution {Style.RESET_ALL}")

    end, evolution = dbs_sample_obj.evolve()
    return {
        "pre_evolution": pre_evolution,
        "evolution": evolution,
        "end": end
    }


@sample_and_save
def one_triple_unequal_mass_sample(m1, m2, m3, alpha, E0, L0, output_file_path, size=1, rejection=True):
    """Sample a triple system N times and save the sample to a file. Use SI units.

    Note: This function uses the sample_and_save decorator to repeat and save the samples to the file. 
          The decorator is defined in probab3/commands/sample_dbs/sample_and_save.py.

    Args:
        m1 (int): The mass of the 1st object.
        m2 (int): The mass of the 2nd object.
        m3 (int): The mass of the 3rd object.
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
                                      mb=MassObj(mass=m2, identity=MassIdentity.POINT.value),
                                      ms=MassObj(mass=m3, identity=MassIdentity.POINT.value),
                                      E0=E0, L0=L0,
                                      alpha=alpha, rejection=rejection)
    
    evolution = dbs_sample_obj.evolve_triple()

    return {
        "evolution": evolution
    }