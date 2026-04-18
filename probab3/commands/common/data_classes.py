import numpy as np
import dataclasses
from typing import Optional

from probab3.commands.common.formulas.general import *

@dataclasses.dataclass
class MassObj():
    mass: int
    identity: MassIdentity

    def radius(self):
        if self.identity == MassIdentity.STAR.value:
            return star_radius(self.mass)
        elif self.identity == MassIdentity.WD.value:
            return WD_radius(self.mass)
        else: # BHs or NSs
            return r_schwarzschild(self.mass)

@dataclasses.dataclass
class TripleState():
    ma: MassObj
    mb: MassObj
    ms: MassObj
    E0: int
    L0: int
    EB: Optional[int] = None
    LB: Optional[int] = None
    CB: Optional[float] = None
    C0: Optional[float] = None

@dataclasses.dataclass
class BinaryState():
    ma: MassObj
    mb: MassObj
    EB: int
    LB: int
    CB: float
    CBd: float

    def aB(self):
        return aB(self.ma.mass, self.mb.mass, self.EB)
    
    def eB(self):
        return eB_from_EBLB(self.ma.mass, self.mb.mass, self.EB, self.LB)

    def qB(self):
        return qB_from_EBLB(self.ma.mass, self.mb.mass, self.EB, self.LB)
    
    def r_TDE(self):
        small_mass = self.mb
        big_mass = self.ma
        # Check if one of the masses is a star
        if self.ma.identity == MassIdentity.STAR.value:
            if self.mb.identity != MassIdentity.STAR.value or self.mb.mass > self.ma.mass:
                small_mass = self.ma
                big_mass = self.mb
        # If not, check if one of the masses is a WD
        elif MassIdentity.STAR.value not in {small_mass.identity, big_mass.identity}:
            if self.ma.identity == MassIdentity.WD.value:
                if self.mb.identity != MassIdentity.WD.value or self.mb.mass > self.ma.mass:
                    small_mass = self.ma
                    big_mass = self.mb
        # If not, check if one of the masses is a NS
        elif MassIdentity.WD.value not in {small_mass.identity, big_mass.identity} and \
                MassIdentity.STAR.value not in {small_mass.identity, big_mass.identity}: 
            if self.ma.identity == MassIdentity.NS.value:
                if self.mb.identity != MassIdentity.NS.value or self.mb.mass > self.ma.mass:
                    small_mass = self.ma
                    big_mass = self.mb
            # Check if both BHs
            elif self.mb.identity == MassIdentity.BH.value and self.mb.mass > self.ma.mass:
                small_mass = self.ma
                big_mass = self.mb

        return r_TDE(small_mass=small_mass.mass, small_radius=small_mass.radius(), big_mass=big_mass.mass)

    def r_collision(self):
        return self.ma.radius() + self.mb.radius()

@dataclasses.dataclass
class TeritiaryState():
    ms: MassObj
    Es: int
    Ls: int
    CBs: float
    Cs: float

@dataclasses.dataclass
class MetastableState():
    dbs_state: Optional[TripleState] = None
    binary_state: Optional[BinaryState] = None
    teritiary_state: Optional[TeritiaryState] = None
    crude_time: Optional[int] = None
    state: Optional[DBSEvolutionState] = None
    ms: Optional[MassObj] = None
    r_TDE: Optional[int] = None
    r_collision: Optional[int] = None
    qB: Optional[int] = None
    aB: Optional[int] = None
    qB_merger: Optional[int] = None
    qB_EM: Optional[int] = None
    period: Optional[int] = None
    scramble_num: Optional[int] = None
    dbs_kick: Optional[int] = None
    ejection_energy: Optional[int] = None
    EB_merger: Optional[int] = None
    EB_EM: Optional[int] = None
    P_dis: Optional[float] = None
    num_of_scrambles: Optional[int] = None
    a_s: Optional[int] = None
    a_s_cluster_tides: Optional[int] = None
    EB_cluster_tides: Optional[int] = None
