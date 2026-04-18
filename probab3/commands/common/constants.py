from enum import Enum
import numpy as np
from matplotlib import pyplot as plt

''' Project Constants '''
PROJECT_NAME = 'probab3'

''' Scientific Constants '''
G = 6.674*(10**(-11))
HUBBLE_TIME = 4.55*(10**(17))
AU = 1.49598*(10**(11))
MSun = 1.98892*(10**(30))
c = 299792458
parsec = 3.086*(10**16)
RSun = 696342000
sinyear = 3.154*(10**7)

F_LIGO = 10


''' Project Enums'''
class MetastableFinalState(Enum):
    IN_CLUSTER_BINARY = "in_cluster_binary"
    IONIZED_BINARY = "ionized_binary"
    IMS_MERGER = "ims_merger"
    IMS_EM_MERGER = "ims_em_merger"
    FS_MERGER = "fs_merger"
    FS_EM_MERGER = "fs_em_merger"
    EJECTED_FS_MERGER = "ej_fs_merger"
    EJECTED_FS_EM = "ej_fs_em"
    EJECTED_BINARY = "ej_binary"
    IMS_TDE = "ims_tde"
    IMS_COLLISION = "ims_collision"
    FS_TDE = "fs_tde"
    FS_COLLISION = "fs_collision"
    EJECTED_FS_TDE = "ej_fs_tde"
    EJECTED_FS_COLLISION = "ej_fs_collision"
    TRIPLE_TD = "triple_td"
    IMS_BINARY_TD = "ims_binary_td"
    FS_BINARY_TD = "fs_binary_td"
    IONIZED_TRIPLE = "ionized_triple"


class EMSensitivityParams(Enum):
    e_min = "e_min"
    f_min = "f_min"

class MassIdentity(Enum):
    STAR = "STAR"
    BH = "BH" # Black Hole
    NS = "NS" # Neutron Star
    WD = "WD" # White Dwarf
    POINT = "POINT"

class DBSEvolutionState(Enum):
    FS = "FS"
    IMS = "IMS"
    EJECTED_FS = "EJECTED_FS"
    NEW_MS = "NEW_MS"
    TRIPLE_TD = "TRIPLE_TD"
    MANUAL_EXCHANGE = "MANUAL_EXCHANGE"

ALL_FINAL_STATES = {MetastableFinalState.IN_CLUSTER_BINARY.value,
                    MetastableFinalState.IONIZED_BINARY.value,
                    MetastableFinalState.IMS_MERGER.value,
                    MetastableFinalState.IMS_EM_MERGER.value,
                    MetastableFinalState.FS_MERGER.value,
                    MetastableFinalState.FS_EM_MERGER.value,
                    MetastableFinalState.EJECTED_FS_MERGER.value,
                    MetastableFinalState.EJECTED_FS_EM.value,
                    MetastableFinalState.EJECTED_BINARY.value,
                    MetastableFinalState.IMS_TDE.value,
                    MetastableFinalState.IMS_COLLISION.value,
                    MetastableFinalState.FS_TDE.value,
                    MetastableFinalState.FS_COLLISION.value,
                    MetastableFinalState.EJECTED_FS_TDE.value,
                    MetastableFinalState.EJECTED_FS_COLLISION.value,
                    MetastableFinalState.TRIPLE_TD.value,
                    MetastableFinalState.IMS_BINARY_TD.value,
                    MetastableFinalState.FS_BINARY_TD.value,
                    MetastableFinalState.IONIZED_TRIPLE.value}

MERGER_ONLY_FINAL_STATES = {MetastableFinalState.IMS_EM_MERGER.value,
                            MetastableFinalState.IMS_MERGER.value,
                            MetastableFinalState.FS_EM_MERGER.value,
                            MetastableFinalState.FS_MERGER.value,
                            MetastableFinalState.EJECTED_FS_EM.value,
                            MetastableFinalState.EJECTED_FS_MERGER.value}

BINARY_FINAL_STATES = {MetastableFinalState.IN_CLUSTER_BINARY.value,
                       MetastableFinalState.IONIZED_BINARY.value,
                       MetastableFinalState.EJECTED_BINARY.value,
                       MetastableFinalState.IMS_BINARY_TD.value}

EQUAL_MASS_FINAL_STATES = {MetastableFinalState.IMS_EM_MERGER.value,
                            MetastableFinalState.IMS_MERGER.value,
                            MetastableFinalState.FS_EM_MERGER.value,
                            MetastableFinalState.FS_MERGER.value,
                            MetastableFinalState.EJECTED_FS_EM.value,
                            MetastableFinalState.EJECTED_FS_MERGER.value,
                            MetastableFinalState.IN_CLUSTER_BINARY.value,
                            MetastableFinalState.IONIZED_BINARY.value,
                            MetastableFinalState.EJECTED_BINARY.value}


OBJECT_IDENTITIES = [
    MassIdentity.STAR.value,
    MassIdentity.WD.value,
    MassIdentity.NS.value,
    MassIdentity.BH.value
]

''' Colors '''

COLORS_PLATE  =  plt.cm.tab20(np.linspace(0, 1, len(ALL_FINAL_STATES) - 4))
COLORS_PLATE2  =  plt.cm.tab20c(np.linspace(0, 1, 20))

FS_COLORS = {
    MetastableFinalState.IN_CLUSTER_BINARY.value: COLORS_PLATE[0],
    MetastableFinalState.IONIZED_TRIPLE.value: COLORS_PLATE[1],
    MetastableFinalState.IMS_MERGER.value: COLORS_PLATE[2],
    MetastableFinalState.IMS_EM_MERGER.value: COLORS_PLATE[3],
    MetastableFinalState.FS_MERGER.value: COLORS_PLATE[4],
    MetastableFinalState.FS_EM_MERGER.value: COLORS_PLATE[5],
    MetastableFinalState.EJECTED_FS_MERGER.value: COLORS_PLATE[6],
    MetastableFinalState.EJECTED_FS_EM.value: COLORS_PLATE[7],
    MetastableFinalState.EJECTED_BINARY.value: COLORS_PLATE[8],
    MetastableFinalState.IMS_TDE.value: COLORS_PLATE[9],
    MetastableFinalState.IMS_COLLISION.value: COLORS_PLATE[10],
    MetastableFinalState.FS_TDE.value: COLORS_PLATE[11],
    MetastableFinalState.FS_COLLISION.value: COLORS_PLATE[12],
    MetastableFinalState.EJECTED_FS_TDE.value: COLORS_PLATE[13],
    MetastableFinalState.EJECTED_FS_COLLISION.value: COLORS_PLATE[14],
    MetastableFinalState.TRIPLE_TD.value: COLORS_PLATE2[-3],
    MetastableFinalState.IMS_BINARY_TD.value: COLORS_PLATE2[-1],
    MetastableFinalState.FS_BINARY_TD.value: COLORS_PLATE2[-2]
}

FS_LABELS = {
    MetastableFinalState.IN_CLUSTER_BINARY.value: "In Cluster Binary",
    MetastableFinalState.IONIZED_BINARY.value: "Ionized Binary",
    MetastableFinalState.IMS_MERGER.value: "IMS Merger",
    MetastableFinalState.IMS_EM_MERGER.value: "IMS OEM",
    MetastableFinalState.FS_MERGER.value: "FS Merger",
    MetastableFinalState.FS_EM_MERGER.value: "FS OEM",
    MetastableFinalState.EJECTED_FS_MERGER.value: "Ejected Merger",
    MetastableFinalState.EJECTED_FS_EM.value: "Ejected OEM",
    MetastableFinalState.EJECTED_BINARY.value: "Ejected Binary",
    MetastableFinalState.IMS_TDE.value: "IMS TDE",
    MetastableFinalState.IMS_COLLISION.value: "IMS Collision",
    MetastableFinalState.FS_TDE.value: "FS TDE",
    MetastableFinalState.FS_COLLISION.value: "FS Collision",
    MetastableFinalState.EJECTED_FS_TDE.value: "Ejected TDE",
    MetastableFinalState.EJECTED_FS_COLLISION.value: "Ejected Collision",
    MetastableFinalState.TRIPLE_TD.value: "Triple TD",
    MetastableFinalState.IMS_BINARY_TD.value: "IMS Binary TD",
    MetastableFinalState.FS_BINARY_TD.value: "FS Binary TD",
    MetastableFinalState.IONIZED_TRIPLE.value: "Ionized Triple"
}
