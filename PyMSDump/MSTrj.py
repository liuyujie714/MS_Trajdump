from PyMSDump_ import TrajLoad
from enum import IntEnum, auto, unique
from dataclasses import dataclass
from typing import Iterator, Union
import numpy as np
import sys

"""
 * Materials Studio .trj use below units:
 * All times are in ps.
 * All energies are in kcal mol-1.
 * All pressure and stress values are in GPa.
 * All volumes are in Å3.
 * All coordinates are in Å.
 * All velocities are in Å ps-1.
 * All forces are in kcal mol-1 Å-1.
 * All temperatures are in K.
 * All logical values are stored as integers (0=FALSE, not 0=TRUE).
"""

@dataclass
class Frame:
    """ A frame structure in a trajectory
    Attributes:
        step (int): current step
        time (float): current time, ps
        positions (np.ndarray): Atom positions in Å
        velocities (Union[np.ndarray, None]): Atom velocities in Å ps-1
        forces (Union[np.ndarray, None]): Atoms forces in kcal mol-1 Å-1
        box (np.ndarray): The converted lower triangular simulation box (same as .gro but unit is Angstrom)。
        crystal (np.ndarray): pdb crystal, unit is angle and angstrom
        ener (np.ndarray): energy information
        pvol (np.ndarray): pressure and volume information
        hasV (bool): if has velocity
        hasF (bool): if has force
    """
    step: int
    time: float
    positions: np.ndarray
    velocities: Union[np.ndarray, None]
    forces: Union[np.ndarray, None]
    box: np.ndarray
    crystal: np.ndarray
    ener: np.ndarray
    pvol: np.ndarray
    hasV: bool
    hasF: bool

@unique
class EnergyType(IntEnum):
    Temp=0  # must from 0
    AvgTemp=auto()
    TimeStep=auto()
    InitialTemp=auto()
    FinalTemp=auto()
    TotalPE=auto()
    BondE=auto()
    AngleE=auto()
    TorsionE=auto()
    InversionE=auto()
    vdWE=auto()
    ElectrostaticE=auto()
    HBondE=auto()
    ConstraintE=auto()
    UreyBradleyE=auto()
    ThreeBodyE=auto()
    TotalCrossTermE=auto()
    BendBendE=auto()
    StretchStretchE=auto()
    StretchBendStretchE=auto()
    StretchTorsionStretchE=auto()
    BendTorsionBendE=auto()
    TorsionBendBendE=auto()
    SeperatedStretchStretchE=auto()
    TorsionStretchE=auto()
    InversionInversionE=auto()
    UserE=auto()
    TotalInternalE=auto()
    TotalNonBondE=auto()
    AvgTotalPE=auto()
    AvgBondE=auto()
    AvgAngleE=auto()
    AvgTorsionE=auto()
    AvgInversionE=auto()
    AvgvdWE=auto()
    AvgElectrostaticE=auto()
    AvgHBondE=auto()
    AvgConstraintE=auto()
    AvgUreyBradleyE=auto()
    AvgThreeBodyE=auto()
    AvgTotalCrossTermE=auto()
    AvgBendBendE=auto()
    AvgStretchStretchE=auto()
    AvgStretchBendStretchE=auto()
    AvgStretchTorsionStretchE=auto()
    AvgBendTorsionBendE=auto()
    AvgTorsionBendBendE=auto()
    AvgSeperatedStretchStretchE=auto()
    AvgTorsionStretchE=auto()
    AvgInversionInversionE=auto()
    AvgUserE=auto()
    AvgTotalInternalE=auto()
    AvgTotalNonBondE=auto()
    TotalE=auto()
    TotalKE=auto()
    AvgTotalE=auto()
    AvgTotalKE=auto()

@unique
class PressVolType(IntEnum):
    Press = 0   # must from 0
    Volume=auto()
    TotalPV=auto()
    KineticStrsPV=auto()
    PotentialStrsPV=auto()
    GyrationRadius=auto()
    AvgPress=auto()
    AvgVolume=auto()
    AvgTotalPV=auto()
    AvgKineticStrsPV=auto()
    AvgPotentialStrsPV=auto()
    AvgGyrationRadius=auto()

class MSTrjReader:
    def __init__(self, ftrj:str, fpdb:str):
        self.trajectory = TrajLoad(ftrj, fpdb)
        self.natoms_ = self.__natoms()
        self.nframes_ = self.__len__()
    def __natoms(self):
        nat = 0
        for fr in self:
            nat = len(fr.positions)
            break
        return  nat
    def __iter__(self) -> Iterator[Frame]:
        try:
            for _ in self.trajectory:
                yield _
        finally:
            # always reset pointer to header, such for and break
            self.trajectory.reset()
    def __len__(self):
        return len([_ for _ in self.trajectory])
    def __str__(self):
        return  f'Total frames: {len(self)}'
    @property  
    def nframes(self):
        """ @brief Total number of frames """
        return self.nframes_
    @property    
    def natoms(self):
        """ @brief Number of atoms in the trajectory """
        return self.natoms_

if __name__ == '__main__':
    trj = MSTrjReader(sys.argv[1], sys.argv[2])
    for ts in trj:
        print(ts)
