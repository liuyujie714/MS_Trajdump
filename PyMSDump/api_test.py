from PyMSDump.MSTrj import MSTrjReader, EnergyType, PressVolType
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import sys

def test_energy(ftrj:str, fpdb:str):
    """ @brif Test all kinds of energies """
    trj = MSTrjReader(ftrj, fpdb)
    for ts in trj:
        print([ts.ener[e] for e in EnergyType])

def test_pressvol(ftrj:str, fpdb:str):
    """ @brif Test all kinds of pressure and volume """
    trj = MSTrjReader(ftrj, fpdb)
    for ts in trj:
        print([ts.pvol[p] for p in PressVolType])

def test_vector(ftrj:str, fpdb:str):
    """ @brif Test all kinds of vectors, such as positions, velocities, forces """
    trj = MSTrjReader(ftrj, fpdb)
    for ts in trj:
        print(ts.positions)
        if ts.hasV:
            print(ts.velocities)
        if ts.hasF:
            print(ts.forces)
            
def test_write_mda(ftrj:str, fpdb:str, fout:str):
    """ @brief Test writing MDAnalysis trajectory from MSTrjReader """
    trj = MSTrjReader(ftrj, fpdb)
    u = mda.Universe.empty(n_atoms=trj.natoms, trajectory=True)
    with XTCWriter(fout, trj.natoms) as w:
        for ts in trj:
            u.atoms.positions = ts.positions
            u.dimensions = ts.crystal
            u.trajectory.ts.frame = ts.step
            u.trajectory.ts.time = ts.time
            w.write(u.atoms)

if __name__ == "__main__":
    test_write_mda(sys.argv[1], sys.argv[2], "test.xtc")
