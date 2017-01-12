from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

def get_energy_ads(x,y,h):

    slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0)

    slab.set_calculator(EMT())
    e_slab = slab.get_potential_energy()

    amol = Atoms('N', positions=[(x, y, 0)])
    amol.set_calculator(EMT())
    e_N = amol.get_potential_energy()

    add_adsorbate(slab, amol, h, 'ontop')
    constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
    slab.set_constraint(constraint)
    #dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
    #dyn.run(fmax=0.05)
    e_N_slab = slab.get_potential_energy()

    return  -e_slab - e_N +e_N_slab

if __name__ == '__main__':
    x=0.1
    y=0.5
    h=2.7
    print get_energy_ads(x,y,h)
