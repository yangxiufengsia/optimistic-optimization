from math import *
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
import random


def get_energy_ads(x,y,h):

    slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0)

    slab.set_calculator(EMT())
    e_slab = slab.get_potential_energy()

    amol = Atoms('N', positions=[(x, y, 0)])
    amol.set_calculator(EMT())
    e_N = amol.get_potential_energy()

    add_adsorbate(slab, amol, h, position=(x,y))
    constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
    slab.set_constraint(constraint)
    #dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
    #dyn.run(fmax=0.05)
    e_N_slab = slab.get_potential_energy()

    return  -(-e_slab - e_N +e_N_slab)
random=[]
#for i in range(500):
random_ini=[]

x1=np.random.uniform(0.0,3.0,5000)
x2=np.random.uniform(0.0,3.0,5000)
x3=np.random.uniform(1.0,4.0,5000)
newx1=np.random.permutation(x1)
newx2=np.random.permutation(x2)
newx3=np.random.permutation(x3)
print x1
for i in range(len(x1)):
    random_ini.append(get_energy_ads(newx1[i],newx2[i],newx3[i]))
random_index=np.argmax(random_ini)
random.append(random_ini[random_index])
print newx1[random_index]
print newx2[random_index]
print newx3[random_index]
print random
print random_ini
