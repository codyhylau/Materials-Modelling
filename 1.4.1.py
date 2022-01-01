from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa
import Morse
from matplotlib import pyplot as plt
from ase.build import bulk
calc = Morse.MorsePotential()
cu = bulk("Cu", "fcc", a=3.6, cubic=True) 
cu.set_calculator(calc)

# Find Shear Modulus G
strain = 0.01
cell0 = cu.get_cell()
cell0[0][1] = cell0[0][0] * strain
cu.set_cell(cell0, scale_atoms=True)
s = cu.get_stress(voigt=False)

tau1 = (s[0][1]+s[1][0])/2 #top mid and left mid
tau2 = np.trace(s)/3 # diagonal
tau3 = (tau1 + tau2)/2 # all non zero
print(s)
# G of copper = 45 GPa
G = tau1/strain
G /= GPa # convert eV/Ang^3 to GPa 
print('\nCalculated Shear Modulus: ', G, 'GPa\n')
