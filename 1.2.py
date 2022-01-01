from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa
import Morse
from matplotlib import pyplot as plt
calc = Morse.MorsePotential()
d = 3*Ang

e = 1e-12 #finite difference
a1 = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d)])
a1.set_calculator(calc) 
a2 = Atoms('2Cu', positions=[(0., 0., 0.), (0., 0., d+e)])
a2.set_calculator(calc) 
f = np.linalg.norm(a1.get_forces()[0])
f2 = -(a1.get_potential_energy()-a2.get_potential_energy())/e #taylor

print('\nForce from function: ' + str(f))
print('Force from derivative: ' + str(f2) +'\n')
print(a1.get_forces())

#if e is too small, it becomes unstable and reaches 0