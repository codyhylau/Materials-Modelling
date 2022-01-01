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

e = np.linspace(-0.1,0.1,50) # strain (start,end,no. of steps)
PE = np.array([])
P = PE
vol = P
cell0 = cu.get_cell()

for i in e:
    cell = cell0*(1-i)
    cu.set_cell(cell, scale_atoms=True)
    vol = np.append(vol,cell[0][0]**3/cu.get_number_of_atoms()) # new volume to plot
    PE = np.append(PE,cu.get_potential_energy()/cu.get_number_of_atoms()) # to plot
    P = np.append(P,-np.trace(cu.get_stress(voigt=False))/3) # to plot

# find Bulk Modulus K 
# K = 130-140 GPa

# derivative method
strain = 0.01
cu = bulk("Cu", "fcc", a=3.6, cubic=True) 
cu.set_calculator(calc)
cell0 = cu.get_cell() # unstrained cell

cu2 = bulk("Cu", "fcc", a=3.6, cubic=True)  # make strained cell
cu2.set_calculator(calc)
cell = cu2.get_cell()
cell *= 1-strain # strain the cell
cu2.set_cell(cell, scale_atoms=True)

v = cell0[0,0]**3
K = v * (np.trace(cu2.get_stress(voigt=False))/3 - np.trace(cu.get_stress(voigt=False))/3 ) / (cell[0,0]**3 - v)  # taylor finite difference to find dP/dV
K /= GPa # convert eV/Ang^3 to GPa 
print('\nMy Bulk modulus: ', K, 'GPa')


from ase.eos import EquationOfState
from ase.units import kJ

# ASE method
eos = EquationOfState(vol, PE, eos="birchmurnaghan")
v0, e0, B = eos.fit()
print('ASE''s Bulk modulus: ', B / kJ * 1.0e24, 'GPa\n')




plt.grid()
plt.plot(vol,PE)
plt.xlabel('Volume / Ang^3')
plt.ylabel('PE / eV')
plt.title('Volume vs PE')
plt.figure(2)
plt.grid()
plt.plot(vol,P,'-r')
plt.xlabel('Volume / Ang^3')
plt.ylabel('Pressure / eV/Ang^3')
plt.title('Volume vs Pressure')
plt.show()




