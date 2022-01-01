from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa
import Morse
from matplotlib import pyplot as plt
calc = Morse.MorsePotential()

start = 2.3
end = 5 # end distance
n = 100 # number of points to plot
oneton = np.linspace(start,end+start,n)
dist = np.zeros((n,3))
dist[:,2] = oneton
PE = np.array([])
F = PE

for d in range(n):
    a = Atoms('2Cu', positions=[(0., 0., 0.), tuple(dist[d])])
    a.set_calculator(calc)
    pe = a.get_potential_energy()
    PE = np.append(PE,pe) # to plot
    f = np.linalg.norm(a.get_forces()[0]) #magnitude of force
    F = np.append(F,f) # to plot

minpe = np.argmin(PE) # find index for min value
print('\ndistance for minimum PE is: ' + str(dist[:,2][minpe]) + ' Ang\n')

plt.grid()
plt.plot(oneton,PE)
plt.xlabel('Distance / Ang')
plt.ylabel('PE / eV')
plt.title('Distance vs PE')
plt.figure(2)
plt.grid()
plt.plot(oneton,F,'-r')
plt.xlabel('Distance / Ang')
plt.ylabel('F / eV/Ang')
plt.title('Distance vs Force')
plt.show()