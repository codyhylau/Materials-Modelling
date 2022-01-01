from __future__ import print_function
import numpy as np
from ase import Atoms
from ase.units import eV, Ang, GPa
import Morse
from matplotlib import pyplot as plt
from ase.build import bulk
calc = Morse.MorsePotential()

# Find Poisson's Ration v (v=0.33)
xstrain = 0.01
cube = 3.6021
cu = bulk("Cu", "fcc", a=cube, cubic=True)
cu.set_calculator(calc)
cell=cu.get_cell()
cell[0][0] *=  1-xstrain # apply strain only in x direction
cu.set_cell(cell, scale_atoms=True)

s0 = cu.get_stress(voigt=False)[1,1]

yzstrain = 0.5
bi_interval = [-yzstrain,yzstrain] #initial bisection interval (start large)
def bisection(interval): # Bisection method to narrow in the yzstrain which leads to 0 stress
    cell[1][1] = cube*(1-interval[0])
    cell[2][2] = cube*(1-interval[0])
    cu.set_cell(cell, scale_atoms=True)
    s1 = -s0-cu.get_stress(voigt=False)[1,1]

    cell[1][1] = cube*(1-interval[1])
    cell[2][2] = cube*(1-interval[1])
    cu.set_cell(cell, scale_atoms=True)
    s2 = -s0-cu.get_stress(voigt=False)[1,1]
    
    mid = sum(interval)/2
    cell[1][1] = cube*(1-mid)
    cell[2][2] = cube*(1-mid)
    cu.set_cell(cell, scale_atoms=True)
    s3 = -s0-cu.get_stress(voigt=False)[1,1]

    # check which one of the interval to be disgarded
    if s3 * s1 < 0: 
        return sorted([mid,interval[0]]),[s1,s2,s3]
    elif s3 * s2 < 0:
        return sorted([mid,interval[1]]),[s1,s2,s3]
    elif s1*s2 > 0: # BAD
        return interval,[s1,s2,s3]

n=0
while n <= 250:

    bi_interval , yzstress = bisection(bi_interval) # performs bisection to find root

    if yzstress[0]*yzstress[1] > 0: # BAD
        print(bi_interval)
        print(yzstress)
        print('\n error: bisection interval doesn''t include the root!\n')
        break
    if abs(yzstress[0]) < 1e-12:
        v = -(sum(bi_interval)/2) / xstrain # calculation of poisson's ration
        print('\nPoisson''s Ratio: ',v,'\n')
        break
    print(bi_interval)
    print(yzstress)
    n += 1


# test using v from above
xstrain = 0.01
cube = 3.6021
cu = bulk("Cu", "fcc", a=cube, cubic=True)
cu.set_calculator(calc)
cell=cu.get_cell()
cell[0,0] *= 1 - xstrain # apply strain only in x direction
cell[1,1] *= 1 + xstrain*v
cell[2,2] *= 1 + xstrain*v
print('Cell sizes: ',cell)
cu.set_cell(cell, scale_atoms=True)
print('Stress Matrix: ',cu.get_stress(voigt=False))

# so why does my value for v lead to essentially 0 y,z pressure but v is too large?