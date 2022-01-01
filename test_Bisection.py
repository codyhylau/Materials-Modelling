import numpy as np

def f(x):
    return x**3 - 27  # root at 3

def bisection(interval): # Bisection method to narrow in the yzstrain which leads to 0 stress

    s1 = f(interval[0])
    s2 = f(interval[1])
    mid = sum(interval)/2
    s3 = f(mid)

    if s3 * s1 < 0:
        return [mid,interval[0]],[s1,s2,s3]
    elif s3 * s2 < 0:
        return [mid,interval[1]],[s1,s2,s3]
    elif s1*s2 > 0:
        return [interval],[s1,s2,s3]


n=0
bi_interval = [-100,100]
while n <= 250:
    bi_interval , yzstress = bisection(bi_interval) # perform bisection to find root
    if yzstress[0]*yzstress[1] > 0:
        print(bi_interval)
        print('\n error: bisection interval doesn''t include the root!\n')
        break
    if abs(yzstress[0]) < 1e-10:
        print('\nPoisson''s Ration: ',bi_interval[0],'\n')
        break
    print(bi_interval)
    print(yzstress)
    n += 1