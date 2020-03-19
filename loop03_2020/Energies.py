from numba import jit
import time
path="D:\Babich\Science\Ostapchuk\LOOP_New_Try_05_2018\Python\\"
import sys
sys.path.append(path)
# from closest import closest

#from Load_Energies import k, ee, ek, dee, dek


'''m= (R - r) * (R - r) + r * r
p = ((R + r) * (R + r) + z * z)  0.5
t = R * R - r * r - z * z

k_2 = 4 * R * r / p ^ 2'''
def closest(array,value):
    '''Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively.'''
    n = len(array)
    # if (value < array[0]):
    #     return -1
    # elif (value > array[n-1]):
    #     return n
    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return array[jl]

# def closest(L,x):
#     dist=abs(L[0]-x)
#     closest=L[0]
#     for el in L:
#         if dist>abs(el-x):
#             dist=abs(el-x)
#             closest=el
#     return closest#, dist
def e(r,z,R,k,ee,ek,P):
    e_start=time.time()
    m = ((R - r) * (R - r) + z * z)
    p = ((R + r) * (R + r) + z * z)** 0.5
    t = (R * R - r * r - z * z)

    #print ('m=',type(m),'p=',p,'t=',t)
    # print '??=',(4. * R * r / ((R + r) * (R + r) + z * z)),'r=',r,'R=','z=',z
    _k = (4. * R * r / ((R + r) * (R + r) + z * z))**0.5


    k=closest(k,_k)
    # print (f'_k={_k}')
    # print (f'k={k}')

    _ee=ee[k]


    _ek = ek[k]

    _e=(t*_ee/m +_ek)/p
    # print(f'Time to calk E={time.time()-e_start}')
    return -_e*P
def er(r,z,R,k,ee,ek,dee,dek, P):
    m = (R - r) * (R - r) + z * z
    p = ((R + r) * (R + r) + z * z) ** 0.5
    t = R * R - r * r - z * z
    _k = (4 * R * r / p ** 2)**0.5
    k2r=4*R*(1-2*r*(R+r)/p**2)/p**2

    k = closest(k, _k)
    _ee = ee[k]
    _ek = ek[k]
    _dee = dee[k]
    _dek=dek[k]
    _er=-(R+r)*(t*_ee/m +_ek)/p**3-(2*r*m-2*(R-r)*t)*_ee/(p*m**2)+k2r*(t*_dee/m+_dek)/p


    return -_er*P


def ez(r,z,R,k,ee,ek,dee,dek, P):
    m = (R - r) * (R - r) + z * z
    p = ((R + r) * (R + r) + z * z) ** 0.5
    t = R * R - r * r - z * z
    _k = (4 * R * r / p ** 2)**0.5
    k2z=-8*z*R*r/p**4

    k = closest(k, _k)
    _ee = ee[k]
    _ek = ek[k]
    _dee = dee[k]
    _dek=dek[k]
    _ez=-z*(t*_ee/m+_ek)/p**3-4*z*(R-r)*R*_ee/(p*m**2)+k2z*(t*_dee/m+_dek)/p


    return -_ez*P
