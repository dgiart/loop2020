from numba import jit

path="D:\Babich\Science\Ostapchuk\LOOP_New_Try_05_2018\Python\\"
import sys
sys.path.append(path)
# from closest import closest

#from Load_Energies import k, ee, ek, dee, dek


'''m= (R - r) * (R - r) + r * r
p = ((R + r) * (R + r) + z * z)  0.5
t = R * R - r * r - z * z

k_2 = 4 * R * r / p ^ 2'''


def closest(L,x):
    dist=abs(L[0]-x)
    closest=L[0]
    for el in L:
        if dist>abs(el-x):
            dist=abs(el-x)
            closest=el
    return closest#, dist
def e(r,z,R,k,ee,ek,P):
    m = ((R - r) * (R - r) + z * z)
    p = ((R + r) * (R + r) + z * z)** 0.5
    t = (R * R - r * r - z * z)

    #print ('m=',type(m),'p=',p,'t=',t)
    # print '??=',(4. * R * r / ((R + r) * (R + r) + z * z)),'r=',r,'R=','z=',z
    _k = (4. * R * r / ((R + r) * (R + r) + z * z))**0.5


    k=closest(k,_k)

    _ee=ee[k]


    _ek = ek[k]

    _e=(t*_ee/m +_ek)/p

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
